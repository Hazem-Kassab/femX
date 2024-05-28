import time


from pathlib import Path

import numpy as np

import csv

from femx.strucutre import Structure
from femx.elements.element import Element

from femx.stress_measures import tresca_stress, von_mises_stress


class Solver:

    def __init__(self, structure: Structure):
        self.structure = structure
        self.dofs = list(self.structure.free_degrees_of_freedom) + list(
            self.structure.restrained_degrees_of_freedom)

    def global_stiffness_matrix(self):
        dofs_count = len(self.dofs)
        matrix = np.zeros((dofs_count, dofs_count))
        self._assemble_global_stiffness_matrix(matrix, self.structure.elements)
        return matrix + matrix.T - np.diag(np.diag(matrix))

    def _assemble_global_stiffness_matrix(self, matrix, elements: list[Element]):
        for element in elements:
            element_stiffness_matrix = element.stiffness_matrix()
            i = 0
            for dof_i in element.degrees_of_freedom:
                j = i
                dof_i_index = self.dofs.index(dof_i)
                for dof_j in element.degrees_of_freedom[i:]:
                    dof_j_index = self.dofs.index(dof_j)
                    matrix[dof_i_index, dof_j_index] += element_stiffness_matrix[i, j]
                    j += 1
                i += 1

    def ff_matrix(self):
        fdofs = len(self.structure.free_degrees_of_freedom)
        return self.global_stiffness_matrix()[:fdofs, :fdofs]

    def sf_matrix(self):
        fdofs = len(self.structure.free_degrees_of_freedom)
        return self.global_stiffness_matrix()[fdofs:, :fdofs]

    def fs_matrix(self):
        return self.sf_matrix().T

    def ss_matrix(self):
        fdofs = len(self.structure.free_degrees_of_freedom)
        return self.global_stiffness_matrix()[fdofs:, fdofs:]

    def force_vector(self):
        return np.array([dof.force for dof in self.dofs])

    def assign_results_to_points(self):
        for element in self.structure.elements:
            for point in element.points:
                sxx = element.stress(point)[0]
                syy = element.stress(point)[1]
                sxy = element.stress(point)[2]

                exx = element.strain(point)[0]
                eyy = element.strain(point)[1]
                exy = element.strain(point)[2]

                point.stress_tensor.tensor = np.array([[sxx, sxy, 0],
                                                       [sxy, syy, 0],
                                                       [0, 0, 0]])

                point.strain_tensor.tensor = np.array([[exx, exy, 0],
                                                       [exy, eyy, 0],
                                                       [0, 0, 0]])

    def solve(self):
        print("Solving...")
        t1 = time.time()
        fdofs = len(self.structure.free_degrees_of_freedom)
        print("Assembling global stiffness matrix...")
        t1b = time.time()
        gmatrix = self.ff_matrix()
        t2b = time.time()
        print(f"Assembly done in {t2b - t1b:.3f} seconds")
        print("*"*30)
        print("inverting stiffness matrix...")
        t1a = time.time()
        inverse_matrix = np.linalg.inv(gmatrix)
        t2a = time.time()
        print(f"Inversion done in {t2a - t1a:.3f} seconds")
        print("*"*30)
        print("Solving for displacements...")
        t1c = time.time()
        displacement_vector = inverse_matrix.dot(self.force_vector()[:fdofs])
        t2c = time.time()
        print(f"Solved in {t2c - t1c:.3f} seconds")
        print("*"*30)
        # reactions_vector = self.sf_matrix().dot(displacement_vector)
        i = 0
        for dof in self.structure.free_degrees_of_freedom:
            dof.displacement = displacement_vector[i]
            i += 1
        self.assign_results_to_points()
        t2 = time.time()
        print(f"Solution time: {t2 - t1:.3f} seconds")

    def average_nodal_values(self):
        for element in self.structure.elements:
            for i, node in enumerate(element.nodes):
                zeta, eta = element.nodal_mapped_coordinates.T[i]
                node.averaged_stress = np.vstack([node.averaged_stress, element.stress(zeta, eta)])
                node.averaged_strain = np.vstack([node.averaged_strain, element.interpolated_strain(zeta, eta)])

    def output(self):
        print("Generating output files...")
        path = Path(__file__).parent
        f0 = open((path / "../output/x_coordinates.txt").resolve(), 'w')
        f1 = open((path / "../output/y_coordinates.txt").resolve(), 'w')
        f2 = open((path / "../output/displacement_x.txt").resolve(), 'w')
        f3 = open((path / "../output/displacement_y.txt").resolve(), 'w')
        f4 = open((path / "../output/stress_x.txt").resolve(), 'w')
        f5 = open((path / "../output/stress_y.txt").resolve(), 'w')
        f6 = open((path / "../output/stress_xy.txt").resolve(), 'w')
        f7 = open((path / "../output/strain_x.txt").resolve(), 'w')
        f8 = open((path / "../output/strain_y.txt").resolve(), 'w')
        f9 = open((path / "../output/strain_xy.txt").resolve(), 'w')
        f10 = open((path / "../output/p_stress_1.txt").resolve(), 'w')
        f11 = open((path / "../output/p_stress_2.txt").resolve(), 'w')
        f12 = open((path / "../output/von_mises_stress.txt").resolve(), 'w')
        f13 = open((path / "../output/element_coordinates.txt").resolve(), 'w')
        f14 = open((path / "../output/element_displaced_coordinates.txt").resolve(), 'w')
        f15 = open((path / "../output/x_displaced_coordinates.txt").resolve(), 'w')
        f16 = open((path / "../output/y_displaced_coordinates.txt").resolve(), 'w')
        f17 = open((path / "../output/number_of_elements.txt").resolve(), 'w')
        f18 = open((path / "../output/averaged_stress_x.txt").resolve(), 'w')
        f19 = open((path / "../output/averaged_stress_y.txt").resolve(), 'w')
        f20 = open((path / "../output/averaged_stress_xy.txt").resolve(), 'w')
        f21 = open((path / "../output/averaged_strain_x.txt").resolve(), 'w')
        f22 = open((path / "../output/averaged_strain_y.txt").resolve(), 'w')
        f23 = open((path / "../output/averaged_strain_xy.txt").resolve(), 'w')

        writer_1 = csv.writer(f13)
        writer_2 = csv.writer(f14)
        f17.write(f"{len(self.structure.elements)}")
        f17.close()
        for element in self.structure.elements:
            element.node_coordinates()
            writer_1.writerow(f"{node.x}" for node in element.nodes)
            writer_2.writerow(
                f"{[str(node.coordinates + np.array([node.x_dof.displacement, node.y_dof.displacement])) for node in element.nodes]}\n")
            for point in element.points:
                x, y = element.get_coordinates(point)
                f0.write(f"{x} ")
                f1.write(f"{y} ")
                f2.write(f"{element.displacement(point)[0]} ")
                f3.write(f"{element.displacement(point)[1]} ")
                f4.write(f"{point.stress_tensor.tensor[0, 0]} ")
                f5.write(f"{point.stress_tensor.tensor[1, 1]} ")
                f6.write(f"{point.stress_tensor.tensor[0, 1]} ")
                f7.write(f"{point.strain_tensor.tensor[0, 0]} ")
                f8.write(f"{point.strain_tensor.tensor[1, 1]} ")
                f9.write(f"{point.strain_tensor.tensor[0, 1]} ")
                f10.write(f"{point.principal_stress_tensor.tensor[0, 0]} ")
                f11.write(f"{point.principal_stress_tensor.tensor[1, 1]} ")
                f12.write(f"{von_mises_stress(point)} ")
                f15.write(f"{x + element.displacement(point)[0]} ")
                f16.write(f"{y + element.displacement(point)[1]} ")

        f0.close()
        f1.close()
        f2.close()
        f3.close()
        f4.close()
        f5.close()
        f6.close()
        f7.close()
        f8.close()
        f9.close()
        f10.close()
        f11.close()
        f12.close()
        f13.close()
        f15.close()
        f16.close()
        f18.close()
        f19.close()
        f20.close()
        f21.close()
        f22.close()
        f23.close()
        print("Output files generated.")
