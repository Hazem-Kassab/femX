import numpy as np

import csv

from strucutre import Structure
from time import process_time

from stress_measures import tresca_stress, von_mises_stress


class Solver:

    def __init__(self, structure: Structure):
        self.structure = structure
        self.dofs = list(self.structure.free_degrees_of_freedom) + list(
            self.structure.restrained_degrees_of_freedom)

    def global_stiffness_matrix(self):
        dofs_count = len(self.dofs)
        matrix = np.zeros((dofs_count, dofs_count))
        for element in self.structure.elements:
            for i, dof_i in enumerate(element.degrees_of_freedom):
                dof_i_index = self.dofs.index(dof_i)
                for j, dof_j in enumerate(element.degrees_of_freedom):
                    dof_j_index = self.dofs.index(dof_j)
                    matrix[dof_i_index, dof_j_index] += element.stiffness_matrix()[i, j]
        return matrix

    def ff_matrix(self):
        fdofs = len(self.structure.free_degrees_of_freedom)
        return self.global_stiffness_matrix()[:fdofs, :fdofs]

    def fs_matrix(self):
        fdofs = len(self.structure.free_degrees_of_freedom)
        return self.global_stiffness_matrix()[:fdofs, fdofs:]

    def sf_matrix(self):
        return self.fs_matrix().T

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
        t1 = process_time()
        fdofs = len(self.structure.free_degrees_of_freedom)
        displacement_vector = np.linalg.inv(self.ff_matrix()).dot(self.force_vector()[:fdofs])
        reactions_vector = self.sf_matrix().dot(displacement_vector)
        i = 0
        for dof in self.structure.free_degrees_of_freedom:
            dof.displacement = displacement_vector[i]
            i += 1
        self.assign_results_to_points()
        t2 = process_time()
        print(f"solution time: {t2 - t1} seconds")

    def average_nodal_values(self):
        for element in self.structure.elements:
            for i, node in enumerate(element.nodes):
                zeta, eta = element.nodal_mapped_coordinates.T[i]
                node.averaged_stress = np.vstack([node.averaged_stress, element.stress(zeta, eta)])
                node.averaged_strain = np.vstack([node.averaged_strain, element.interpolated_strain(zeta, eta)])

    def output(self):
        # self.average_nodal_values()
        f0 = open(r"output/x_coordinates.txt", 'w')
        f1 = open(r"output/y_coordinates.txt", 'w')
        f2 = open(r"output/displacement_x.txt", 'w')
        f3 = open(r"output/displacement_y.txt", 'w')
        f4 = open(r"output/stress_x.txt", 'w')
        f5 = open(r"output/stress_y.txt", 'w')
        f6 = open(r"output/stress_xy.txt", 'w')
        f7 = open(r"output/strain_x.txt", 'w')
        f8 = open(r"output/strain_y.txt", 'w')
        f9 = open(r"output/strain_xy.txt", 'w')
        f10 = open(r"output/p_stress_1.txt", 'w')
        f11 = open(r"output/p_stress_2.txt", 'w')
        f12 = open(r"output/von_mises_stress.txt", 'w')
        f13 = open(r"output/element_coordinates.csv", 'w')
        f14 = open(r"output/element_displaced_coordinates.csv", 'w')
        f15 = open(r"output/x_displaced_coordinates.txt", 'w')
        f16 = open(r"output/y_displaced_coordinates.txt", 'w')
        f17 = open(r"output/number_of_elements.txt", 'w')
        f18 = open(r"output/averaged_stress_x.txt", 'w')
        f19 = open(r"output/averaged_stress_y.txt", 'w')
        f20 = open(r"output/averaged_stress_xy.txt", 'w')
        f21 = open(r"output/averaged_strain_x.txt", 'w')
        f22 = open(r"output/averaged_strain_y.txt", 'w')
        f23 = open(r"output/averaged_strain_xy.txt", 'w')

        writer_1 = csv.writer(f13)
        writer_2 = csv.writer(f14)
        f17.write(f"{len(self.structure.elements)}")
        f17.close()
        for element in self.structure.elements:
            element.node_coordinates()
            writer_1.writerow(f"{node.x}" for node in element.nodes)
            writer_2.writerow(
                f"{[str(node.coordinates + np.array([node.x_dof.displacement, node.y_dof.displacement])) for node in element.nodes]}\n")
            zetas = etas = np.linspace(-1, 1, 20, True)
            for point in element.points:
                # print(point.stress_tensor.tensor)
                # for eta in etas:
                zeta = point.zeta
                eta = point.eta
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
                f10.write(f"{point.principle_stress_tensor.tensor[0, 0]} ")
                f11.write(f"{point.principle_stress_tensor.tensor[1, 1]} ")
                f12.write(f"{von_mises_stress(point)} ")
                f15.write(f"{x + element.displacement(point)[0]} ")
                f16.write(f"{y + element.displacement(point)[1]} ")
                # f18.write(f"{element.averaged_stress(zeta, eta)[0]} ")
                # f19.write(f"{element.averaged_stress(zeta, eta)[1]} ")
                # f20.write(f"{element.averaged_stress(zeta, eta)[2]} ")
                # f21.write(f"{element.averaged_strain(zeta, eta)[0]} ")
                # f22.write(f"{element.averaged_strain(zeta, eta)[1]} ")
                # f23.write(f"{element.averaged_strain(zeta, eta)[2]} ")

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
