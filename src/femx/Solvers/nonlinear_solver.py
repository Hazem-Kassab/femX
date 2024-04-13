import numpy as np

from Solvers.solver import Solver
from degree_of_freedom import DOF
from stress_measures import von_mises_stress


class NonlinearSolver(Solver):

    def solve_incrementally(self, dof_recorded_force: DOF, dof_recorded_displacement: DOF, preforce, steps=10):
        force = [0, abs(preforce)]
        displacement = [0, abs(dof_recorded_displacement.displacement)]
        incremental_force_vector = np.array([dof.force for dof in self.structure.free_degrees_of_freedom]) / steps
        dof_index = self.structure.free_degrees_of_freedom.index(dof_recorded_force)
        for i in range(steps):
            print(i)
            # print(f"external force = {incremental_force_vector[dof_index] * i + preforce}")
            # print(f"internal_force = {self.internal_force_vector()[dof_index]}")
            force.append(abs(incremental_force_vector[dof_index] * i + preforce))
            displacement.append(abs(dof_recorded_displacement.displacement))
            try:
                incremental_displacement_vector = np.linalg.inv(self.ff_matrix()).dot(incremental_force_vector)
            except np.linalg.LinAlgError:
                print(incremental_force_vector[dof_index] * i)
                print(dof_recorded_displacement.displacement)
                break
            for index, dof in enumerate(self.structure.free_degrees_of_freedom):
                dof.displacement += incremental_displacement_vector[index]
            self.assign_results_to_points()
            self.update_yield_status()
        return displacement, force

    def update_yield_status(self):
        for element in self.structure.elements:
            for point in element.gauss_points:
                element.material.test_yield(von_mises_stress, point)

    def internal_force_vector(self):
        disp = np.array([dof.displacement for dof in self.structure.free_degrees_of_freedom])
        return self.ff_matrix().dot(disp)
