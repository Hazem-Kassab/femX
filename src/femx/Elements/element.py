from abc import ABC, abstractmethod

import numpy as np

from Materials.material import Material
from degree_of_freedom import DOF
from gauss_point import GaussPoint
from node import Node
from point import Point


class Element(ABC):
    counter = 0

    def __init__(self, thickness, material: Material):
        self.thickness = thickness
        self.material = material
        Element.counter += 1
        self.id = Element.counter

    @property
    def thickness(self):
        return self._thickness

    @thickness.setter
    def thickness(self, value):
        self._thickness = value

    @property
    def material(self) -> Material:
        return self._material

    @material.setter
    def material(self, value: Material):
        self._material = value

    @abstractmethod
    def jacobian(self, point: GaussPoint):
        raise NotImplementedError

    # @property
    # def area(self):
    #     area = 0
    #     for zeta in self.sampling_points:
    #         for eta in self.sampling_points:
    #             area += np.linalg.det(self.jacobian(zeta, eta))
    #     return area

    @property
    @abstractmethod
    def nodes(self) -> list[Node]:
        raise NotImplementedError

    @property
    @abstractmethod
    def nodal_points(self):
        raise NotImplementedError
    @property
    def degrees_of_freedom(self) -> list[DOF]:
        dofs = []
        for node in self.nodes:
            for dof in node.degrees_of_freedom:
                dofs.append(dof)
        return dofs

    @property
    def points(self) -> list[Point]:
        return self.nodal_points + self.gauss_points

    @property
    def nodal_mapped_coordinates(self):
        return np.array([point.coordinates for point in self.nodal_points]).T

    @abstractmethod
    def shape_function(self, point: Point):
        raise NotImplementedError

    @abstractmethod
    def B_matrix(self, gauss_point: Point):
        raise NotImplementedError

    def stiffness_matrix(self):
        n = len(self.degrees_of_freedom)
        matrix = np.zeros((n, n))
        for gp in self.gauss_points:
            matrix += gp.weight * self.B_matrix(gp).T.dot(
                self.elasticity_matrix().dot(self.B_matrix(gp))) * self.thickness * np.linalg.det(
                self.jacobian(gp))

        return matrix

    def node_coordinates(self):
        return np.array([node.coordinates for node in self.nodes]).T

    def nodal_displacement_vector(self):
        return np.array([dof.displacement for dof in self.degrees_of_freedom])

    def elasticity_matrix(self):
        e = self.material.elasticity_modulus
        v = self.material.poissons_ratio
        return e / (1 - v ** 2) * np.array([[1, v, 0],
                                            [v, 1, 0],
                                            [0, 0, (1 - v) / 2]])

    # @abstractmethod
    def strain_shape_function(self, zeta, eta):
        raise NotImplementedError

    def strain(self, point: Point):
        return self.B_matrix(point).dot(self.nodal_displacement_vector())

    def stress(self, point: Point):
        return self.elasticity_matrix().dot(self.strain(point))

    @property
    @abstractmethod
    def gauss_points(self) -> list[GaussPoint]:
        raise NotImplementedError

    def gauss_points_strain(self):
        gauss_points_strain = []
        for gp in self.gauss_points:
            strain = self.strain(gp.zeta, gp.eta)
            gauss_points_strain.append(strain)

        return np.array(gauss_points_strain).T

    def interpolated_strain(self, zeta, eta):
        return self.gauss_points_strain().dot(self.strain_shape_function(zeta, eta))

    def averaged_nodal_stress(self):
        return np.array([np.average(node.averaged_stress, axis=0) for node in self.nodes]).T

    def averaged_stress(self, zeta, eta):
        return self.averaged_nodal_stress().dot(self.shape_function(zeta, eta))

    def averaged_nodal_strain(self):
        return np.array([np.average(node.averaged_strain, axis=0) for node in self.nodes]).T

    def averaged_strain(self, zeta, eta):
        return self.averaged_nodal_strain().dot(self.shape_function(zeta, eta))

    def stress_tensor(self, zeta, eta):
        sxx, syy, txy = self.stress(zeta, eta)
        return np.array([[sxx, txy],
                         [txy, syy]])

    def principle_stress(self, zeta, eta):
        eigenvalues, eigenvectors = np.linalg.eig(self.stress_tensor(zeta, eta))
        return eigenvalues

    def von_mises_stress(self, zeta, eta):
        s1, s2 = self.principle_stress(zeta, eta)
        return (s1**2 - s1*s2 + s2**2)**0.5

    def get_coordinates(self, point: Point):
        return self.node_coordinates().dot(self.shape_function(point))

    def displacement(self, point: Point):
        nodal_displacements = self.nodal_displacement_vector().reshape(self.node_coordinates().T.shape).T
        return nodal_displacements.dot(self.shape_function(point))

    def reset(self):
        pass

    def __repr__(self):
        s = f"#### Element {self.id} ####"
        for node in self.nodes:
            s += f"\n node {node.id} at {node.coordinates}"
        return s
