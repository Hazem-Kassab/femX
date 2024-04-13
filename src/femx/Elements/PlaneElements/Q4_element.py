import math

import numpy as np

from gauss_point import GaussPoint
from node import Node
from Materials.material import Material
from Elements.element import Element
from point import Point


class Q4Element(Element):

    def __init__(self, node_1: Node, node_2: Node, node_3: Node, node_4: Node, thickness, material: Material,
                 plane_stress=True):
        super().__init__(thickness, material)
        self.node_1 = node_1
        self.node_2 = node_2
        self.node_3 = node_3
        self.node_4 = node_4
        self.gauss_point_1 = GaussPoint(-1/math.sqrt(3), -1/math.sqrt(3), 1)
        self.gauss_point_2 = GaussPoint(1/math.sqrt(3), -1/math.sqrt(3), 1)
        self.gauss_point_3 = GaussPoint(1 / math.sqrt(3), 1 / math.sqrt(3), 1)
        self.gauss_point_4 = GaussPoint(-1 / math.sqrt(3), 1 / math.sqrt(3), 1)
        self.point_1 = Point(-1, -1)
        self.point_2 = Point(1, -1)
        self.point_3 = Point(1, 1)
        self.point_4 = Point(-1, 1)

    @property
    def nodal_points(self):
        return [self.point_1, self.point_2, self.point_3, self.point_4]

    @property
    def gauss_points(self):
        return [self.gauss_point_1, self.gauss_point_2, self.gauss_point_3, self.gauss_point_4]

    @property
    def nodes(self) -> list[Node]:
        return [self.node_1, self.node_2, self.node_3, self.node_4]

    def shape_function(self, point: Point):
        zeta, eta = point.coordinates
        n1 = 0.25 * (1 - zeta) * (1 - eta)
        n2 = 0.25 * (1 + zeta) * (1 - eta)
        n3 = 0.25 * (1 + zeta) * (1 + eta)
        n4 = 0.25 * (1 - zeta) * (1 + eta)
        return np.array([n1, n2, n3, n4])

    def jacobian(self, gauss_point: GaussPoint):
        zeta, eta= gauss_point.coordinates
        derivative = 0.25 * np.array([[-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
                                      [-(1 - zeta), -(1 + zeta), (1 + zeta), (1 - zeta)]])

        return derivative.dot(self.node_coordinates().T)

    def B_matrix(self, gauss_point: GaussPoint):
        zeta, eta = gauss_point.coordinates
        m1 = np.array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 1, 0]])

        m2 = np.zeros((4, 4))

        inv_jacobian = np.linalg.inv(self.jacobian(gauss_point))
        m2[:2, :2] = inv_jacobian
        m2[2:, 2:] = inv_jacobian

        ne1 = -0.25 * (1 - zeta)
        ne2 = -0.25 * (1 + zeta)
        ne3 = 0.25 * (1 + zeta)
        ne4 = 0.25 * (1 - zeta)

        nz1 = -0.25 * (1 - eta)
        nz2 = 0.25 * (1 - eta)
        nz3 = 0.25 * (1 + eta)
        nz4 = -0.25 * (1 + eta)

        m3 = np.array([[nz1, 0, nz2, 0, nz3, 0, nz4, 0],
                       [ne1, 0, ne2, 0, ne3, 0, ne4, 0],
                       [0, nz1, 0, nz2, 0, nz3, 0, nz4],
                       [0, ne1, 0, ne2, 0, ne3, 0, ne4]])

        return m1.dot(m2.dot(m3))
