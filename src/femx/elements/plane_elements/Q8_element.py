# 8-node plane-stress element with 9 Gauss Points

import numpy as np

from femx.elements.element import Element
from femx.materials.material import Material
from femx.gauss_point import GaussPoint
from femx.node import Node
from femx.point import Point


class Q8Element(Element):

    def __init__(self, node_1: Node, node_2: Node, node_3: Node, node_4: Node, node_5: Node, node_6: Node, node_7: Node,
                 node_8: Node, thickness, material: Material):
        super().__init__(thickness, material)
        self.node_1 = node_1
        self.node_2 = node_2
        self.node_3 = node_3
        self.node_4 = node_4
        self.node_5 = node_5
        self.node_6 = node_6
        self.node_7 = node_7
        self.node_8 = node_8
        self.gauss_point_1 = GaussPoint(-0.6**0.5, -0.6**0.5, 5/9)
        self.gauss_point_2 = GaussPoint(0, -0.6**0.5, 5/9)
        self.gauss_point_3 = GaussPoint(0.6**0.5, -0.6**0.5, 5/9)
        self.gauss_point_4 = GaussPoint(-0.6**0.5, 0, 5/9)
        self.gauss_point_5 = GaussPoint(0, 0, 8/9)
        self.gauss_point_6 = GaussPoint(0.6**0.5, 0, 5/9)
        self.gauss_point_7 = GaussPoint(-0.6**0.5, 0.6**0.5, 5/9)
        self.gauss_point_8 = GaussPoint(0, 0.6**0.5, 5/9)
        self.gauss_point_9 = GaussPoint(0.6**0.5, 0.6**0.5, 5/9)
        self.point_1 = Point(-1, 1)
        self.point_2 = Point(0, -1)
        self.point_3 = Point(1, -1)
        self.point_4 = Point(1, 0)
        self.point_5 = Point(1, 1)
        self.point_6 = Point(0, 1)
        self.point_7 = Point(-1, -1)
        self.point_8 = Point(-1, 0)

    @property
    def gauss_points(self):
        return [self.gauss_point_1, self.gauss_point_2, self.gauss_point_3, self.gauss_point_4,
                self.gauss_point_5, self.gauss_point_6, self.gauss_point_7, self.gauss_point_8, self.gauss_point_9]

    @property
    def nodal_points(self):
        return [self.point_1, self.point_2, self.point_3, self.point_4,
                self.point_5, self.point_6, self.point_7, self.point_8]

    def jacobian(self, gauss_point: GaussPoint):
        zeta, eta = gauss_point.coordinates

        n2z = 0.5 * (-2 * zeta) * (1 - eta)
        n2e = -0.5 * (1 - zeta ** 2)

        n4z = 0.5 * (1 - eta ** 2)
        n4e = 0.5 * (1 + zeta) * (-2 * eta)

        n6z = 0.5 * (-2 * zeta) * (1 + eta)
        n6e = 0.5 * (1 - zeta ** 2)

        n8z = -0.5 * (1 - eta ** 2)
        n8e = 0.5 * (1 - zeta) * (-2 * eta)

        n1z = -0.25 * (1 - eta) - 0.5 * (n8z + n2z)
        n1e = -0.25 * (1 - zeta) - 0.5 * (n8e + n2e)

        n3z = 0.25 * (1 - eta) - 0.5 * (n2z + n4z)
        n3e = -0.25 * (1 + zeta) - 0.5 * (n2e + n4e)

        n5z = 0.25 * (1 + eta) - 0.5 * (n4z + n6z)
        n5e = 0.25 * (1 + zeta) - 0.5 * (n4e + n6e)

        n7z = -0.25 * (1 + eta) - 0.5 * (n6z + n8z)
        n7e = 0.25 * (1 - zeta) - 0.5 * (n6e + n8e)

        derivative = np.array([[n1z, n2z, n3z, n4z, n5z, n6z, n7z, n8z],
                               [n1e, n2e, n3e, n4e, n5e, n6e, n7e, n8e]])
        return derivative.dot(self.node_coordinates().T)

    @property
    def nodes(self) -> list[Node]:
        return [self.node_1, self.node_2, self.node_3, self.node_4, self.node_5, self.node_6, self.node_7, self.node_8]

    def shape_function(self, gauss_point: GaussPoint):
        zeta, eta = gauss_point.coordinates
        n2 = 0.5 * (1 - zeta ** 2) * (1 - eta)
        n4 = 0.5 * (1 + zeta) * (1 - eta ** 2)
        n6 = 0.5 * (1 - zeta ** 2) * (1 + eta)
        n8 = 0.5 * (1 - zeta) * (1 - eta ** 2)
        n1 = 0.25 * (1 - zeta) * (1 - eta) - 0.5 * (n8 + n2)
        n3 = 0.25 * (1 + zeta) * (1 - eta) - 0.5 * (n2 + n4)
        n5 = 0.25 * (1 + zeta) * (1 + eta) - 0.5 * (n4 + n6)
        n7 = 0.25 * (1 - zeta) * (1 + eta) - 0.5 * (n6 + n8)

        return np.array([n1, n2, n3, n4, n5, n6, n7, n8])

    def B_matrix(self, gauss_point: GaussPoint):
        zeta, eta = gauss_point.coordinates

        m1 = np.array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 1, 0]])

        m2 = np.zeros((4, 4))

        inv_jacobian = np.linalg.inv(self.jacobian(gauss_point))
        m2[:2, :2] = inv_jacobian
        m2[2:, 2:] = inv_jacobian

        n2z = 0.5 * (-2 * zeta) * (1 - eta)
        n2e = -0.5 * (1 - zeta ** 2)

        n4z = 0.5 * (1 - eta ** 2)
        n4e = 0.5 * (1 + zeta) * (-2 * eta)

        n6z = 0.5 * (-2 * zeta) * (1 + eta)
        n6e = 0.5 * (1 - zeta ** 2)

        n8z = -0.5 * (1 - eta ** 2)
        n8e = 0.5 * (1 - zeta) * (-2 * eta)

        n1z = -0.25 * (1 - eta) - 0.5 * (n8z + n2z)
        n1e = -0.25 * (1 - zeta) - 0.5 * (n8e + n2e)

        n3z = 0.25 * (1 - eta) - 0.5 * (n2z + n4z)
        n3e = -0.25 * (1 + zeta) - 0.5 * (n2e + n4e)

        n5z = 0.25 * (1 + eta) - 0.5 * (n4z + n6z)
        n5e = 0.25 * (1 + zeta) - 0.5 * (n4e + n6e)

        n7z = -0.25 * (1 + eta) - 0.5 * (n6z + n8z)
        n7e = 0.25 * (1 - zeta) - 0.5 * (n6e + n8e)

        m3 = np.array([[n1z, 0, n2z, 0, n3z, 0, n4z, 0, n5z, 0, n6z, 0, n7z, 0, n8z, 0],
                       [n1e, 0, n2e, 0, n3e, 0, n4e, 0, n5e, 0, n6e, 0, n7e, 0, n8e, 0],
                       [0, n1z, 0, n2z, 0, n3z, 0, n4z, 0, n5z, 0, n6z, 0, n7z, 0, n8z],
                       [0, n1e, 0, n2e, 0, n3e, 0, n4e, 0, n5e, 0, n6e, 0, n7e, 0, n8e]])

        return m1.dot(m2.dot(m3))

    def strain_shape_function(self, zeta, eta):
        r = zeta / 3**0.5
        s = eta / 3**0.5
        n1 = 0.25 * (1 - r) * (1 - s)
        n2 = 0.25 * (1 + r) * (1 - s)
        n3 = 0.25 * (1 + r) * (1 + s)
        n4 = 0.25 * (1 - r) * (1 + s)
        return np.array([n1, n2, n3, n4])

