import numpy as np

from femx.point import Point


def von_mises_stress(point: Point):
    J2 = point.stress_deviator_tensor.second_invariant
    return (3*abs(J2))**0.5


def tresca_stress(point: Point):
    zeta, rho, theta = point.westergaard_coordinates
    return 2 ** 0.5 * rho * np.sin(theta + np.pi / 3)