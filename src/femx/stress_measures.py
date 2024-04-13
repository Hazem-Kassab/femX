import numpy as np

from gauss_point import GaussPoint


def von_mises_stress(point: GaussPoint):
    J2 = point.stress_deviator_tensor.second_invariant
    return (3*abs(J2))**0.5


def tresca_stress(point: GaussPoint):
    zeta, rho, theta = point.westergaard_coordinates
    return 2 ** 0.5 * rho * np.sin(theta + np.pi / 3)