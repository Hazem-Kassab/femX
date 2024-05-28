from abc import ABC, abstractmethod

import numpy as np

from femx.materials.material import Material
from femx.point import Point


class FailureCriterion(ABC):

    def __init__(self, material: Material):
        self.material = material

    @staticmethod
    def westergaard_coordinates(point: Point):
        I1 = point.stress_tensor.first_invariant
        J2 = point.stress_deviator_tensor.second_invariant
        J3 = point.stress_deviator_tensor.third_invariant
        zeta = I1 / 3 ** 0.5
        rho = (2 * J2) ** 0.5
        theta = np.arccos(3 * 3 ** 0.5 / 2 * J3 / J2**(3/2))/3
        return np.array([zeta, rho, theta])


    @property
    @abstractmethod
    def function(self, *args):
        raise NotImplementedError
