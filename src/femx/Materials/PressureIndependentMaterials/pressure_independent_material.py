from abc import ABC

import numpy as np

from Materials.material import Material
from gauss_point import GaussPoint
from point import Point
from stress_measures import tresca_stress, von_mises_stress


class PressureIndependentMaterial(Material, ABC):

    def __init__(self, elasticity_modulus, poissons_ratio, yield_strength, rupture_strength):
        super().__init__(elasticity_modulus, poissons_ratio)
        self.yield_strength = yield_strength
        self.rupture_strength = rupture_strength

    @property
    def yield_strength(self):
        return self._yield_strength

    @yield_strength.setter
    def yield_strength(self, value):
        self._yield_strength = value

    @property
    def rupture_strength(self):
        return self._rupture_strength

    @rupture_strength.setter
    def rupture_strength(self, value):
        self._rupture_strength = value

    # def test_tresca_yield(self, point: GaussPoint):
    #     so = self.yield_strength
    #     if tresca_stress(point) >= so:
    #         point.yielded = True
    #     # zeta, rho, theta = point.westergaard_coordinates
    #     # so = self.yield_strength
    #     # return 2 ** 0.5 * rho * np.sin(theta + np.pi / 3) - so
    #
    # def test_von_mises_yield(self, point: GaussPoint):
    #     so = self.yield_strength
    #     if von_mises_stress(point) >= so:
    #         point.yielded = True
        # J2 = point.deviator_tensor.second_invariant
        # K = self.yield_strength / 3**0.5
        # return J2 - K**2

    def test_yield(self, func, point: Point):
        so = self.yield_strength
        if func(point) >= so:
            print("Yielded!")
            point.yielded = True


