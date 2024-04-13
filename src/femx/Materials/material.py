from abc import ABC, abstractmethod

from point import Point


class Material(ABC):

    def __init__(self, elasticity_modulus, poissons_ratio):
        self.elasticity_modulus = elasticity_modulus
        self.poissons_ratio = poissons_ratio

    @property
    def elasticity_modulus(self):
        return self._elasticity_modulus

    @elasticity_modulus.setter
    def elasticity_modulus(self, value):
        self._elasticity_modulus = value

    @property
    def poissons_ratio(self):
        return self._poissons_ratio

    @poissons_ratio.setter
    def poissons_ratio(self, value):
        self._poissons_ratio = value

    @property
    def shear_modulus(self):
        return self.elasticity_modulus / (2 * (1 + self.poissons_ratio))

    @abstractmethod
    def test_yield(self, func, point: Point):
        pass
