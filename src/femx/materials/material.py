from abc import ABC, abstractmethod

from femx.point import Point


class Material(ABC):

    def __init__(self, elasticity_modulus, poissons_ratio):
        self.elasticity_modulus = elasticity_modulus
        self.poissons_ratio = poissons_ratio

    @property
    def shear_modulus(self):
        return self.elasticity_modulus / (2 * (1 + self.poissons_ratio))

    @abstractmethod
    def test_yield(self, func, point: Point):
        pass
