import numpy as np

from femx.point import Point


class GaussPoint(Point):

    def __init__(self, zeta, eta, weight):
        super().__init__(zeta, eta)
        self.weight = weight
        self.yielded = False

    @property
    def yielded(self):
        return self._yielded

    @yielded.setter
    def yielded(self, value: bool):
        if value is True:
            self.weight = 0
        self._yielded = value
