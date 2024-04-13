import numpy as np


class Vector:
    def __init__(self, array: np.ndarray):
        self.vector = array

    @property
    def vector(self) -> np.ndarray:
        return self._vector

    @vector.setter
    def vector(self, array: np.ndarray):
        self._vector = array

    @property
    def magnitude(self):
        return np.sum([x ** 2 for x in self.vector]) ** 0.5

    @property
    def normalized(self) -> np.ndarray:
        return self.vector/self.magnitude
