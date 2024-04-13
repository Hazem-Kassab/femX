import numpy as np


class Tensor:

    def __init__(self, tensor: np.ndarray):
        self.tensor = tensor

    @property
    def tensor(self) -> np.ndarray:
        return self._tensor

    @tensor.setter
    def tensor(self, value: np.ndarray):
        self._tensor = value

    def _eigenpairs(self):
        return np.linalg.eig(self.tensor)

    @property
    def principle_directions(self) -> np.ndarray:
        return self._eigenpairs()[1]

    @property
    def diagonalized(self) -> np.ndarray:
        return self.principle_directions.T.dot(self.tensor.dot(self.principle_directions))

    @property
    def first_invariant(self):
        return self.tensor.trace()

    @property
    def second_invariant(self):
        s11 = self.diagonalized[0, 0]
        s22 = self.diagonalized[1, 1]
        s33 = self.diagonalized[2, 2]
        return s11 * s22 + s22 * s33 + s33 * s11

    @property
    def third_invariant(self):
        s11 = self.diagonalized[0, 0]
        s22 = self.diagonalized[1, 1]
        s33 = self.diagonalized[2, 2]
        return s11 * s22 * s33
