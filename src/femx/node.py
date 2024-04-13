import numpy as np

from degree_of_freedom import DOF


class Node:
    counter = 0

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.x_dof = DOF()
        self.y_dof = DOF()
        self.averaged_stress = np.zeros((1, 3))
        self.averaged_strain = np.zeros((1, 3))
        Node.counter += 1
        self.id = Node.counter

    @property
    def degrees_of_freedom(self):
        return [self.x_dof, self.y_dof]

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x = value

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        self._y = value

    @property
    def coordinates(self):
        return np.array([self.x, self.y])
    
    @property
    def averaged_stress(self) -> np.ndarray:
        return self._averaged_stress

    @averaged_stress.setter
    def averaged_stress(self, value: np.ndarray):
        self._averaged_stress = value

    @property
    def averaged_strain(self) -> np.ndarray:
        return self._averaged_strain

    @averaged_strain.setter
    def averaged_strain(self, value: np.ndarray):
        self._averaged_strain = value

    def __repr__(self):
        return f"Node {self.id} at ({self.x}, {self.y})"
