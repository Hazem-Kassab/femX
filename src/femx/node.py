import numpy as np

from femx.degree_of_freedom import DOF


class Node:
    counter = 0

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.x_dof = DOF()
        self.y_dof = DOF()
        Node.counter += 1
        self.id = Node.counter

    @property
    def degrees_of_freedom(self):
        return [self.x_dof, self.y_dof]

    @property
    def coordinates(self):
        return np.array([self.x, self.y])

    def __repr__(self):
        return f"Node {self.id} at ({self.x}, {self.y})"
