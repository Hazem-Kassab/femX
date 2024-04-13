import numpy as np

from tensor import Tensor
from vector import Vector


class Point:
    def __init__(self, zeta, eta):
        self.zeta = zeta
        self.eta = eta
        self.coordinates = np.array([zeta, eta])
        self.stress_tensor = Tensor(np.zeros((3, 3)))
        self.strain_tensor = Tensor(np.zeros((3, 3)))

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, array: np.ndarray):
        self._coordinates = array

    @property
    def stress_tensor(self) -> Tensor:
        return self._stress_tensor

    @stress_tensor.setter
    def stress_tensor(self, value: Tensor):
        self._stress_tensor = value

    @property
    def strain_tensor(self) -> Tensor:
        return self._strain_tensor

    @strain_tensor.setter
    def strain_tensor(self, value: Tensor):
        self._strain_tensor = value

    @staticmethod
    def _normalize(vector: np.ndarray) -> np.ndarray:
        return vector/(vector.dot(vector))**0.5

    @staticmethod
    def _magnitude(vector: np.ndarray):
        return np.sum([x**2 for x in vector])**0.5

    @property
    def principle_stress_tensor(self) -> Tensor:
        return Tensor(self.stress_tensor.diagonalized)

    @property
    def principle_strain_tensor(self) -> Tensor:
        return Tensor(self.strain_tensor.diagonalized)

    def stress_vector(self, normal_vector: Vector) -> Vector:
        return Vector(self.stress_tensor.tensor.dot(normal_vector.normalized))

    @property
    def mean_hydrostatic_stress_tensor(self) -> Tensor:
        p = self.stress_tensor.first_invariant / 3
        return Tensor(p*np.identity(3))

    @property
    def mean_hydrostatic_strain_tensor(self):
        sp = self.strain_tensor.first_invariant / 3
        return Tensor(sp*np.identity(3))

    @property
    def stress_deviator_tensor(self) -> Tensor:
        deviator_tensor = self.stress_tensor.tensor - self.mean_hydrostatic_stress_tensor.tensor
        return Tensor(deviator_tensor)

    @property
    def strain_deviator_tensor(self) -> Tensor:
        deviator_tensor = self.strain_tensor.tensor - self.mean_hydrostatic_strain_tensor.tensor
        return Tensor(deviator_tensor)

    @property
    def oct_stress_vector(self) -> Vector:
        unit_normal = Vector(np.array([1, 1, 1])).normalized
        return Vector(self.stress_tensor.diagonalized.dot(unit_normal))

    @property
    def oct_normal_stress_vector(self) -> Vector:
        unit_normal = Vector(np.array([1, 1, 1])).normalized
        return Vector(self.oct_stress_vector.vector.dot(unit_normal)*unit_normal)

    @property
    def oct_shear_stress_vector(self) -> Vector:
        return Vector(self.oct_stress_vector.vector - self.oct_normal_stress_vector.vector)

    @property
    def westergaard_coordinates(self):
        I1 = self.stress_tensor.first_invariant
        J2 = self.stress_deviator_tensor.second_invariant
        J3 = self.stress_deviator_tensor.third_invariant
        zeta = I1 / 3 ** 0.5
        rho = (2 * J2) ** 0.5
        theta = np.arccos(3 * 3 ** 0.5 / 2 * J3 / J2**(3/2))/3
        return np.array([zeta, rho, theta])
