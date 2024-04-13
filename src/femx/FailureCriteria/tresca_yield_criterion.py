import numpy as np

from FailureCriteria.failure_criterion import FailureCriterion


class TrescaYieldCriterion(FailureCriterion):

    def test_yield(self, point):
        zeta, rho, theta = FailureCriterion.westergaard_coordinates(point)
        so = self.material.yield_strength
        return 2**0.5 * rho * np.sin(theta + np.pi/3) - so
