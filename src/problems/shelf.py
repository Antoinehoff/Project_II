from . import base
import numpy as np

class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):
        """
        Left and right down corner are fixed
        """
        return [(0, nely - 1, 0), (0, nely - 1, 1), (nelx - 1, nely - 1, 0), (nelx - 1, nely - 1, 1)]

    @staticmethod
    def set_forces(nelx, nely, f, params):
        """
        force on the uper side only
        """
        for x in range(nelx):
            f[2 * ((nely + 1) * x + 1*nely/4) + 1] = -1.0 / nelx
            f[2 * ((nely + 1) * x + 2*nely/4) + 1] = -1.0 / nelx
            f[2 * ((nely + 1) * x + 3*nely/4) + 1] = -1.0 / nelx

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        for x in range(nelx):
            lb[(nely) * x + 1*nely/4] = ub[(nely) * x + 1*nely/4]
            lb[(nely) * x + 2*nely/4] = ub[(nely) * x + 2*nely/4]
            lb[(nely) * x + 3*nely/4] = ub[(nely) * x + 3*nely/4]
        return None
