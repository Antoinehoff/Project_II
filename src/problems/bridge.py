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
            f[2 * ((nely + 1) * x + 0) + 1] = -1.0 / nelx
        print("AMP = " + str(-1.0/nelx))

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        """
        uper side is set to be passive
        """
        for x in range(nelx):
            lb[0 + x * nely] = ub[0 + x * nely]
#        for y in range(nely/2,nely):
#            for x in range(nelx/5,4*nelx/5):
#                ub[y + x * nely] = lb[y + x * nely]
