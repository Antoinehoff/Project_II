from . import base
import numpy as np

class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):
        """
        Left and right down corner are fixed
        """
        return [(0, nely-1, 0), (0, nely-1, 1), (3*nelx/3, nely - 1, 0), (3*nelx/3, nely - 1, 1)]

    @staticmethod
    def set_forces(nelx, nely, f, params):
        """
        force on the uper side only
        """
        for x in range(nelx):
            for s in range(0,params.Nloads+1):
                f[2 * ((nely + 1) * x + s*nely/(params.Nloads+1)) + 1] = -1.0 /(params.Nloads-1)/nelx

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        """
        for x in range(nelx):
            for s in range(1,params.Nloads+1):
                lb[(nely) * x + s*nely/(params.Nloads+1)] = ub[(nely) * x + 1*nely/4]
                """
        for x in range(nelx):
            for s in range(0,params.Nloads+1):
                for y in range(s*nely/(params.Nloads+1),s*nely/(params.Nloads+1)+nely*1/40):
                    lb[(nely) * x + y] = ub[(nely) * x + y]

        return None
