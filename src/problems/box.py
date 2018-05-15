from . import base
import numpy as np

class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):
        """
        Fixation only on the down side
        """
        tab0 = []
        tab1 = []
        tab  = []
        tab0 = [(x, nely - 1, 0) for x in range(nelx)]
        tab1 = [(x, nely - 1, 1) for x in range(nelx)]
        for i in range(nelx):
            tab.append(tab0[i])
            tab.append(tab1[i])
        return tab

    @staticmethod
    def set_forces(nelx, nely, f, params):
        """
        force only on the upper side
        """
        for x in range(nelx):
            f[2 * ((nely + 1) * x + 0) + 1, 0] = -1. / (nelx)
        print("||F||**2 = "+str(np.dot(np.array(f[:,0]),np.array(f[:,0]))))
#			f[2 * ((nely + 1) * x    + nely) + 1, 0] = 1./(nelx)
#		for y in range(nely):
#			f[2 * ((nely + 1) * 0    + y   ) + 0, 0] = 1./(nely)
#			f[2 * ((nely + 1) * nelx + y   ) + 0, 0] = -1./(nely)

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        """
        the 4 sides are set to be passive
        """
        for x in range(nelx):
            lb[0 + x * nely] = ub[0 + x * nely]
            lb[nely + x * nely - 1] = ub[nely + x * nely - 1]
        for y in range(nely):
            lb[y + 0 * nely] = ub[y + 0 * nely]
            lb[y + (nelx - 1) * nely] = ub[y + (nelx - 1) * nely]
