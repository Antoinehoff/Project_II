from . import base


class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):
        return [(0,nely,0),(0,nely,1),(nelx,nely,0),(nelx,nely,1)]

    @staticmethod
    def set_forces(nelx, nely, f, params):
        for x in range(nelx):
			f[2 * ((nely + 1) * x + 0)+1] = -1.0/nelx

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        for x in range(nelx):
            lb[0 + x * nely] = ub[0 + x * nely]
#        for y in range(nely):
#            lb[y + 0 * nely] = ub[y + 0 * nely]
#            lb[y + (nelx-1) * nely] = ub[y + (nelx-1) * nely]
