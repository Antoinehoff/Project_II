from . import base


class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):
        return [(0,nely,0),(0,nely,1),(nelx,nely,0),(nelx,nely,1)\
               ,(nelx,0,0),(nelx,0,1),(0,0,0),(0,0,1)]

    @staticmethod
    def set_forces(nelx, nely, f, params):
		for x in range(nelx):
			f[2 * ((nely + 1) * x    + 0   ) + 1, 0] = -1./(nelx)
			f[2 * ((nely + 1) * x    + nely) + 1, 0] = 1./(nelx)
		for y in range(nely):
			f[2 * ((nely + 1) * 0    + y   ) + 0, 0] = 1./(nely)
			f[2 * ((nely + 1) * nelx + y   ) + 0, 0] = -1./(nely)

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        """
        for x in range(nelx):
            lb[0 + x * nely] = ub[0 + x * nely]
            lb[nely + x * nely - 1] = ub[nely + x * nely - 1]
        for y in range(nely):
            lb[y + 0 * nely] = ub[y + 0 * nely]
            lb[y + (nelx-1) * nely] = ub[y + (nelx-1) * nely]
        """
