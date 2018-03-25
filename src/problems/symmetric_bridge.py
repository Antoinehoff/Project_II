from . import base


class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):
        return ([0, nely, 0], [0, nely, 1], [nelx, nely, 0],[nelx, nely, 1])
        
    @staticmethod
    def get_symmetry_plane(nelx,nely):
		a = [1, 0]
		c = int(nelx/2)
		return (numpy.array(a), c)

    @staticmethod
    def set_forces(nelx, nely, f, params):
		for x in range(nelx):
			f[2 * int((nely + 1) * x + nely *0) + 1, 0] = -2/nelx

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        pass
