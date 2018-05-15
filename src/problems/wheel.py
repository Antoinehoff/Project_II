from . import base
import numpy as np
import sys

check_config = False

class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):
        """
        Center and ground only are fixed
        """
        return [(nelx/2,nely/2,0),(nelx/2,nely/2,1)]

    @staticmethod
    def set_forces(nelx, nely, f, params):
        """
        Force on all the peripheral elements
        """
        peripheral_elements = onCircle_elements(nelx,nely)
        amp = 1.0/len(peripheral_elements)
        print("AMP = " +str(amp))
        for P in peripheral_elements:
            x = P[0]
            y = P[1]
            R = radial_direction(P, nelx, nely)
            f[2 * ((nely+1) * x + y) + 0, 0] = (amp) * R[0]
            f[2 * ((nely+1) * x + y) + 1, 0] = (amp) * R[1]
            if check_config:
                print(str(P[0])+','+str(P[1])+','+str(f[1 * ((nely+1) * x + y) + 1, 0])+\
              ','+str(f[2 * ((nely+1) * x + y) + 1, 0]))
#        print("||F||**2 = "+str(np.dot(np.array(f[:,0]),np.array(f[:,0]))))
        if check_config : sys.exit()

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        """
        Circle external elements are empty and border are full
        Center and ground points are set to be 1
        """
        lb[nely/2 + nelx/2 * nely] = ub[(nely/2+1) + (nelx/2+1) * nely]
#        lb[(nely-1) + (nelx/2+1) * nely] = ub[(nely-1) + (nelx/2+1) * nely]
        for x in range(nelx):
            for y in range(nely):
                if ((x-  nelx/2)**2 + (y - nely/2)**2)**.5 > nelx/2.0 + .5 :
                    ub[y + x * nely] = lb[y + x * nely]
                if abs(((x-(nelx-1)/2.0)**2+(y-(nely-1)/2.0)**2)**.5 - (nelx-1)/2.0) < .5 :
                    lb[y + x * nely] = ub[y + x * nely]


def onCircle_elements(nelx,nely):
    """
    return a list of the element centroid located on the border
    of a circle of radius nelx
    """
    list = []
    for x in range(nelx):
        for y in range(nely):
            if abs(((x-(nelx-1)/2.0)**2+(y-(nely-1)/2.0)**2)**.5 - (nelx-1)/2.0) < .5 :
                list.append((x,y))
    return list

def radial_direction(point, nelx, nely):
    """
    compute the force normal to the circle
    """
    (x,y) = point
    f = np.array([(nelx-1)/2.0 - x, (nely-1)/2.0 - y])
    f = f/np.sqrt(np.dot(f,f))
    return f
