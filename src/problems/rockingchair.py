from . import base
import numpy as np
import sys

check_config = False

class BoundaryConditions(base.BaseProblem):

    @staticmethod
    def get_fixed_nodes(nelx, nely, params):

        tab  = []

        #Contact point is fixed
#        tab.append((nelx/2.0,nely/2.0,0))
#        tab.append((nelx/2.0,nely/2.0,1))
        #Plateau is fixed
#        y    = nely/2
#        for x in range(nelx):
#            tab.append((x, y, 0))
#            tab.append((x, y, 1))

        #Arc in contact is fixed
        """
        peripheral_elements = onCircle_elements(nelx,nely)
        Nsector=params.Nsector
        for P in peripheral_elements:
            x = P[0]
            y = P[1]
            phi = np.arctan2(y-nely/2,x-nelx/2)
            if  phi < (1.0/2.0 + 1.0/(2*Nsector))*np.pi and phi > (1.0/2.0 - 1.0/(2*Nsector))*np.pi :
               tab.append((x, y, 0))
               tab.append((x, y, 1))
        """
        #The border circle is fixed
        peripheral_elements = onCircle_elements(nelx,nely)
        Nsector=params.Nsector
        for P in peripheral_elements:
            x = P[0]
            y = P[1]
            phi = np.arctan2(y-nely/2,x-nelx/2)
            if  phi < np.pi and phi > 0 :
#            or  phi < -3*np.pi/4 and phi > -np.pi :
               tab.append((x, y, 0))
               tab.append((x, y, 1))
        #Small circle inside only
#        for x in range(nelx):
#            for y in range(nely):
#                if abs(((x-(nelx-1)/2.0)**2+(y-(nely-1)/2.0)**2)**.5 - (nelx-1)/25.0) < .5 :
#                    tab.append((x, y, 0))
#                    tab.append((x, y, 1))
        return tab

    @staticmethod
    def set_forces(nelx, nely, f, params):
        """
        Force on all the peripheral elements

        peripheral_elements = onCircle_elements(nelx,nely)
        print(params.problemType)
        if params.problemType == "ProblemType.AppearanceWithMaxComplianceAndSymmetry" \
        or params.problemType == "ProblemType.ComplianceWithSymmetry":
            Nsector = params.Nsector
            ampx = 1.0/len(peripheral_elements)*Nsector
            ampy = 1.0/len(peripheral_elements)*Nsector + 1.0/len(peripheral_elements)
            #print("AMP = " +str(amp))
            for P in peripheral_elements:
                x = P[0]
                y = P[1]
                R = radial_direction(P, nelx, nely)
                phi = np.arctan2(y-nely/2,x-nelx/2)
                if  phi < (1.0/2.0 + 1.0/(2*Nsector))*np.pi and phi > (1.0/2.0 - 1.0/(2*Nsector))*np.pi :
                    f[2 * ((nely+1) * x + y) + 0, 0] += (ampx) * R[0]
                    f[2 * ((nely+1) * x + y) + 1, 0] += (ampy) * R[1]
                    if check_config:
                        print(str(P[0])+','+str(P[1])+','+str(f[2 * ((nely+1) * x + y) + 0, 0])+\
                      ','+str(f[2 * ((nely+1) * x + y) + 1, 0]))

        else :
            ampx = 1.0/len(peripheral_elements)
            ampy = 2.0/len(peripheral_elements)
            for P in peripheral_elements:
                x = P[0]
                y = P[1]
                R = radial_direction(P, nelx, nely)
                phi = np.arctan2(y-nely/2,x-nelx/2)
                if  phi < np.pi and phi > 0 :
                    f[2 * ((nely+1) * x + y) + 0, 0] += (ampx) * R[0]
                    f[2 * ((nely+1) * x + y) + 1, 0] += (ampy) * R[1]
                    if check_config:
                        print(str(P[0])+','+str(P[1])+','+str(f[2 * ((nely+1) * x + y) + 0, 0])+\
                      ','+str(f[2 * ((nely+1) * x + y) + 1, 0]))
        """
        #Force on the Plateau
        y = nely/2
        for x in range(nelx):
            f[2 * ((nely+1) * x + y) + 1, 0] += -1.0/nelx
            if check_config:
                print(str(x)+','+str(y)+','+str(f[2 * ((nely+1) * x + y) + 0, 0])+\
              ','+str(f[2 * ((nely+1) * x + y) + 1, 0]))

        #Force on the arms
#        y = 2*nely/5
#        for x in range(2*nelx/6,5*nelx/6):
#            f[2 * ((nely+1) * x + y) + 1, 0] += -0.5/nelx
#            if check_config:
#                print(str(x)+','+str(y)+','+str(f[2 * ((nely+1) * x + y) + 0, 0])+\
#              ','+str(f[2 * ((nely+1) * x + y) + 1, 0]))

        #Force on the back
        for y in range(nely/2):
            x = y
            if ((x + .5 -  nelx/2)**2 + (y + .5 - nely/2)**2)**.5 < nelx/2.0 + .5 :
                f[2 * ((nely+1) * x + y) + 0, 0] += 0.5/nelx
                f[2 * ((nely+1) * x + y) + 1, 0] += 0.5/nelx
                if check_config:
                    print(str(x)+','+str(y)+','+str(f[2 * ((nely+1) * x + y) + 0, 0])+\
                  ','+str(f[2 * ((nely+1) * x + y) + 1, 0]))
        if check_config : sys.exit()

    @staticmethod
    def set_passive_elements(nelx, nely, lb, ub, params):
        """
        Circle external elements are empty and border are full
        Center and ground points are set to be 1
        """
        #upper region
#        for x in range(nelx):
#            for y in range(nely/2):
#                ub[y + x * nely] = lb[y + x * nely]
        #Outter circle
        for x in range(nelx):
            for y in range(nely/2,nely):
                if ((x + .5 -  nelx/2)**2 + (y + .5 - nely/2)**2)**.5 > nelx/2.0 + .5 :
                    ub[y + x * nely] = lb[y + x * nely]
        #Plateau of the wheelchair
#        for x in range(nelx):
#            lb[nely/2 + x * nely] = ub[nely/2 + x * nely]

        peripheral_elements = onCircle_elements(nelx,nely)
        Nsector=params.Nsector
        for P in peripheral_elements:
            x = P[0]
            y = P[1]
            phi = np.arctan2(y-nely/2,x-nelx/2)
            if  phi < np.pi and phi > 0 \
            or  phi < -3*np.pi/4 and phi > -np.pi :
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
