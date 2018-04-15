# -*- coding: utf-8 -*-
# Third party libs
import numpy
import nlopt

# Local libs
import images
from common import ProblemType
from filter import Filter
from topopt import TopoptProblem
from appearance import AppearanceCL

### Added by Antoine Hoffmann EPFL 2018
# For parallel computing and perf. measurement
from joblib import Parallel, delayed
import multiprocessing
import time
###

def unwrap_self_parallelProcess(tab):
    return Solver.parallelProcess(tab)

class Solver(object):

    def __init__(self, nelx, nely, params, problem_type, bc, gui=None, sym_bool=False, epsilon_1=0.0001, epsilon_2=0.0001):
        """
        Allocate and initialize internal data structures.
        """
        n = nelx * nely
        self.nelx = nelx
        self.nely = nely
        self.opt = nlopt.opt(nlopt.LD_MMA, n)
        self.problem_type = problem_type

        # Alloc arrays
        self.x_phys = numpy.ones(n)
        self.x_rgb = numpy.ones((nelx, nely, 3))

        # Set bounds
        lb = numpy.array(params.densityMin * numpy.ones(n, dtype=float))
        ub = numpy.array(1.0 * numpy.ones(n, dtype=float))
        bc.set_passive_elements(nelx, nely, lb, ub, params.problemOptions)
        self.opt.set_upper_bounds(ub)
        self.opt.set_lower_bounds(lb)

        # Set stopping criteria
        self.opt.set_maxeval(params.maxSolverStep)
        self.opt.set_ftol_rel(0.0)

        # Setup topopt problem
        if problem_type.involves_compliance():
            self.problem = TopoptProblem(nelx, nely, params, bc)

        # Setup filter
        self.filtering = Filter(nelx, nely, params, problem_type)

        # Setup appearance matcher
        if problem_type.involves_appearance():
            self.exemplar_rgb = images.load_exemplar(params)
            self.appearance_cl = AppearanceCL(params.lambdaOccurrenceMap,
                                              params.exponentSimilarityMetric,
                                              params.appearanceNormWeight)

        # Setup user parameters
        self.params = params

        # Setup functional right-hand sides
        self.volume_max = self.params.volumeFracMax * n
        self.volume_min = self.params.volumeFracMin * n
        self.appearance_max = 0
        self.compliance_max = 0
        if problem_type.involves_compliance():
            if 'complianceMax' in params:
                self.compliance_max = params.complianceMax

        # Define optimization problem
        self.init_problem()

        # Set GUI callback
        self.gui = gui

        #### Added by Antoine Hoffmann EPFL 2018
        self.sym_bool = sym_bool #boolean to know if the problem has symmetry
        self.connection_table = numpy.ones(self.nelx*self.nely)*-1 #map symmetry
        self.a_array = self.params.a #sym. planes normal vectors
        self.c_array = self.params.c #sym. planes origin shifts
        self.epsilon_1 = epsilon_1
        self.epsilon_2 = epsilon_2
    	self.Xmin = numpy.array([.5, .5])
        self.Xmax = numpy.array([self.nelx-.5, self.nely-.5])
        self.rmax = numpy.sqrt(numpy.dot(self.Xmax-self.Xmin,self.Xmax-self.Xmin)) #Length scale of mesh
        ####

    def init_problem(self):
        """
        Define objective function, constraint functions, and constraints rhs.
        """
        # Set objective and constraint functions
        self.opt.remove_inequality_constraints()
        if self.problem_type == ProblemType.Appearance:
            self.opt.set_min_objective(self.appearance_function)
        elif self.problem_type == ProblemType.Compliance:
            self.opt.set_min_objective(self.compliance_function)
            self.opt.add_inequality_constraint(self.volume_max_function)
            self.opt.add_inequality_constraint(self.volume_min_function)
        elif self.problem_type == ProblemType.AppearanceWithMaxCompliance:
            self.opt.set_min_objective(self.appearance_function)
            self.opt.add_inequality_constraint(self.compliance_function)
            self.opt.add_inequality_constraint(self.volume_max_function)
            self.opt.add_inequality_constraint(self.volume_min_function)

            ####Added by Antoine Hoffmann EPFL 2018
        elif self.problem_type == ProblemType.ComplianceWithSymmetry:
            self.sym_bool = True
            self.opt.set_min_objective(self.compliance_function)
            self.opt.add_inequality_constraint(self.volume_max_function)
            self.opt.add_inequality_constraint(self.volume_min_function)
        elif self.problem_type == ProblemType.AppearanceWithMaxComplianceAndSymmetry:
            self.sym_bool = True
            self.opt.set_min_objective(self.appearance_function)
            self.opt.add_inequality_constraint(self.compliance_function)
            self.opt.add_inequality_constraint(self.volume_max_function)
            self.opt.add_inequality_constraint(self.volume_min_function)
            ##############

    def guess_initial_state(self):
        # Lower and upper bounds
        lb = self.opt.get_lower_bounds()
        ub = self.opt.get_upper_bounds()
        # Material budget
        passive = (lb == ub)
        num_passive = numpy.sum(passive)
        num_active = self.nelx * self.nely - num_passive
        budget = self.volume_max - numpy.sum(lb)
        active_frac = budget / num_active
        if active_frac < self.params.densityMin or active_frac > 1:
            assert False, "Unsatisfiable volume constraint"
        # Initial guess for non passive elements
        x = numpy.ones(self.nely * self.nelx, dtype=float) * active_frac
        return numpy.clip(x, lb, ub)

    def enforce_volume_constraint(self, x, volume_function):
        """
        Given a problem instance with a violated volume constraint g_v, an objective f,
        and possibly other constraints g_i, tries to minimize g_v under the constraint that
        the other functionals (excluding appearance) should not increase.
        """
        # If there is no volume constraint, or if it is satisfied already, do nothing
        if not self.problem_type.has_volume_constraint() or volume_function(x) <= 0:
            return x

        print("Enforce volume constraint")

        # Otherwise, define a new subproblem to minimize
        self.opt.remove_inequality_constraints()
        old_stopval = self.opt.get_stopval()
        self.opt.set_stopval(0)

        # Objective: volume
        self.opt.set_min_objective(volume_function)

        # Constraint: compliance
        old_compliance_max = self.compliance_max
        if self.problem_type.involves_compliance():
            self.compliance_max = 0
            self.compliance_max = self.compliance_function(x)
            self.opt.add_inequality_constraint(self.compliance_function)

        # Solve subproblem
        x = self.opt.optimize(x)

        # Restore previous state
        self.compliance_max = old_compliance_max
        self.opt.set_stopval(old_stopval)
        self.init_problem()
        return x

    def enforce_compliance_constraint(self, x):
        """
        Given a problem instance with a violated compliance constraint g_c, an objective f,
        and possibly other constraints g_i, tries to minimize g_c under the constraint that
        the other functionals (excluding appearance) should not increase.
        """

        # If there is no volume constraint, or if it is satisfied already, do nothing
        if not self.problem_type.has_compliance_constraint() or self.compliance_function(x) <= 0:
            return x

        print("Enforce compliance constraint")

        # Otherwise, define a new subproblem to minimize
        self.opt.remove_inequality_constraints()
        old_stopval = self.opt.get_stopval()
        self.opt.set_stopval(0)

        # Objective: compliance
        self.opt.set_min_objective(self.compliance_function)

        # Constraint: volume
        old_volume_max = self.volume_max
        old_volume_min = self.volume_min
        if self.problem_type.involves_volume():
            self.volume_max = 0
            self.volume_max = self.volume_max_function(x)
            self.opt.add_inequality_constraint(self.volume_max_function)
            self.volume_min = 0
            self.volume_min = self.volume_min_function(x)
            self.opt.add_inequality_constraint(self.volume_min_function)

        # Solve subproblem
        x = self.opt.optimize(x)

        # Restore previous state
        self.volume_max = old_volume_max
        self.volume_min = old_volume_min
        self.opt.set_stopval(old_stopval)
        self.init_problem()
        return x

    def optimize(self, x, enforce_constraints=True):
        print("* " + str(self.problem_type))
        lb = self.opt.get_lower_bounds()
        ub = self.opt.get_upper_bounds()

        # If first attempt, guess initial state
        if x is None:
            x = self.guess_initial_state()

        # Make sure bounds are strictly satisfied (avoid Nlopt crash)
        x = numpy.clip(x, lb, ub)
        print("* Initial volume = " + str(sum(x) / float(len(x))))

        # Enforce violated constraint by solving alternative subproblems
        if enforce_constraints:
            self.enforce_volume_constraint(x, self.volume_max_function)
            self.enforce_volume_constraint(x, self.volume_min_function)
            x = self.enforce_compliance_constraint(x)
            print("Enforcing constraints: done")

        ###Added by Antoine Hoffmann EPFL 2018
        if self.sym_bool:
            self.enforce_symmetry_constraint(x,self)
        ###

        # Launch optimization
        x = self.opt.optimize(x)
        print("* Last optimum value = " + str(self.opt.last_optimum_value()))
        return x

    def last_optimum_value(self):
        return self.opt.last_optimum_value()

    def compliance_function(self, x, dc=None):
        """
        Compliance function. When used as a constraint, corresponds to the inequality:
        >>> u.f - c_max ≤ 0
        """
        # Filter design variables
        self.filtering.filter_variables(x, self.x_phys)

        # Display physical variables
        if self.gui:
            self.gui.update(self.x_phys)

        # Setup and solve FE problem
        self.problem.compute_displacements(self.x_phys)

        # Compliance and sensitivity
        obj = self.problem.compute_compliance(self.x_phys, dc)

        # Sensitivity filtering
        if dc is not None:
            self.filtering.filter_compliance_sensitivities(self.x_phys, dc)

        print("- Compliance = %.3f" % (obj))

        return obj - self.compliance_max

    def volume_max_function(self, x, dv=None):
        """
        Maximum volume function. When used as a constraint, corresponds to the inequality:
        >>> volume(x_phys) - v_max ≤ 0
        """
        # Filter design variables
        self.filtering.filter_variables(x, self.x_phys)

        if dv is not None:
            # Volume sensitivities
            dv[:] = 1.0

            # Sensitivity filtering
            self.filtering.filter_volume_sensitivities(self.x_phys, dv)

        vol = sum(self.x_phys)

        print("- Volume = %.3f %%" % (vol * 100.0 / float(len(x))))

        return vol - self.volume_max

    def volume_min_function(self, x, dv=None):
        """
        Minimum volume function. When used as a constraint, corresponds to the inequality:
        >>> -volume(x_phys) + v_min ≤ 0
        """
        # Filter design variables
        self.filtering.filter_variables(x, self.x_phys)

        if dv is not None:
            # Volume sensitivities
            dv[:] = -1.0

            # Sensitivity filtering
            self.filtering.filter_volume_sensitivities(self.x_phys, dv)

        vol = sum(self.x_phys)

        print("- Volume = %.3f %%" % (vol * 100.0 / float(len(x))))

        return -vol + self.volume_min

    def appearance_function(self, x, da=None):
        """
        Appearance function. When used as a constraint, corresponds to the inequality:
        >>> app(x_phys) - a_max ≤ 0
        """
        # Filter design variables
        self.filtering.filter_variables(x, self.x_phys)

        # Display physical variables
        if self.gui:
            self.gui.update(self.x_phys)

        # Appearance and its gradient
        self.x_rgb[:] = numpy.reshape(
            numpy.repeat(self.x_phys.reshape((self.nelx, self.nely)), 3, axis=1), (self.nelx, self.nely, 3))
        patch_size = (self.params.neighborhoodSize, self.params.neighborhoodSize)
        pm_iter = self.params.patchMatchIter
        if da is not None:
            grad_reshape = da.reshape((self.nelx, self.nely))
        else:
            grad_reshape = None
        sim = self.appearance_cl.compute(self.x_rgb, self.exemplar_rgb, grad_reshape, patch_size, pm_iter)

        if da is not None:
            # Sensitivity filtering
            self.filtering.filter_appearance_sensitivities(self.x_phys, da)

        print("- Appearance = %.3f" % (sim))
        return sim - self.appearance_max

    ### Symmetry functions added by Antoine Hoffmann EPFL 2018

    def index_to_position(self, index):
    	"""
    	Convert the index of a element to the centroid of the element
    	"""
    	return numpy.array([(index % self.nelx)+.5, int(index/self.nelx)+.5])

    def add_symmetry_planes(self, a, c):
    	"""
    	Use to add a new symmetry plane to the problem
    	"""
        a = a/numpy.sqrt(numpy.dot(a,a))#normalize the vector
    	if len(self.a_array)==0 or len(self.c_array)==0 :
    		self.a_array = [numpy.array(a)]
    		self.c_array = [numpy.array(c)]
    	else :
    		self.a_array = numpy.append(self.a_array,[a],0)
    		self.c_array = numpy.append(self.c_array,[c],0)

    def in_domain(self, X):
    	"""
    	Check is a given point is inside or outside the design domain
    	"""
    	flag = True
    	for n in range(numpy.shape(self.a_array)[0]):
    		a = self.a_array[n]
    		c = self.c_array[n]
    		if(numpy.dot(X-c,a) > self.epsilon_1*self.rmax):
    			flag = False
    	return flag

    def parallelProcess(self, tab):
    	"""
    	Task to be done in parallel for construct_connection_table_parallel
    	"""
        n = tab[0]
        i = tab[1]
        a = self.a_array[n]
        c = self.c_array[n]
    	nmast=-1
    	X_i = self.index_to_position(i) #Centroid of element
    	if( not self.in_domain(X_i)):
    			X_i_proj = X_i - (numpy.dot(X_i-c,a))*a
    			dmin = self.rmax
    			for j in range(self.nelx*self.nely):
    				X_j = self.index_to_position(j)
    				if(i!=j and numpy.dot(X_j-c,a) < self.rmax * self.epsilon_1):
    					temp1 = numpy.dot(X_i_proj - X_i, X_i_proj - X_j)
    					temp2 = numpy.linalg.norm(X_i_proj - X_i)*numpy.linalg.norm(X_i_proj - X_j) + self.epsilon_2
    					alpha_ij = temp1 / temp2
    					if(alpha_ij < self.epsilon_2-1): #Not the same as in Kosaka and Swan paper	if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper
    						if(abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))<dmin):
    							nmast=j
    							dmin=abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))
    	return nmast

    def construct_connection_table_parallel(self):
    	"""
    	Adaptation of the function construct_connection_table but with a parallelization
    	of the main loop.
    	"""
    	nmast = 0
    	alpha_ij = 0
    	num_cores = multiprocessing.cpu_count()
    	temp=[]
        nel=self.nelx*self.nely
    	for n in range(numpy.shape(self.a_array)[0]): #For each symmetry plane...
    		a = self.a_array[n]
    		c = self.c_array[n]
#    		temp = Parallel(n_jobs=num_cores)\
#    		                  (delayed(unwrap_self_parallelProcess)  \
#    		                  (self,i)\
#                              for i in zip([n]*nel,range(nel)))
#    		                  for i in range(self.nelx*self.nely)\
#                             for i in zip([self]*nel,zip([n]*nel,range(nel))))

                for i in range(nel):
    			if (temp[i]>-1) : self.connection_table[i]=temp[i]

    def construct_connection_table(self):
        a_array = self.params.a
        c_array = self.params.c
        for i in range(len(self.params.a)):
            self.add_symmetry_planes(a_array[i],c_array[i])
    	alpha_ij = 0
    	for n in range(numpy.shape(self.a_array)[0]): #For each symmetry plane...
    		a = self.a_array[n]
    		c = self.c_array[n]
    		for i in range(self.nelx*self.nely): #1st Loop over all elements
    			X_i = self.index_to_position(i) #Centroid of element
    			if(not self.in_domain(X_i)):
    					X_i_proj = X_i - (numpy.dot(X_i-c,a))*a
    					dmin = self.rmax
    					for j in range(self.nelx*self.nely):
    						X_j = self.index_to_position(j)
    						if(numpy.dot(X_j-c,a) < self.rmax * self.epsilon_1):
    							temp1 = numpy.dot(X_i_proj - X_i, X_i_proj - X_j)
    							temp2 = numpy.linalg.norm(X_i_proj - X_i)*numpy.linalg.norm(X_i_proj - X_j) + self.epsilon_2
    							alpha_ij = temp1 / temp2
    							if(alpha_ij < self.epsilon_2-1): #Not the same as in Kosaka and Swan paper
    								if(abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))<dmin):
    									nmast=j
    									dmin=abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))
    									self.connection_table[i]=nmast

    def enforce_symmetry_constraint(x,self):
        for i in range(self.n):
            if self.connection_table[i]> -1:
                self.x[i]=self.x[self.connection_table[i]]


    def get_symmetry_image(self):
    	"""
    	Create an image with line showing the symmetry planes
    	"""
    	image = numpy.ones(self.nelx*self.nely)
    	for i in range(self.nelx*self.nely):
    		X_i = index_to_position(i)
    #		if(numpy.dot(X_i-c,a) > 0):
    		if(not in_domain(X_i, self.a_array, self.c_array, 0.000001, rmax)):
    			image[i] = .5;
    	return (image.reshape(self.nely,self.nelx))
    ###
