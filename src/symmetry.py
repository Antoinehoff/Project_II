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
import solver

class Symmetric_Solver(Solver): 
    
	def design_variable_merge(self, x, nelx, nely):
		"""
		Implementation of the algorithm from Kosaka and Swan 1999
		"""
		return x
	
	def index_to_position(index, nelx, nely):
		"""
		Convert the index of a element to a position on the grid
		"""
		
		return np.array([index % nely, int(index / nelx)])
