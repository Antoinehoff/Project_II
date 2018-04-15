#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions to create the connection table in view to enforce symmetry in the
topological optimization
Antoine Hoffmann EPFL 2018
"""
# System libs
import os
import json
import argparse
import importlib

# Third party libs
import numpy
import scipy.misc
import scipy

# For parallel computing
from joblib import Parallel, delayed
import multiprocessing
import time


def index_to_position(index, nelx, nely):
	"""
	Convert the index of a element to the centroid of the element
	"""
	return numpy.array([(index % nelx)+.5, int(index/nelx)+.5])

def add_symmetry_planes(a_array, c_array, a, c, empty = 0):
	"""
	Use to add a new symmetry plane to the problem
	"""
	if len(a_array)==0 or len(c_array)==0 :
		a_array = [numpy.array(a/numpy.sqrt(numpy.dot(a,a)))]
		c_array = [numpy.array(c)]
	else :
		a = a/numpy.sqrt(numpy.dot(a,a))
		a_array = numpy.append(a_array,[a],0)
		c_array = numpy.append(c_array,[c],0)
	return a_array, c_array

def get_symmetry_image(a_array, c_array, nelx, nely):
	"""
	Create an image with line showing the symmetry planes
	"""
	image = numpy.ones(nelx*nely)
	Xmin = index_to_position(0, nelx, nely)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	rmax = numpy.sqrt(numpy.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
	for i in range(nelx*nely):
		X_i = index_to_position(i, nelx, nely)
		if(not in_domain(X_i, a_array, c_array, 0.000001, rmax)):
			image[i] = .5;
	return (image.reshape(nely,nelx))

def in_domain(X, a_array, c_array, epsilon, rmax):
	"""
	Check is a given point is inside or outside the design domain
	"""
	flag = True
	for n in range(numpy.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		if(numpy.dot(X-c,a) > epsilon*rmax):
			flag = False
	return flag


def construct_connection_table(a_array, c_array, nelx, nely):
	"""
	Implementation of the algorithm from Kosaka and Swan 1999 to construct the
	table containing for an element i, outside the design domain, its symmetric
	element j inside the design domain. connection_table[i]=j if i is outside the
	domain, =-1 if it is inside the domain or have no symmetric element.
	"""
	Xmin = index_to_position(0, nelx, nely)
	print 'Xmin = {0}'.format(Xmin)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	print 'Xmax = {0}'.format(Xmax)
	rmax = numpy.sqrt(numpy.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
	print 'rmax = {0}'.format(rmax)
	epsilon_1 = 0.000000
	epsilon_2 = 0.0000001
	nmast = 0
	k=0 #Counter on elements in the interior of the plane
	m=0
	alpha_ij = 0
	connection_table = numpy.ones(nelx*nely)*-1 #initilaized at -1 everywhere
	for n in range(numpy.shape(a_array)[0]): #For each symmetry plane...
		a = a_array[n]
		c = c_array[n]
		for i in range(nelx*nely): #1st Loop over all elements
			X_i = index_to_position(i, nelx, nely) #Centroid of element
			if(not in_domain(X_i, a_array, c_array, epsilon_1, rmax)):
					k=k+1
					X_i_proj = X_i - (numpy.dot(X_i-c,a))*a
					dmin = rmax
					for j in range(nelx*nely):
						X_j = index_to_position(j, nelx, nely)
						if(numpy.dot(X_j-c,a) < rmax * epsilon_1):
							temp1 = numpy.dot(X_i_proj - X_i, X_i_proj - X_j)
							temp2 = numpy.linalg.norm(X_i_proj - X_i)*numpy.linalg.norm(X_i_proj - X_j) + epsilon_2
							alpha_ij = temp1 / temp2
							if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper
								if(abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))<dmin):
									nmast=j
									dmin=abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))
									m=m+1
									connection_table[i]=nmast
	print '#elements : {0} '.format(nelx*nely)
	print '#inside elements : {0}({1}%)'.format(k,int(100.0*k/nelx/nely))
	print '#dependant elements : {0}({1}%)'.format(m,int(100.0*m/nelx/nely))
	return connection_table.astype(int)

def construct_connection_table_parallel(a_array, c_array, nelx, nely):
	"""
	Adaptation of the function construct_connection_table but with a parallelization
	of the main loop.
	"""
	#Initialization phase
	start = time.time()
	Xmin = index_to_position(0, nelx, nely)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	rmax = numpy.sqrt(numpy.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
	epsilon_1 = 0.0000001
	epsilon_2 = 0.0000001
	nmast = 0
	alpha_ij = 0
	connection_table = numpy.ones(nelx*nely)*-1
	num_cores = multiprocessing.cpu_count()
	temp=[]
	for n in range(numpy.shape(a_array)[0]): #For each symmetry plane...
		a = a_array[n]
		c = c_array[n]
		temp = Parallel(n_jobs=num_cores)\
		                  (delayed(parallelProcess)  \
		                  (a_array, c_array, n, nelx, nely, rmax, epsilon_1, epsilon_2, i)\
		                  for i in range(nelx*nely))
		for i in range(nelx*nely):
			if (temp[i]>-1) : connection_table[i]=temp[i]
	end = time.time()
	print('Ellapsed time : ' + str(end-start)+'[s]')
	return connection_table.astype(int)

def parallelProcess(a_array, c_array, n, nelx, nely, rmax, epsilon_1, epsilon_2, i):
	"""
	Task to be done in parallel for construct_connection_table_parallel
	"""
	a = a_array[n]
	c = c_array[n]
	nmast=-1
	X_i = index_to_position(i, nelx, nely) #Centroid of element
	if( not in_domain(X_i, a_array, c_array, epsilon_1, rmax)):
		X_i_proj = X_i - numpy.dot(X_i-c,a)*a
		dmin = rmax
		for j in range(nelx*nely):
			X_j = index_to_position(j, nelx, nely)
			if(i!=j and numpy.dot(X_j-c,a) < rmax * epsilon_1):
				temp1 = numpy.dot(X_i_proj - X_i, X_i_proj - X_j)
				temp2 = numpy.linalg.norm(X_i_proj - X_i)*numpy.linalg.norm(X_i_proj - X_j) + epsilon_2
				alpha_ij = temp1 / temp2
				if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper	if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper
					if(abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))<dmin):
						nmast=j
						dmin=abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))
	return nmast
