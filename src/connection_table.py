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
import sys

# Third party libs
import numpy as np
import scipy.misc
import scipy

def index_to_position(index, nelx, nely):
	"""
	Convert the index of a element to the centroid of the element
	"""
	pos = np.array([(index / (nely))+.5, int(index%(nely))+.5])
	return pos

def position_to_index(pos, nelx, nely):
	"""
	Convert a position vector to the index of the element containing it
	"""
	index = int(pos[0])*nely + int(pos[1])
	return index

def add_symmetry_planes(a_array, c_array, a, c, empty = 0):
	"""
	Use to add a new symmetry plane to the problem
	"""
	if len(a_array)==0 or len(c_array)==0 :
		a_array = [np.array(a/np.sqrt(np.dot(a,a)))]
		c_array = [np.array(c)]
	else :
		a = a/np.sqrt(np.dot(a,a))
		a_array = np.append(a_array,[a],0)
		c_array = np.append(c_array,[c],0)
	return a_array, c_array

def get_symmetry_image(a_array, c_array, nelx, nely):
	"""
	Create an image with line showing the symmetry planes
	"""
	image = np.ones(nelx*nely)
	Xmin = index_to_position(0, nelx, nely)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	rmax = np.sqrt(np.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
	for i in range(nelx*nely):
		X_i = index_to_position(i, nelx, nely)
		if(in_domain(X_i, a_array, c_array)==0):
			image[i] = .5;
		if(in_domain(X_i, a_array, c_array)==2):
			image[i]=0;
	return (image.reshape(nely,nelx))

def in_domain(X, a_array, c_array):
	"""
	Check is a given point is inside or outside the design domain
	"""
	flag = 1
	for n in range(np.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		dist = np.dot(X-c,a)
		if(dist > 0.7):
			flag = 0
		if(abs(dist)<0.7):
			flag = 2
	return flag

def get_symmetric_element(index, a, c, nelx, nely):
	"""
	Return the index of the symmetric element w.r.t. a plane (a,c)
	"""
	x_i = index_to_position(index,nelx,nely)
	dist = np.dot(x_i-c,a)
	x_proj = x_i - dist * a
	x_sym = x_proj - dist * a
	index_sym = position_to_index(x_sym,nelx,nely)
	return index_sym

def construct_connection_table(a_array, c_array, nelx, nely):
	"""
	Simple algorithm O(nelx*nely) to construct the table containing
	for an element i, outside the design domain, its symmetric
	element j inside the design domain. Symmetry planes variant.
	"""
	connection_table = np.array(range(nelx*nely))
	for n in range(np.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		for i in range(nelx*nely):
			X_i = index_to_position(i,nelx,nely)
			index_sym = connection_table[i]
			if in_domain(X_i, [a], [c])!=1:
				index_sym = connection_table[get_symmetric_element(i,a,c,nelx,nely)]
			connection_table[i]=index_sym
#	print(connection_table.reshape(nelx,nely))
#	sys.exit('Stop')
	return connection_table


def construct_mapping_vector_wheel(Nsector, nelx, nely):
	"""
	Simple algorithm O(nelx*nely) to construct the table containing
	for an element i, outside the design domain, its symmetric
	element j inside the design domain. Symmetry sectors variant.
	"""
	mapping_vector = [[i] for i in range(nelx*nely)]
	phi = 2*np.pi/Nsector
	c = np.array([nelx/2,nely/2]); #Center of the circle
	for i in range(nelx*nely):
		X_i = index_to_position(i,nelx,nely) - c #All positions are centered to the circle
		r_i = np.sqrt(np.dot(X_i,X_i))
		if r_i <= nelx/2 :
			phi_i = np.arctan2(X_i[1],X_i[0])
			if phi_i < phi/2 :
				for j in range(Nsector):
					theta = phi_i + phi*j
					#Symmetry with others sectors
					if j > 0 :
						x_sym = r_i*np.cos(theta)
						y_sym = r_i*np.sin(theta)
						sym_pos = np.array([x_sym,y_sym]) + c
						sym_idx = position_to_index(sym_pos,nelx,nely)
						mapping_vector[i].append(sym_idx)
					#Middle sector symmetry
					x_sym = r_i*np.cos(theta + phi - 2 * phi_i)
					y_sym = r_i*np.sin(theta + phi - 2 * phi_i)
					sym_pos = np.array([x_sym,y_sym]) + c
					sym_idx = position_to_index(sym_pos,nelx,nely)
					mapping_vector[i].append(sym_idx)


	#print(np.array(mapping_vector).reshape(nelx,nely))
	#sys.exit('Stop')
	return mapping_vector

def construct_mapping_vector_rockingchair(Nsector, nelx, nely):
	"""
	Simple algorithm O(nelx*nely) to construct the table containing
	for an element i, outside the design domain, its symmetric
	element j inside the design domain. Symmetry sectors variant.
	"""
	if Nsector%2 == 0:
		sys.exit('ABORT : Please put an odd number of sector in json file')

	mapping_vector = [[i] for i in range(nelx*nely)]
	phi = np.pi/Nsector
	c = np.array([nelx/2,nely/2]); #Center of the circle
	for i in range(nelx*nely):
		X_i = index_to_position(i,nelx,nely) - c #All positions are centered to the circle
		r_i = np.sqrt(np.dot(X_i,X_i))
		if r_i <= nelx/2 :
			phi_i = np.arctan2(X_i[1],X_i[0])
			if phi_i > 0 and phi_i < phi:
				for j in range(Nsector):
					theta = phi_i + phi*j
					#Symmetry with others sectors
					if j > 0 :
						x_sym = r_i*np.cos(theta)
						y_sym = r_i*np.sin(theta)
						sym_pos = np.array([x_sym,y_sym]) + c
						sym_idx = position_to_index(sym_pos,nelx,nely)
						mapping_vector[i].append(sym_idx)
					#middle sector symmetry
					x_sym = r_i*np.cos(theta + phi - 2 * phi_i)
					y_sym = r_i*np.sin(theta + phi - 2 * phi_i)
					sym_pos = np.array([x_sym,y_sym]) + c
					sym_idx = position_to_index(sym_pos,nelx,nely)
					#mapping_vector[i].append(sym_idx)


	#print(np.array(mapping_vector).reshape(nelx,nely))
	#sys.exit('Stop')
	return mapping_vector

def construct_mapping_vector_shelf(params,nelx,nely):
	"""
	simple method to impose sector copy between the rows of
	the shelf vertically
	"""
	N=params.Nsector
	mapping_vector = [[i] for i in range(nelx*nely)]
	for i in range(nelx*nely):
		X_i = index_to_position(i,nelx,nely)
		if X_i[1]< nely*1.0/N :
			for j in range(N):
				#Symmetry with others sectors
				if j > 0 :
					x_sym = X_i[0]
					y_sym = X_i[1] + j * nely*1.0/N
					sym_pos = np.array([x_sym,y_sym])
					sym_idx = position_to_index(sym_pos,nelx,nely)
					if sym_idx not in mapping_vector[i] : mapping_vector[i].append(sym_idx)
			for index in mapping_vector[i] : mapping_vector[index] = mapping_vector[i]
		#Vertical symmetry
		if params.vert_sym :
			if X_i[0] < nelx/2.0 :
				x_sym = nelx - X_i[0]
				y_sym = X_i[1]
				sym_pos = np.array([x_sym,y_sym])
				sym_idx = position_to_index(sym_pos,nelx,nely)
				mapping_vector[i].append(sym_idx)
				mapping_vector[sym_idx].append(i)
	return mapping_vector

def construct_mapping_vector(connection_table):
	mapping_vector = [[] for i in range(len(connection_table))]
	for i in range(len(connection_table)):
		mapping_vector[connection_table[i]].append(i)
	mapping_vector = [mapping_vector[i] for i in range(len(connection_table)) if len(mapping_vector[i])>0]
	return mapping_vector

def get_sym_indep_indices(connection_table):
	indices=[]
	for i in range(len(connection_table)):
		if connection_table[i]==i:
			indices.append(i)
	return indices

def apply_symm_to_grad(grad,mapping_vector):
	grad_mean = grad
	for i in range(len(mapping_vector)):
		sublist=mapping_vector[i]
		S = len(sublist)
		mean = sum(grad[j] for j in sublist)/S
		for index in sublist:
			grad_mean[index]=mean
	return grad_mean

""" OLD VERSION (KOSAKA AND SWANN
	Implementation of the algorithm from Kosaka and Swan 1999 to construct the
	table containing for an element i, outside the design domain, its symmetric
	element j inside the design domain. connection_table[i]=j if i is outside the
	domain, =-1 if it is inside the domain or have no symmetric element.

def in_domain(X, a_array, c_array, epsilon, rmax):
	flag = True
	for n in range(np.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		if(np.dot(X-c,a) > epsilon*rmax):
			flag = False
	return flag

def in_domain(X, a_array, c_array):
	flag = 1
	for n in range(np.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		dist = np.dot(X-c,a)
		if(dist > 0.7):
			flag = 0
		if(abs(dist)<0.7):
			flag = 2
	return flag

def construct_connection_table(a_array, c_array, nelx, nely):
	Xmin = index_to_position(0, nelx, nely)
	print 'Xmin = {0}'.format(Xmin)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	print 'Xmax = {0}'.format(Xmax)
	rmax = np.sqrt(np.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
	print 'rmax = {0}'.format(rmax)
	epsilon_1 = 0.000000
	epsilon_2 = 0.0000001
	nmast = 0
	k=0 #Counter on elements in the interior of the plane
	m=0
	alpha_ij = 0
	connection_table = np.ones(nelx*nely)*-1 #initilaized at -1 everywhere
	for n in range(np.shape(a_array)[0]): #For each symmetry plane...
		a = a_array[n]
		c = c_array[n]
		for i in range(nelx*nely): #1st Loop over all elements
			X_i = index_to_position(i, nelx, nely) #Centroid of element
			if(not in_domain(X_i, a_array, c_array, epsilon_1, rmax)):
					k=k+1
					X_i_proj = X_i - (np.dot(X_i-c,a))*a
					dmin = rmax
					for j in range(nelx*nely):
						X_j = index_to_position(j, nelx, nely)
						if(np.dot(X_j-c,a) < rmax * epsilon_1):
							temp1 = np.dot(X_i_proj - X_i, X_i_proj - X_j)
							temp2 = np.linalg.norm(X_i_proj - X_i)*np.linalg.norm(X_i_proj - X_j) + epsilon_2
							alpha_ij = temp1 / temp2
							if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper
								if(abs(np.linalg.norm(X_i_proj - X_i)-np.linalg.norm(X_i_proj - X_j))<dmin):
									nmast=j
									dmin=abs(np.linalg.norm(X_i_proj - X_i)-np.linalg.norm(X_i_proj - X_j))
									m=m+1
									connection_table[i]=nmast
	print '#elements : {0} '.format(nelx*nely)
	print '#inside elements : {0}({1}%)'.format(k,int(100.0*k/nelx/nely))
	print '#dependant elements : {0}({1}%)'.format(m,int(100.0*m/nelx/nely))
	return connection_table.astype(int)

def construct_connection_table_parallel(a_array, c_array, nelx, nely):
	#Initialization phase
	start = time.time()
	Xmin = index_to_position(0, nelx, nely)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	rmax = np.sqrt(np.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
	epsilon_1 = 0.0000001
	epsilon_2 = 0.0000001
	nmast = 0
	alpha_ij = 0
	connection_table = np.ones(nelx*nely)*-1
	num_cores = multiprocessing.cpu_count()
	temp=[]
	for n in range(np.shape(a_array)[0]): #For each symmetry plane...
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
	a = a_array[n]
	c = c_array[n]
	nmast=-1
	X_i = index_to_position(i, nelx, nely) #Centroid of element
	if( not in_domain(X_i, a_array, c_array, epsilon_1, rmax)):
		X_i_proj = X_i - np.dot(X_i-c,a)*a
		dmin = rmax
		for j in range(nelx*nely):
			X_j = index_to_position(j, nelx, nely)
			if(i!=j and np.dot(X_j-c,a) < rmax * epsilon_1):
				temp1 = np.dot(X_i_proj - X_i, X_i_proj - X_j)
				temp2 = np.linalg.norm(X_i_proj - X_i)*np.linalg.norm(X_i_proj - X_j) + epsilon_2
				alpha_ij = temp1 / temp2
				if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper	if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper
					if(abs(np.linalg.norm(X_i_proj - X_i)-np.linalg.norm(X_i_proj - X_j))<dmin):
						nmast=j
						dmin=abs(np.linalg.norm(X_i_proj - X_i)-np.linalg.norm(X_i_proj - X_j))
	return nmast

"""
