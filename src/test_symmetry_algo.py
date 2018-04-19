#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test script for additional functions
Antoine Hoffmann EPFL 2018
"""

# System libs
import os
import os.path
import io
import json
import pickle
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
	return numpy.array([(index % nelx)+.5, index/nelx+.5])

def position_to_index(v, nelx, nely):
	"""
	Convert a position vector to the index of the element containing it
	"""
	return int(v[0]) + int(v[1])*nelx


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
	for n in range(numpy.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		dist = numpy.dot(X-c,a)
		if(dist > 0.7):
			flag = 0
		if(abs(dist)<0.7):
			flag = 2
	return flag

def get_symmetric_element(index, a, c, nelx, nely):
	"""
	Return the index of the symmetric element w.r.t. a plane a,c
	"""
	x_i = index_to_position(index,nelx,nely)
	dist = numpy.dot(x_i-c,a)
	x_proj = x_i - dist * a
	x_sym = x_proj - dist * a
	index_sym = position_to_index(x_sym,nelx,nely)
	return index_sym

def construct_connection_table(a_array, c_array, nelx, nely):
	"""
	Implementation of the algorithm from Kosaka and Swan 1999 to construct the
	table containing for an element i, outside the design domain, its symmetric
	element j inside the design domain. connection_table[i]=j if i is outside the
	domain, =-1 if it is inside the domain or have no symmetric element.
	"""
	filename = 'connection_tables/connection_table'
	for i in range(a_array.shape[0]):
		filename += '_' + str(a_array[i])
	filename = '.txt'
	print(filename)

	if not os.path.exists(filename):
		Xmin = index_to_position(0, nelx, nely)
		print 'Xmin = {0}'.format(Xmin)
		Xmax = index_to_position(nelx*nely-1, nelx, nely)
		print 'Xmax = {0}'.format(Xmax)
		rmax = numpy.sqrt(numpy.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
		print 'rmax = {0}'.format(rmax)
		epsilon_1 = 0.001
		epsilon_2 = 0.001
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
				if(not in_domain(X_i, a_array, c_array)):
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
		# Write Pickle file
		with open(filename, 'w') as f:
			pickle.dump(connection_table, f)
	else :
		# Read Pickle file
		with open(filename) as f:
		    connection_table = pickle.load(f)
	return connection_table.astype(int)

def construct_connection_table_parallel(a_array, c_array, nelx, nely):
	"""
	Adaptation of the function construct_connection_table but with a parallelization
	of the main loop.
	"""
	filename = 'connection_tables/connection_table_'+str(nelx)+'x'+str(nely)
	for i in range(len(a_array)):
		filename += '_' + str(a_array[i])
	filename += '.json'
	print(filename)
	if not os.path.exists(filename):
		Xmin = index_to_position(0, nelx, nely)
		print 'Xmin = {0}'.format(Xmin)
		Xmax = index_to_position(nelx*nely-1, nelx, nely)
		print 'Xmax = {0}'.format(Xmax)
		rmax = numpy.sqrt(numpy.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
		print 'rmax = {0}'.format(rmax)
		epsilon_1 = 0.001
		epsilon_2 = 0.001
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
				"""
		# Write Pickle file
		with open(filename, 'w') as f:
			pickle.dump(connection_table, f)
	else :
		# Read Pickle file
		with open(filename) as f:
		    connection_table = pickle.load(f)
			"""
	return connection_table.astype(int)

def parallelProcess(a_array, c_array, n, nelx, nely, rmax, epsilon_1, epsilon_2, i):
	"""
	Task to be done in parallel for construct_connection_table_parallel
	"""
	a = a_array[n]
	c = c_array[n]
	nmast=-1
	X_i = index_to_position(i, nelx, nely) #Centroid of element
	X_i_proj = X_i - numpy.dot(X_i-c,a)*a # projection to the plane
	if i==0:
		print('Element '+ str(i) + ' position : '+ str(X_i))
		print('Element '+ str(i) + ' projection : '+ str(X_i_proj))
	if(in_domain(X_i, a_array, c_array, epsilon_1, rmax)==0):
			dmin = rmax
			for j in range(nelx*nely):
				X_j = index_to_position(j, nelx, nely)
				if((i!=j and numpy.dot(X_j-c,a) < rmax * epsilon_1)):
					temp1 = numpy.dot(X_i_proj - X_i, X_i_proj - X_j)
					temp2 = numpy.linalg.norm(X_i_proj - X_i)*numpy.linalg.norm(X_i_proj - X_j) + epsilon_2
					alpha_ij = temp1 / temp2
					if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper	if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper
						if(abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))<dmin):
							nmast=j
							dmin=abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))
	return nmast

def construct_connection_table_parallel_dummy(a_array, c_array, nelx, nely):
	"""
	Simple algorithm O(nelx*nely) to construct the table containing
	for an element i, outside the design domain, its symmetric
	element j inside the design domain. connection_table[i]=j if i is outside the
	domain, =-1 if it is inside the domain or have no symmetric element.
	"""
	num_cores = multiprocessing.cpu_count()
	connection_table = numpy.array(range(nelx*nely))
	with Parallel(n_jobs=1) as parallel:
		for n in range(numpy.shape(a_array)[0]):
			print(len(a_array))
			a = a_array[n]
			c = c_array[n]
			connection_table = parallel(delayed(parallelProcess_dummy)\
			               (a, c, nelx, nely, connection_table, i)\
			                for i in range(nelx*nely))
	return connection_table

def parallelProcess_dummy(a, c, nelx, nely, connection_table, i):
	X_i = index_to_position(i,nelx,nely)
	index_sym = connection_table[i]
	if in_domain(X_i, [a], [c])!=1:
			index_sym = connection_table[get_symmetric_element(i,a,c,nelx,nely)]
	return index_sym

def construct_connection_table_dummy(a_array, c_array, nelx, nely):
	"""
	Simple algorithm O(nelx*nely) to construct the table containing
	for an element i, outside the design domain, its symmetric
	element j inside the design domain. connection_table[i]=j if i is outside the
	domain, =-1 if it is inside the domain or have no symmetric element.
	"""
	num_cores = multiprocessing.cpu_count()
	connection_table = numpy.array(range(nelx*nely))
	for n in range(numpy.shape(a_array)[0]):
		print(len(a_array))
		a = a_array[n]
		c = c_array[n]
		for i in range(nelx*nely):
			X_i = index_to_position(i,nelx,nely)
			index_sym = connection_table[i]
			if in_domain(X_i, [a], [c])!=1:
					index_sym = connection_table[get_symmetric_element(i,a,c,nelx,nely)]
			connection_table[i]=index_sym
	return connection_table

def get_sym_indep_indices(connection_table):
	indices=[]
	for i in range(len(connection_table)):
		if connection_table[i]==i:
			indices.append(i)
	return indices

def construct_mapping_vector(connection_table):
	mapping_vector = [[] for i in range(len(connection_table))]
	for i in range(len(connection_table)):
		mapping_vector[connection_table[i]].append(i)

	mapping_vector = [mapping_vector[i] for i in range(len(connection_table)) if len(mapping_vector[i])>0]
	return mapping_vector

"""TESTS"""
def impose_symmetry_constraint(connection_table, a, c, nelx, nely):
	num_cores = multiprocessing.cpu_count()
	temp = Parallel(n_jobs=num_cores)\
		                  (delayed(parallelProcess_dummy)  \
		                  (a, c, nelx, nely, connection_table, i)\
		                  for i in range(nelx*nely))
	connection_table = temp
	return connection_table

def parallelProcess_test(a, c, nelx, nely, connection_table, i):
	X_i = index_to_position(i,nelx,nely)
	index_sym = connection_table[i]
	if in_domain(X_i, [a], [c])!=1:
			index_sym = connection_table[get_symmetric_element(i,a,c,nelx,nely)]
	return index_sym
""""""

def construct_test_image(nelx,nely):
		img_list = numpy.ones(nelx*nely)
		numpy.put(img_list,range(0*nelx/10 + nely/10 * nelx,10*nelx/10 + nely/10 * nelx),0)
		for i in range(nely):
			img_list[5*nelx/10 + i*nelx]=0
			img_list[5*nelx/10 + i*nelx-1]=0
			img_list[2*nelx/10 + i*nelx]=0
		img_list[5*nelx/10 + 4*nelx/10 * nelx]=0
		img_list[7*nelx/10 + nely/10 * nelx]=0
		img_list[6*nelx/10 + 2*nely/10 * nelx]=0
		img_list[9*nelx/10 + 9*nely/10 * nelx]=0
		return img_list
""""""

def main():
	"""
	Main to test the performance and robustness of the functions above before
	adding them to the topological optimization code of Martinez et al.
	"""
	from scipy import misc
	#Parallel flag :
	parallel_bool = False
	dummy_bool = True
	img_to_read = False
	img_path='test_image/'
#	img_name='frog'
#	img_format='_216x216.gif'
#	img_name='acropolisi'
	img_name='pascha'
	img_format='_64x64.png'
#	img_format='_32x32.png'
	extension_name = ''
	if img_to_read : extension_name += '_'+img_name +'_'
	if parallel_bool : extension_name += '_parallel_'

	#Construction of a benchmark image :
	if not img_to_read:
		(nelx, nely) = (50,25)
		img_list = construct_test_image(nelx,nely)
		img_array = img_list.reshape(nely,nelx)
		print '(nelx, nely) = {0}'.format(img_array.shape)
		scipy.misc.toimage(img_array, cmin=0.0, cmax=1).save('output/'+str(nelx) + 'x'+str(nely) + 'image_test.png')

	if img_to_read:
		#Using a test image in RGB format
		img_array = scipy.misc.imread(img_path+img_name+img_format)
		(nelx,nely) =[img_array.shape[i] for i in range(2)]
		print '(nelx, nely) = {0}'.format((nelx,nely))
		img_list = numpy.array(img_array).reshape(nelx*nely,4)

	#Symmetry planes :
	a1 = [-nely, nelx] # '\'plane
	c1 = [int(nelx/2), int(nely/2)]
	a2 = [nely,nelx]# '/'plane
	c2 = [nelx/2, nely/2]
	a3 = [ 1, 0]# '|'plane
	c3 = [int(nelx/2), int(nely/2)]
	a4 = [ 0, 1]# '-'plane
	c4 = [int(nelx/2), int(nely/2)]
	a_array = []
	c_array = []
#	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a1, c1)
#	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a2, c2)
	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a3, c3)
#	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a4, c4)
	start = time.time()
	if parallel_bool :
		connection_table = construct_connection_table_parallel(a_array, c_array, nelx,nely)
	elif dummy_bool :
		connection_table = construct_connection_table_dummy(a_array, c_array, nelx,nely)
#		connection_table = numpy.array(range(nelx*nely))
#		for n in range(len(a_array)):
#			a = a_array[n]
#			c = c_array[n]
#			connection_table = impose_symmetry_constraint(connection_table,a,c,nelx,nely)
	else:
		connection_table = construct_connection_table(a_array, c_array, nelx,nely)
	end = time.time()
	print('Ellapsed time : ' + str(end-start)+'[s]')

	plane_image = get_symmetry_image(a_array, c_array, nelx, nely)
	scipy.misc.toimage(plane_image, cmin=0.0, cmax=1).save('output/'+str(nelx) + 'x'+str(nely) + 'plane_image.png')

#	for i in range(nelx*nely):
#		if(connection_table[i]>-1):
#			img_list[i] = img_list[connection_table[i]]

	mapping_vector = construct_mapping_vector(connection_table)
	for i in range(len(mapping_vector)):
		sublist = mapping_vector[i]
		for index in sublist :
			img_list[index]=img_list[sublist[0]]
	"""
	for i in range(len(mapping_vector)):
		S=len(mapping_vector[i])
		mean = sum([img_list[j] for j in mapping_vector[i]])/S
		for index in mapping_vector[i]:
			img_list[index]=mean
	"""
	test= numpy.array(connection_table)
	print(test.reshape(nelx,nely))
	if img_to_read :
		scipy.misc.toimage(img_list.reshape(nely,nelx,4), cmin=0.0, cmax=256).save('output/'+str(nelx) + 'x'+str(nely) + extension_name + 'results.png')
	else:
		scipy.misc.toimage(img_list.reshape(nely,nelx), cmin=0.0, cmax=1).save('output/'+str(nelx) + 'x'+str(nely) + extension_name + 'results.png')

if __name__ == "__main__":
    main()
