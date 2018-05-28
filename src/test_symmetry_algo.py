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
import time

import connection_table as ct

def index_to_position(index, nelx, nely):
	"""
	Convert the index of a element to the centroid of the element
	"""
	pos = numpy.array([(index / (nely))+.5, int(index%(nely))+.5])
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
	Return the index of the symmetric element w.r.t. a plane (a,c)
	"""
	x_i = index_to_position(index,nelx,nely)
	dist = numpy.dot(x_i-c,a)
	x_proj = x_i - dist * a
	x_sym = x_proj - dist * a
	index_sym = position_to_index(x_sym,nelx,nely)
	return index_sym

def construct_connection_table(a_array, c_array, nelx, nely):
	"""
	Simple algorithm O(nelx*nely) to construct the table containing
	for an element i, outside the design domain, its symmetric
	element j inside the design domain. connection_table[i]=j if i is outside the
	domain, =-1 if it is inside the domain or have no symmetric element.
	"""
	connection_table = numpy.array(range(nelx*nely))
	for n in range(numpy.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		for i in range(nelx*nely):
			X_i = index_to_position(i,nelx,nely)
			index_sym = connection_table[i]
			if in_domain(X_i, [a], [c])!=1:
				index_sym = connection_table[get_symmetric_element(i,a,c,nelx,nely)]
			connection_table[i]=index_sym
	return connection_table

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


def main():
	"""
	Main to test the performance and robustness of the functions above before
	adding them to the topological optimization code of Martinez et al.
	"""
	from scipy import misc
	img_path='test_image/'
#	img_name='frog'
#	img_format='_216x216.gif'
#	img_name='acropolisi'
	img_name='pascha'
	img_format='_64x64.png'
#	img_format='_32x32.png'
	#Using a test image in RGB format
	img_array = scipy.misc.imread(img_path+img_name+img_format)
	(nelx,nely) =[img_array.shape[i] for i in range(2)]
	img_list = numpy.array(img_array).reshape(nelx*nely,4)
	for i in range(nelx):
		for j in range(nely):
			img_list[nely*j + i]=img_array[i][j]
	extension_name=''
	#Symmetry planes :
	a1 = [-nely, nelx] # '\'plane
	c1 = [int(nelx/2), int(nely/2)]
	a2 = [nely,nelx]# '/'plane
	c2 = [nelx/2, nely/2]
	a3 = [ 0, 1]# '-'plane
	c3 = [int(nelx/2), int(nely/2)]
	a4 = [ 1, 0]# '|'plane
	c4 = [int(nelx/2), int(nely/2)]
	a_array = []
	c_array = []
#	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a1, c1)
	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a2, c2)
#	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a3, c3)
#	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a4, c4)
	start = time.time()
	connection_table = construct_connection_table(a_array, c_array, nelx,nely)
	end = time.time()
	print('Ellapsed time : ' + str(end-start)+'[s]')

#	plane_image = get_symmetry_image(a_array, c_array, nelx, nely)
#	scipy.misc.toimage(plane_image, cmin=0.0, cmax=1).save('output/'+str(nelx) + 'x'+str(nely) + 'plane_image.png')

#	print(connection_table.reshape(nelx,nely))
	mapping_vector = ct.construct_mapping_vector_wheel(3,nelx,nely)
	print(mapping_vector)
	for sublist in mapping_vector:
		temp = img_list[sublist[0]]
		for index in sublist :
			img_list[index]=temp

	scipy.misc.toimage(img_list.reshape(nely,nelx,4), cmin=0.0, cmax=256).save('output/'+str(nelx) + 'x'+str(nely) + extension_name + 'results.png')
	print('ouput saved at : ' +'output/'+str(nelx) + 'x'+str(nely) + extension_name + 'results.png')
if __name__ == "__main__":
    main()
