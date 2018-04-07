#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test script for additional functions
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


def design_variable_merge(x, nelx, nely):
	"""
	Implementation of the algorithm from Kosaka and Swan 1999
	"""			
	return x
	
def index_to_position(index, nelx, nely):
	"""
	Convert the index of a element to the centroid of the element
	"""
	return numpy.array([(index % nelx)+.5, int(index/nelx)+.5])

def add_symmetry_planes(a_array, c_array, a, c, empty = 0):
	if empty :
		a_array = [numpy.array(a/numpy.sqrt(numpy.dot(a,a)))]
		c_array = [numpy.array(c)]
	else :
		a = a/numpy.sqrt(numpy.dot(a,a))
		a_array = numpy.append(a_array,[a],0)
		c_array = numpy.append(c_array,[c],0)
	return a_array, c_array
	
def get_symmetry_image(a_array, c_array, nelx, nely):
	image = numpy.ones(nelx*nely)
	
	for n in range(numpy.shape(a_array)[0]):
		a = a_array[n]
		c = c_array[n]
		for i in range(nelx*nely):
			X_i = index_to_position(i, nelx, nely)
			if(numpy.dot(X_i-c,a) > 0):
				image[i] = .5;
	return (image.reshape(nely,nelx))
	
def in_domain(X, a_array, c_array, epsilon, rmax):
	#Check if the point X is in the domain defined by the planes
	# Interior is for negative dot product
	flag = True
	for n in range(numpy.shape(a_array)[0]): 
		a = a_array[n]
		c = c_array[n]
		if(numpy.dot(X-c,a) > epsilon*rmax):
			flag = False
	
	return flag

def construct_connection_table(a_array, c_array, nelx, nely):	
	#Initialization phase
	Xmin = index_to_position(0, nelx, nely)
	print 'Xmin = {0}'.format(Xmin)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	print 'Xmax = {0}'.format(Xmax)
	rmax = numpy.sqrt(numpy.dot(Xmax-Xmin,Xmax-Xmin)) #Length scale of mesh
	print 'rmax = {0}'.format(rmax)
	epsilon_1 = 0.0000001
	epsilon_2 = 0.0000001
	nmast = 0
	k=0 #Counter on elements in the interior of the plane
	m=0
	alpha_ij = 0
	connection_table = numpy.ones(nelx*nely)*-1
		
	for n in range(numpy.shape(a_array)[0]): #For each symmetry plane...
		a = a_array[n]
		c = c_array[n]
		for i in range(nelx*nely): #1st Loop over all elements
			X_i = index_to_position(i, nelx, nely) #Centroid of element
			
			if(numpy.dot(X_i-c,a) < rmax * epsilon_1):
#			if(in_domain(X_i, a_array, c_array, epsilon_1, rmax)):
				for n in range(numpy.shape(a_array)[0]):
					a = a_array[n]
					c = c_array[n]
					k=k+1
					X_i_proj = X_i - (numpy.dot(X_i-c,a))*a
					dmin = rmax		
					for j in range(nelx*nely):
						X_j = index_to_position(j, nelx, nely)
						if(i!=j and not in_domain(X_j, a_array, c_array, epsilon_1, rmax)):
							temp1 = numpy.dot(X_i_proj - X_i, X_i_proj - X_j)
							temp2 = numpy.linalg.norm(X_i_proj - X_i)*numpy.linalg.norm(X_i_proj - X_j) + epsilon_2
							alpha_ij = temp1 / temp2
		#					print("i = ",i," j = ",j," temp1 = ",temp1, " temp2 = ", temp2, "alpha_ij = ", alpha_ij)
							if(alpha_ij < epsilon_2-1): #Not the same as in Kosaka and Swan paper
								if(abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))<dmin):
									nmast=j
									dmin=abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))
									m=m+1
		#							print("Replacing coordinate ",i," by ",nmast," (j=",j,")")
									connection_table[nmast]=i	
	print '#elements : {0} '.format(nelx*nely)
	print '#inside elements : {0}({1}%)'.format(k,int(100.0*k/nelx/nely))
	print '#dependant elements : {0}({1}%)'.format(m,int(100.0*m/nelx/nely))
	return connection_table.astype(int)

def main():
	from scipy import misc
	#Construction of a benchmark image :
	(nelx, nely) = (64,64)
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
	img_array = img_list.reshape(nely,nelx)
	print '(nelx, nely) = {0}'.format(img_array.shape)
	scipy.misc.toimage(img_array, cmin=0.0, cmax=1).save('output/'+str(nelx) + 'x'+str(nely) + 'image_test.png')
	
	#Symmetry planes :
	a1 = [-nely, nelx]
	c1 = [int(nelx/2), int(nely/2)]
	a2 = [nely,nelx]
	c2 = [int(nelx/2), int(nely/2)]
	a3 = []
	a_array = []
	c_array = []
	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a2, c2,1)
	(a_array, c_array) = add_symmetry_planes(a_array, c_array, a1, c1)
	
	connection_table = construct_connection_table(a_array, c_array, nelx,nely)
	
	plane_image = get_symmetry_image(a_array, c_array, nelx, nely)
	scipy.misc.toimage(plane_image, cmin=0.0, cmax=1).save('output/'+str(nelx) + 'x'+str(nely) + 'plane_image.png')
	
	#print(connection_table)
	
	for i in range(nelx*nely):
		if(connection_table[i]>-1):
			img_list[i] = img_list[connection_table[i]]
	
	scipy.misc.toimage(img_list.reshape(nely,nelx), cmin=0.0, cmax=1).save('output/'+str(nelx) + 'x'+str(nely) + 'results.png')


if __name__ == "__main__":
    main()
