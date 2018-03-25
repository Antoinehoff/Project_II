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
	return numpy.array([index % nelx+1/2, int(index/nelx)+1/2])
	
def get_symmetry_plane(nelx,nely):
	a = numpy.array([0, 1])
	a = a/(numpy.sqrt(numpy.dot(a,a)))
	c = [int(nelx/2)-1, int(nely/2)-1]
	return (a, c)

def main():
	nelx = 10
	nely = 10
	
	x = numpy.arange(0,nelx*nely,1)
	(a,c) = get_symmetry_plane(nelx,nely)
	Xmin = index_to_position(0, nelx, nely)
	print("Xmin = ", Xmin)
	Xmax = index_to_position(nelx*nely-1, nelx, nely)
	print("Xmax = ", Xmax)
	rmax = numpy.sqrt(numpy.dot(Xmax-Xmin,Xmax-Xmin))
	print("rmax = ",rmax)

	epsilon_1 = 0.00001
	epsilon_2 = 0.001
	nmast = 0
	k=0
	alpha_ij = 0
	image = numpy.arange(nelx*nely).reshape(nelx,nely)
	
	Connection_table = numpy.ones(nelx*nely)

	for i in range(nelx*nely):
		X_i = index_to_position(i, nelx, nely)
		if(numpy.dot(X_i-c,a) < epsilon_1*rmax):
		#	image[index_to_position(i,nelx,nely)[0]][index_to_position(i,nelx,nely)[1]]=1
			k=k+1
			X_i_proj = X_i - (numpy.dot(X_i-c,a))*a
		#	image[X_proj_index[0]][X_proj_index[1]]=0
		#else:
		#	image[index_to_position(i,nelx,nely)[0]][index_to_position(i,nelx,nely)[1]]=0.5
		dmin = rmax
		
		for j in range(nelx*nely):
			X_j = index_to_position(j, nelx, nely)
			if(i!=j and numpy.dot(X_j-c,a) > epsilon_1*rmax):
				temp1 = numpy.dot(X_i_proj - X_i, X_i_proj - X_j)
				temp2 = numpy.linalg.norm(X_i_proj - X_i)*numpy.linalg.norm(X_i_proj - X_j)+ epsilon_2
				alpha_ij = temp1 / temp2
				if(alpha_ij > 1-epsilon_2):
					if(abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))<dmin):
						nmast=j
						dmin=abs(numpy.linalg.norm(X_i_proj - X_i)-numpy.linalg.norm(X_i_proj - X_j))
						Connection_table[i]=nmast	
	
	print(Connection_table.reshape(nelx,nely))					
	print(k)
	scipy.misc.toimage(image, cmin=0.0, cmax=nelx*nely).save("output.png")


if __name__ == "__main__":
    main()
