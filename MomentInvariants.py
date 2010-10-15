#!/usr/bin/env python
# encoding: utf-8
"""
MomentInvariants.py
This scripts computes four moment invariants of the given data 
and returns order 3, 4, 5 and FL invariant moments as a numpy array.
It should be noted, current script is based on the assumption that 
the background intensity value is zero, otherwise a shape mask 
multiplication is required before running this script.  

Created by Ray Razlighi on 2010-10-11.
Copyright (c) 2010 __Columbia University__. All rights reserved.
"""

__version__ = "$Revision: 1 $"
# $Source$

import numpy as np


def moment_invarients(input_data, verbosity=0):
	
	input_data_size = input_data.shape
	temp_zero = np.zeros(input_data_size)
	#inpu_data_mask = np.array(np.array(input_data, dtype=bool), dtype=int)
	
	
	# overal summation of the input data
	M000 = input_data.sum()
	
	# The center point of the moment along x direction: mean_x
	x_coordinate_value = temp_zero.copy()
	for i in range(input_data_size[0]):  x_coordinate_value[i,:,:] = i
	M100 = x_coordinate_value * input_data
	M100 = M100.sum()
	mean_x = M100/M000
	
	# The center point of the moment along y direction: mean_y
	y_coordinate_value = temp_zero.copy()
	for i in range(input_data_size[1]):  y_coordinate_value[:,i,:] = i
	M010 = y_coordinate_value * input_data
	M010 = M010.sum()
	mean_y = M010/M000
	
	# The center point of the moment along z direction: mean_z
	z_coordinate_value = temp_zero.copy()
	for i in range(input_data_size[2]):  z_coordinate_value[:,:,i] = i
	M001 = z_coordinate_value * input_data
	M001 = M001.sum()
	mean_z = M001/M000
	
	mean = [M000, mean_x, mean_y, mean_z]
	
	# Compute the Central Moments (CM's) up to order 3 	
	x_coordinate_value = x_coordinate_value - mean[1]
	y_coordinate_value = y_coordinate_value - mean[2]
	z_coordinate_value = z_coordinate_value - mean[3]
	
	order = [2, 0, 0]
	CM200 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM200 = CM200.sum()
	
	order = [0, 2, 0]
	CM020 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM020 = CM020.sum()
	
	order = [0, 0, 2]
	CM002 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM002 = CM002.sum()
	
	order = [1, 1, 0]
	CM110 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM110 = CM110.sum()
	
	order = [1, 0, 1]
	CM101 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM101 = CM101.sum()
	
	order = [0, 1, 1]
	CM011 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM011 = CM011.sum()
	
	order = [1, 1, 1]
	CM111 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM111 = CM111.sum()
	
	order = [3, 0, 0]
	CM300 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM300 = CM300.sum()
	
	order = [0, 3, 0]
	CM030 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM030 = CM030.sum()
	
	order = [0, 0, 3]
	CM003 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM003 = CM003.sum()
	
	order = [0, 1, 2]
	CM012 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM012 = CM012.sum()
	
	order = [0, 2, 1]
	CM021 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM021 = CM021.sum()
		
	order = [1, 0, 2]
	CM102 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM102 = CM102.sum()
	
	order = [1, 2, 0]
	CM120 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM120 = CM120.sum()
	
	order = [2, 0, 1]
	CM201 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM201 = CM201.sum()
	
	order = [2, 1, 0]
	CM210 = pow(x_coordinate_value,order[0]) * pow(y_coordinate_value,\
	  order[1]) * pow(z_coordinate_value,order[2]) * input_data
	CM210 = CM210.sum()
	
	# Compute Rotational Invariant Moments of order 3, 4, 5 and FL
	invariant_moment_3 = (CM200 + CM020 + CM002) / pow(mean[0],5/3.0)
	
	invariant_moment_4 = (CM200*CM020 - pow(CM110,2) + CM200*CM002 - \
	  pow(CM101,2) + CM020*CM002 - pow(CM011,2)) / pow(mean[0],10/3.0)
	
	invariant_moment_5 = (CM200*CM020*CM002 + 2*CM110*CM101*CM011 - \
	  CM200*pow(CM011,2) - CM020*pow(CM101,2) - \
	  CM002*pow(CM110,2)) / pow(mean[0],15/3.0)
	
	invariant_moment_FL = (pow(CM003,2) + 6*pow(CM012,2) + \
	  6*pow(CM021,2) + pow(CM030,2) + 6*pow(CM102,2) + 15*pow(CM111,2) - \
	  3*CM102*CM120 + 6*pow(CM120,2) - 3*CM021*CM201 + 6*pow(CM201,2) - \
	  3*CM003*(CM021+CM201) - 3*CM030*CM210 + 6*pow(CM210,2) - \
	  3*CM012*(CM030+CM210) - 3*CM102*CM300 - 3*CM120*CM300 + \
	  pow(CM300,2)) / pow(mean[0],10/3.0)
	
	if verbosity:
		print('Rotational Invariant Moments of order 3: %g' %invariant_moment_3)
		print('Rotational Invariant Moments of order 4: %g' %invariant_moment_4)
		print('Rotational Invariant Moments of order 5: %g' %invariant_moment_5)
		print('Rotational FL Invariant Moments: %g' %invariant_moment_FL)

	return np.array([[invariant_moment_3], [invariant_moment_4], 
	  [invariant_moment_5], [invariant_moment_FL]])

