#!/usr/bin/env python
# encoding: utf-8
"""
AffineTransformationMatrix.py

This function optains 4X4 affine transformation matrix for the given
Translation, Rotation and scales. The order of the transformations are 
important. Assumed order is translation, rotation and then scaling. If
different order of transformation is required then this function is 
required to be called three times and multiply the result matrix in
the desired order. It is strongly recommended that this function to 
be used only for one form of transformation at the time.  
   
Created by Ray Razlighi on 2010-10-15.
Copyright (c) 2010 __Columbia University__. All rights reserved.
"""

from pylab import *
import numpy as np



def affine_transform_matrix(xyz_translation, rotation_angle, size_scales):
	
	x_translation = -xyz_translation[0]
	y_translation = -xyz_translation[1]
	z_translation = -xyz_translation[2]
	
	phi  = radians(rotation_angle[0])
	teta = radians(rotation_angle[1])
	psi  = radians(rotation_angle[2])
	
	x_scale = 1.0/size_scales[0]
	y_scale = 1.0/size_scales[1]
	z_scale = 1.0/size_scales[2]
	
	translation_matrix = np.matrix([
	[1,  0, 0, x_translation], 
	[0,  1, 0, y_translation], 
	[0,  0, 1, z_translation], 
	[0,  0, 0, 1            ]]) 
	
	rotation_matrix = np.matrix([
	[cos(teta)*cos(psi), -cos(phi)*sin(psi) + sin(phi)*sin(teta)*cos(psi),
	sin(phi)*sin(psi) + cos(phi)*sin(teta)*cos(psi), 0], 
	[cos(teta)*sin(psi),  cos(phi)*cos(psi) + sin(phi)*sin(teta)*sin(psi),
	-sin(phi)*cos(psi) + cos(phi)*sin(teta)*sin(psi), 0], 
	[sin(teta),           sin(phi)*cos(teta),                               
	cos(phi)*cos(teta),                              0],
	[0,                   0,                                                
	0,                                               1]])

	scale_matrix = np.matrix([
	[x_scale,  0,       0,        0], 
	[0,        y_scale, 0,        0], 
	[0,        0,       z_scale,  0], 
	[0,        0,       0,        1]]) 
		
	return rotation_matrix * scale_matrix * translation_matrix

