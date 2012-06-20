#!/usr/bin/env python
# encoding: utf-8
"""
MomentInvariants_Test.py
This script tests the computation of 4 moment invariants by function:
moment_invariants()

Created by Ray Razlighi on 2010-10-11.
Copyright (c) 2010 __Columbia University__. All rights reserved.
"""

__version__ = "$Revision: 1 $"
# $Source$

import sys
import numpy as np
from pylab import *
import scipy.ndimage.interpolation as spip
import mpl_toolkits.mplot3d.axes3d as p3


from MomentInvariants import moment_invarients
from AffineTransformMatrix import affine_transform_matrix



def main(argv=None):
	
	# create a random image volum inside a zero space
	input_image = np.zeros([160,160,10])
	internal_shape = np.round(np.random.random([100,100,10])*100)
	input_image[30:130,30:130, :] = internal_shape
	
	# create affine transformation matrix for x_translation(20), z_rotation(10) 
	# and x_scaling(.5)
	affine_matrix = affine_transform_matrix([20, 0, 0], [0, 0, 10], [.5, 1, 1])
	
	# obtains the transformed image volume
	output_image = np.round(spip.affine_transform(input_image, 
	  affine_matrix[0:3,0:3], offset=[-20, 0, 0], output_shape=None, 
	  order=3, mode='constant', cval=0.0, prefilter=True))
	
	# computes and prints the moments
	moment_invarients(input_image, 1)
	moment_invarients(output_image, 1)
	
	#visualization of the two image volume
	(x, y, z) = input_image.nonzero()
	fig=figure(1)
	ax = p3.Axes3D(fig)
	ax.scatter3D(ravel(x),ravel(y),ravel(z))
	ax.set_xlim3d([0,160])
	ax.set_ylim3d(0,160)
	ax.set_zlim3d(0,10)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	fig.add_axes(ax)
	
	(x, y, z) = output_image.nonzero()
	fig=figure(2)
	ax = p3.Axes3D(fig)
	ax.scatter3D(ravel(x),ravel(y),ravel(z))
	ax.set_xlim3d([0,160])
	ax.set_ylim3d(0,160)
	ax.set_zlim3d(0,10)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	fig.add_axes(ax)
	show()
	
	
	
if __name__ == '__main__':
	main()





