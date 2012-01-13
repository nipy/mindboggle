""" Laplace Beltrami Operator
--------------------------------------------------------
Parameters: VTK: file location
					File containing nodes and meshes on which LBO will be performed
			which: string (Choice of which graph Laplacian.)
					basic - non-normalized Laplacian (Lap = D - W)
					norm1 - normalized Laplacian (Lap = ddmi_sq * L * ddmi_sq) - recovers definition
					norm2 - normalized Laplacian (Lap = ddmi_sq * W * ddmi_sq)
					norm3 - normalized Laplacian (Lap = inv(D) * L)
					random_walk - random walk Laplacian (Lap = inv(D) * W) 
			kernel: function (Choice of kernel to compute edge weights)
					rbf_kernel - Gaussian kernel with parameter sigma
					cotangent_kernel - typical weight calculation for laplace_beltrami_Operator
					inverse_distance - weight is the inverse of the disance between two nodes
			sigma:	float
					parameter for use in gaussian kernel
			components: float
					number of eigenvalues/vectors to compute
			export: boolean
					Would you like to write the eigenvalues to a file for independent analysis?

Returns:	Nodes: numpy array
			
			Meshes: numpy array 
			
			aff_mat: scipy sparse matrix in csr format
			
			Laplacian: scipy sparse matrix in csr format
			
			evals: array of eigenvalues
			
			evecs: array of eigenvectors
			
			vtk files
---------------------------------------------------------			
"""
import numpy as np
from scipy.sparse.linalg.eigen.arpack import *
from scipy.sparse import lil_matrix, csc_matrix
import pickle

from vtk_operations import *
from compute_weights import *
from graph_operations import *

default_VTK = '/Users/eli/Documents/PythonFiles/Data/testdatalabels.txt'
default_kernel = cotangent_kernel

################################################################## 
### Main Function:
def laplace_beltrami_operator(VTK=default_VTK, which='norm1', 
							  kernel=default_kernel, sigma=10, components=50, export=False):
	""" laplace_beltrami_operator."""
	
	# Extract Nodes and Meshes from VTK file:
	vtk = convert_vtk_to_matrix(VTK, return_files=True, labels=False, fundi=False, delete=False)
	Nodes, Meshes, A = vtk['Nodes'], vtk['Meshes'], vtk['main_file']
	
	# Compute Weight Matrix:
	aff_mat = compute_weights(Nodes, Meshes, kernel=kernel, add_to_graph=False, sigma=sigma)
	
	# Compute Graph Laplacian:
	Laplacian = graph_laplacian(aff_mat, which=which)
	
	# Obtain Eigen Spectrum:
	print 'Finding the eigen spectrum...'
	Laplacian = Laplacian.tocsc()
	if components >= Nodes.shape[0]:
		components = Nodes.shape[0] - 1
	evals, evecs = eigsh(Laplacian, components, which='SM') # does sigma=0 work?
	
	if export:
		file_name = VTK[:-4] + '.npy'
		h = open(file_name, 'w')
		pickle.dump(evals, h)
		h.close()

	return Nodes, Meshes, aff_mat, Laplacian, evals, evecs
	
	
	
## Piece of Old Code:
	'''
	# Linearize to visualize.
	
	ones = np.ones(Nodes.shape[0])
	natural = np.cumsum(ones)
	print 'Writing results to a file...'

	for i in xrange(1, components):
		vector = np.argsort(evecs[:,i])
		components = np.zeros(Nodes.shape[0])
		components[vector] = natural
		file_name = '/Users/eli/Desktop/Eigen' + str(i) + '.vtk'
		h = open(file_name, 'w')
		line_num = 0
		total_lines_to_copy = 9 + Nodes.shape[0] + Meshes.shape[0] 
		for line in A:
			if line_num < total_lines_to_copy:
				h.write(line)
				line_num += 1
		for i in xrange(len(evecs[:,i])):
			h.write('{0}\n'.format(str(components[i])))
	
		h.close()
	'''
