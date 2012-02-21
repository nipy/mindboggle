""" Laplace Beltrami Operator
--------------------------------------------------------
Parameters: file_name: file location
					VTK File containing nodes and meshes on which LBO will be performed
			which: string (Choice of which graph Laplacian.)
					basic - non-normalized Laplacian (Lap = D - W)
					norm1 - normalized Laplacian (Lap = ddmi_sq * L * ddmi_sq) - recovers definition
					norm2 - normalized Laplacian (Lap = ddmi_sq * W * ddmi_sq)
					norm3 - normalized Laplacian (Lap = inv(D) * L)
					random_walk - random walk Laplacian (Lap = inv(D) * W) 
			kernel: function (Choice of kernel to compute edge weights)
					rbf_kernel - Gaussian kernel with parameter sigma
					cotangent_kernel - typical weight calculation for Laplace_Beltrami_Operator
					inverse_distance - weight is the inverse of the disance between two nodes
			sigma:	float
					parameter for use in gaussian kernel
			components: float
					number of eigenvalues/vectors to compute
			export: boolean
					Would you like to write the eigenvalues to a file for independent analysis?
			return_header: boolean
					Would you like the header of the vtk file, extracted by pyvtk?
					
Returns:	Nodes: numpy array
			
			Meshes: numpy array 
			
			aff_mat: scipy sparse matrix in csr format
			
			Laplacian: scipy sparse matrix in csr format
			
			evals: numpy array of eigenvalues
			
			evecs: numpy array of eigenvectors
			
			header: string
---------------------------------------------------------			
"""
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import lil_matrix, csc_matrix
import pickle
from time import time
import pyvtk

import vtk_operations as vo
from compute_weights import *
from graph_operations import *

default_VTK = '/home/eli/PythonFiles/Data/testdatalabels.vtk'

################################################################## 
################### Main Function: ##############################

def Laplace_Beltrami_Operator(file_name=default_VTK, which='norm1', kernel=cotangent_kernel, sigma=10, components=200, export=True, return_header=False):
	""" Laplace_Beltrami_Operator."""
	
	##############################################################
	##### Extract Nodes and Meshes from VTK file: 
	
	Data = pyvtk.VtkData(file_name)
	Nodes = np.asarray(Data.structure.points)
	Meshes = np.asarray(Data.structure.polygons)
	
	##############################################################
	##### Compute Weight Matrix:
	
	t0 = time()
	aff_mat = compute_weights(Nodes, Meshes, kernel=kernel, add_to_graph=False, sigma=sigma)
	print("done in %0.3fs" % (time() - t0))

	##############################################################
	##### Compute Graph Laplacian:

	t0 = time()
	Laplacian = graph_laplacian(aff_mat, which=which)
	print("done in %0.3fs" % (time() - t0))

	##############################################################
	##### Obtain Eigen Spectrum:
	
	print 'Finding the eigen spectrum...'
	t0 = time()
	Laplacian = Laplacian.tocsc()
	if components >= Nodes.shape[0]:
		components = Nodes.shape[0] - 1
	evals, evecs = eigsh(Laplacian, components, sigma=0)
	print("done in %0.3fs" % (time() - t0))
	
	if export:
		file_name2 = file_name[:-4] + '.npy'
		h = open(file_name2, 'w')
		pickle.dump(evals, h)
		h.close()

	################################################################
	##### Write New VTK File With Feidler Vector Overlayed On Brain:
	
	# f = vo.write_all(file_name[:-4]+'_eigen.vtk', Nodes, Meshes, evecs[:,1], label_type='Eigenvector')

	#################################################################
	##### Return Results:
	
	if return_header:
		return Nodes, Meshes, aff_mat, Laplacian, evals, evecs, Data.header
	else:
		return Nodes, Meshes, aff_mat, Laplacian, evals, evecs
