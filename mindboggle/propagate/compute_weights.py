""" Compute Weights - module for computing weights for affinity matrix.
------------------------------------------------
Input: Matrices containing nodes and meshes
Output: Sparse affinity matrix, (graph G)
------------------------------------------------
Parameters:
-----------
	Nodes: numpy array
	Meshes: numpy array
	kernel: function which determines weights of edges
	add_to_graph: boolean if weights should be added to graph
	G: networkx Graph
	sigma: float (parameter for rbf_kernel)

Features:
---------
rbf_kernel - Gaussian kernel, with parameter sigma
cotangent_kernel - weight calculation for Laplace_Beltrami_Operator
inverse_distance - additional kernel where the weight is the inverse of the disance between two nodes
--------------------------------------------------------------------
"""

import numpy as np
import networkx as nx
from scipy.sparse import lil_matrix
from kernels import *

# Default Parameters:
default_sigma = 20
default_graph = nx.Graph()
default_kernel = rbf_kernel

def compute_weights(Nodes, Meshes, kernel=default_kernel, add_to_graph=True,
					G=default_graph, sigma=default_sigma): 
	"""Constructs weighted edges of graph. Also computes affinity matrix for use later."""

	if kernel is rbf_kernel or kernel is inverse_distance:
		print 'Computing weights using specified kernel with parameter =', sigma
		
		# Construct matrix of edge lines by breaking triangle into three edges.
		if Meshes.shape[1] == 3:
			edge_mat = np.vstack((Meshes.T[0:2].T, Meshes.T[1:3].T, Meshes.T[:3:2].T))
		elif Meshes.shape[1] == 2:
			edge_mat = Meshes
		# Augment matrix to contain edge weight in the third column 
		weighted_edges = np.asarray([[i, j, kernel(Nodes[i], Nodes[j], sigma)] 
									  for [i, j] in edge_mat])
	
		# Add weights to graph
		if add_to_graph:
			print 'Adding weighted edges to the graph...'
			G.add_weighted_edges_from(weighted_edges)
	
		# Construct affinity matrix
		print 'Constructing sparse affinity matrix...'
		aff_mat = lil_matrix((Nodes.shape[0], Nodes.shape[0]))
		for [i, j, edge_weight] in weighted_edges:
			aff_mat[i, j] = aff_mat[j, i] = edge_weight
	
	elif kernel is cotangent_kernel:
		print 'Computing weights using cotangents...'
		aff_mat = cotangent_kernel(Nodes, Meshes)

		# Add weights to graph
		if add_to_graph:
			edges = np.nonzero(aff_mat)
			edge_mat = np.hstack((edges[0].T[:, np.newaxis], edges[1].T[:, np.newaxis]))
			weighted_edges = np.asarray([[edge_mat[i,0], edge_mat[i,1],aff_mat[edge_mat[i]]] for i in xrange(aff_mat.shape[0])])
			print 'Adding weighted edges to the graph...'
			G.add_weighted_edges_from(weighted_edges)			
	
	if add_to_graph:
		return G, aff_mat.tocsr()
	else:
		return aff_mat.tocsr()
