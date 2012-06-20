###		Graph Operations:

import numpy as np
from scipy import sparse

###############################################################################
# Diagonal Degree Matrix

"""  Function: Computes Diagonal Degree Matrix:
------------------------------------------------------
Parameters: W: N x N sparse matrix in csr format
				Affinity matrix
			inverse: boolean
				Compute inverse of diagonal degree matrix?
			square_root: boolean
				Compute square root of diagonal degree matrix?

Returns:	ddm: N x N sparse matrix in csr format
				Diagonal matrix.
----------------------------------------------------  """

def compute_diagonal_degree_matrix(W, inverse=False, square_root=False):
	
	ddm = sparse.lil_matrix((W.shape[0], W.shape[0]))
	
	if inverse:
		if not square_root:
			ddm.setdiag(1 / W.sum(axis=1))	
		else:
			ddm.setdiag(np.sqrt(1 / W.sum(axis=1)))
	
	else:
		ddm.setdiag(W.sum(axis=1))		
	
	return ddm.tocsr()
	
###############################################################################
# Graph laplacian

""" Function: Computes normalized and unnormalized graph laplacians
----------------------------------------------------------
Parameters: W: N x N sparse matrix
				Matrix in csr format, affinity matrix
			which: string
				basic - non-normalized Laplacian (Lap = D - W)
				norm1 - normalized Laplacian (Lap = ddmi_sq * L * ddmi_sq) - recovers definition
				norm2 - normalized Laplacian (Lap = ddmi_sq * W * ddmi_sq)
				norm3 - normalized Laplacian (Lap = inv(D) * L)
				random_walk - random walk Laplacian (Lap = inv(D) * W) 

Returns:	Laplacian: N x N sparse matrix
				Matrix in csr format, Graph Laplacian of affinity matrix 

---------------------------------------------------------- """

def graph_laplacian(W, which='norm1'):
	
	if which is 'basic':
		print 'Calculating Unnormalized Laplacian...'
		Laplacian = compute_diagonal_degree_matrix(W) - W	
		return Laplacian
	
	elif which is 'norm1':
		print "Normalizing the Laplacian..."
		ddmi_sq = compute_diagonal_degree_matrix(W, inverse=True, square_root=True)
		Laplacian = ddmi_sq * (compute_diagonal_degree_matrix(W, inverse=False, square_root=False) - W) * ddmi_sq
		return Laplacian
		
	elif which is 'norm2':
		print "Normalizing the Laplacian..."
		ddmi_sq = compute_diagonal_degree_matrix(W, inverse=True, square_root=True)
		Laplacian = ddmi_sq * W * ddmi_sq
		return Laplacian
		
	elif which is 'norm3':
		print "Normalizing the Laplacian..."
		ddmi = compute_diagonal_degree_matrix(W, inverse=True, square_root=False)
		Laplacian = ddmi * (compute_diagonal_degree_matrix(W, inverse=False, square_root=False) - W)
		return Laplacian
		
	elif which is 'random_walk':
		print "Computing Random Walk Laplacian..."
		ddmi = compute_diagonal_degree_matrix(W, inverse=True, square_root=False)
		Laplacian = ddmi * W
	
	else:
		print 'Option is not available'
		return 0


		