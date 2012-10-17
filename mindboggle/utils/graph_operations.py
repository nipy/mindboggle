#!/usr/bin/python
"""
Graph operations:

    - Diagonal degree matrix
    - Matrix weights and affinity matrix
    - Graph Laplacian

Authors:
    - Eliezer Stavsky (eli.stavsky@gmail.com)

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import numpy as np
import networkx as nx
#from scipy.sparse import lil_matrix
from mindboggle.utils.kernels import rbf_kernel, cotangent_kernel, inverse_distance

###############################################################################
# -----------------------------------------------------------------------------
#     Diagonal degree matrix
# -----------------------------------------------------------------------------
###############################################################################

def diagonal_degree_matrix(W, inverse=False, square_root=False):
    """
    Compute diagonal degree matrix.

    Parameters
    ----------
    W : N x N sparse matrix in csr format (affinity matrix)
    inverse : boolean (compute inverse of diagonal degree matrix?)
    square_root : boolean (compute square root of diagonal degree matrix?)

    Returns
    -------
    ddm : N x N sparse matrix in csr format (diagonal matrix)
         "csr" stands for "compressed sparse row" matrix
         (http://docs.scipy.org/doc/scipy/reference/sparse.html)
    """
    import math
    from scipy.sparse import lil_matrix

    stability_term = 0.000001

    ddm = lil_matrix((W.shape[0], W.shape[0]))

    if inverse:
        if not square_root:
            ddm.setdiag(1 / (W.sum(axis=1) + stability_term))
        else:
            ddm.setdiag(math.sqrt(1 / (W.sum(axis=1) + stability_term)))

    else:
        ddm.setdiag(W.sum(axis=1))

    return ddm.tocsr()

###############################################################################
# -----------------------------------------------------------------------------
#     Matrix weights and affinity matrix
# -----------------------------------------------------------------------------
###############################################################################

def weight_graph(Nodes, Indices, Meshes, kernel=rbf_kernel, add_to_graph=True,
                 G=nx.Graph(), sigma=20):
    """
    Construct weighted edges of a graph and compute an affinity matrix.

    Parameters
    ----------
    Nodes : numpy array
    Indices : numpy array
    Meshes : numpy array
    kernel : function which determines weights of edges
        - rbf_kernel: Gaussian kernel, with parameter sigma
        - cotangent_kernel: weight calculation for Laplace_Beltrami_Operator
        - inverse_distance: additional kernel where the weight is the inverse
          of the distance between two nodes
    add_to_graph :  boolean (add to graph?)
    G :  networkx graph
    sigma :  float (parameter for rbf_kernel)

    Returns
    -------
    G :  networkx graph
    affinity_matrix :  numpy array (sparse affinity matrix)

    """
    import numpy as np
    from scipy.sparse import lil_matrix
    from mindboggle.utils.kernels import rbf_kernel, cotangent_kernel, inverse_distance

    if kernel is rbf_kernel or kernel is inverse_distance:
        print('Computing weights using {0} kernel with parameter = {1}'.format(
              kernel, sigma))

        # Construct matrix of edge lines by breaking triangle into three edges.
        if Meshes.shape[1] == 3:
            edge_mat = np.vstack((Meshes.T[0:2].T, Meshes.T[1:3].T, Meshes.T[:3:2].T))
        elif Meshes.shape[1] == 2:
            edge_mat = Meshes
        # Augment matrix to contain edge weight in the third column
        weighted_edges = np.asarray([[Indices[i], Indices[j],
            kernel(Nodes[Indices[i]], Nodes[Indices[j]], sigma)]
            for [i, j] in edge_mat])

        # Add weights to graph
        if add_to_graph:
            print('Adding weighted edges to the graph...')
            G.add_weighted_edges_from(weighted_edges)

        # Construct affinity matrix
        print('Constructing sparse affinity matrix...')
        affinity_matrix = lil_matrix((Nodes.shape[0], Nodes.shape[0]))
        for [i, j, edge_weight] in weighted_edges:
            affinity_matrix[i, j] = affinity_matrix[j, i] = edge_weight

    elif kernel is cotangent_kernel:
        print('Computing weights using cotangents...')
        affinity_matrix = cotangent_kernel(Nodes, Meshes)

        # Add weights to graph
        if add_to_graph:
            edges = np.nonzero(affinity_matrix)
            edge_mat = np.hstack((edges[0].T[:, np.newaxis],
                                  edges[1].T[:, np.newaxis]))
            weighted_edges = np.asarray([[edge_mat[i,0],
                                          edge_mat[i,1],
                                          affinity_matrix[edge_mat[i]]]
                                          for i in xrange(affinity_matrix.shape[0])])
            print('Adding weighted edges to the graph...')
            G.add_weighted_edges_from(weighted_edges)

    # Return the affinity matrix as a "compressed sparse row" matrix
    # (http://docs.scipy.org/doc/scipy/reference/sparse.html)
    if add_to_graph:
        return G, affinity_matrix.tocsr()
    else:
        return affinity_matrix.tocsr()

###############################################################################
# -----------------------------------------------------------------------------
#     Graph Laplacian
# -----------------------------------------------------------------------------
###############################################################################

def graph_laplacian(W, type_of_laplacian='norm1'):
    """
    Compute normalized and unnormalized graph Laplacians.

    Parameters
    ----------
    W : N x N sparse matrix in csr format (affinity matrix)
       "csr" stands for "compressed sparse row" matrix
       (http://docs.scipy.org/doc/scipy/reference/sparse.html)
    type_of_laplacian : string
        - basic: non-normalized Laplacian (Lap = D - W)
        - norm1: normalized Laplacian (Lap = ddmi_sq * L * ddmi_sq) - recovers definition
        - norm2: normalized Laplacian (Lap = ddmi_sq * W * ddmi_sq)
        - norm3: normalized Laplacian (Lap = inv(D) * L)
        - random_walk: random walk Laplacian (Lap = inv(D) * W)

    Returns
    -------
    Laplacian : N x N sparse matrix in csr format
               (Graph Laplacian of affinity matrix)

    """

    if type_of_laplacian is 'basic':
        print 'Calculating unnormalized Laplacian...'
        Laplacian = diagonal_degree_matrix(W) - W
        return Laplacian

    elif type_of_laplacian is 'norm1':
        print "Normalizing the Laplacian..."
        ddmi_sq = diagonal_degree_matrix(W, inverse=True, square_root=True)
        Laplacian = ddmi_sq * (diagonal_degree_matrix(W, inverse=False, square_root=False) - W) * ddmi_sq
        return Laplacian

    elif type_of_laplacian is 'norm2':
        print "Normalizing the Laplacian..."
        ddmi_sq = diagonal_degree_matrix(W, inverse=True, square_root=True)
        Laplacian = ddmi_sq * W * ddmi_sq
        return Laplacian

    elif type_of_laplacian is 'norm3':
        print "Normalizing the Laplacian..."
        ddmi = diagonal_degree_matrix(W, inverse=True, square_root=False)
        Laplacian = ddmi * (diagonal_degree_matrix(W, inverse=False, square_root=False) - W)
        return Laplacian

    elif type_of_laplacian is 'random_walk':
        print "Computing Random Walk Laplacian..."
        ddmi = diagonal_degree_matrix(W, inverse=True, square_root=False)
        Laplacian = ddmi * W

    else:
        print 'Option is not available'
        return 0
