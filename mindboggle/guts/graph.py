#!/usr/bin/env python
"""
Graph operations:

    - Diagonal degree matrix
    - Matrix weights and affinity matrix
    - Graph Laplacian

Authors:
    - Eliezer Stavsky, 2012 (eli.stavsky@gmail.com)
    - Arno Klein, 2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
import networkx as nx

from mindboggle.guts.kernels import rbf_kernel


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

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.graph import diagonal_degree_matrix
    >>> W = np.array([[10,2.3,3], [0,0,3], [0,1.5,0]])
    >>> tocsr = diagonal_degree_matrix(W, inverse=False, square_root=False)
    >>> tocsr.data
    array([ 15.3,   3. ,   1.5])

    """
    import numpy as np
    from scipy.sparse import lil_matrix

    stability_term = 0.000001

    ddm = lil_matrix((W.shape[0], W.shape[0]))

    if inverse:
        if not square_root:
            ddm.setdiag(1 / (W.sum(axis=1) + stability_term))
        else:
            #ddm.setdiag(math.sqrt(1 / (W.sum(axis=1) + stability_term)))
            ddm.setdiag(np.sqrt(1 / (W.sum(axis=1) + stability_term)))

    else:
        ddm.setdiag(W.sum(axis=1))

    return ddm.tocsr()


def weight_graph(Nodes, Indices, Meshes, kernel=rbf_kernel, add_to_graph=True,
                 G=nx.Graph(), sigma=20, verbose=False):
    """
    Construct weighted edges of a graph and compute an affinity matrix.

    Parameters
    ----------
    Nodes : numpy array
    Indices : list of integers
    Meshes : numpy array
    kernel : function which determines weights of edges
        - rbf_kernel: Gaussian kernel, with parameter sigma
        - cotangent_kernel: weight calculation for Laplace_Beltrami_Operator
          (NOTE: option removed until it can be tested)
        - inverse_distance: additional kernel where the weight is the inverse
          of the distance between two nodes
    add_to_graph :  boolean (add to graph?)
    G :  networkx graph
    sigma :  float (parameter for rbf_kernel)
    verbose : Boolean
        print statements?

    Returns
    -------
    G :  networkx graph
    affinity_matrix :  numpy array (sparse affinity matrix)

    Examples
    --------
    >>> import numpy as np
    >>> import networkx as nx
    >>> from mindboggle.guts.kernels import rbf_kernel
    >>> from mindboggle.guts.graph import weight_graph
    >>> Nodes = np.array([0,1,2,3,4])
    >>> Indices = [0,1,2,3,4]
    >>> Meshes = np.array([[1,2,3],[0,1,2],[0,1,3],[0,1,4],[0,2,3],[0,3,4]])
    >>> kernel = rbf_kernel
    >>> add_to_graph = True
    >>> G = nx.Graph()
    >>> sigma = 20
    >>> verbose = False
    >>> G, affinity_matrix = weight_graph(Nodes, Indices, Meshes, kernel,
    ...                                   add_to_graph, G, sigma, verbose)
    >>> G.size()
    9
    >>> G.degree()
    {0.0: 4, 1.0: 4, 2.0: 3, 3.0: 4, 4.0: 3}

    """
    import numpy as np
    from scipy.sparse import lil_matrix
    from mindboggle.guts.kernels import rbf_kernel, inverse_distance
                                        #cotangent_kernel

    if kernel is rbf_kernel or kernel is inverse_distance:
        if verbose:
            if kernel is rbf_kernel:
                print('Compute weights using rbf kernel (sigma={0})'.
                      format(sigma))
            else:
                print('Compute weights using inverse distance kernel '
                      '(sigma={0})'.format(sigma))

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
            if verbose:
                print('Add weighted edges to the graph')
            G.add_weighted_edges_from(weighted_edges)

        # Construct affinity matrix
        if verbose:
            print('Construct sparse affinity matrix of size {0}'.
                format(Nodes.shape[0]))
        affinity_matrix = lil_matrix((Nodes.shape[0], Nodes.shape[0]))
        for [i, j, edge_weight] in weighted_edges:
            affinity_matrix[i, j] = affinity_matrix[j, i] = edge_weight

    # elif kernel is cotangent_kernel:
    #     if verbose:
    #         print('Compute weights using cotangents')
    #     affinity_matrix = cotangent_kernel(Nodes, Meshes)
    #
    #     # Add weights to graph
    #     if add_to_graph:
    #         edges = np.nonzero(affinity_matrix)
    #         edge_mat = np.hstack((edges[0].T[:, np.newaxis],
    #                               edges[1].T[:, np.newaxis]))
    #         weighted_edges = np.asarray([[edge_mat[i,0],
    #                                       edge_mat[i,1],
    #                                       affinity_matrix[edge_mat[i]]]
    #                                       for i in range(affinity_matrix.shape[0])])
    #         if verbose:
    #             print('Add weighted edges to the graph')
    #         G.add_weighted_edges_from(weighted_edges)

    # Return the affinity matrix as a "compressed sparse row" matrix
    # (http://docs.scipy.org/doc/scipy/reference/sparse.html)
    if add_to_graph:
        return G, affinity_matrix.tocsr()
    else:
        return affinity_matrix.tocsr()


def graph_laplacian(W, type_of_laplacian='norm1', verbose=False):
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
    verbose : Boolean
        print statements?

    Returns
    -------
    Laplacian : N x N sparse matrix in csr format
               (Graph Laplacian of affinity matrix)

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.sparse as sparse
    >>> from mindboggle.guts.graph import graph_laplacian
    >>> row = np.array([0, 0, 1, 2, 2, 2])
    >>> col = np.array([0, 2, 2, 0, 1, 2])
    >>> data = np.array([1, 2, 3, 4, 5, 6])
    >>> W = sparse.csr_matrix((data, (row, col)), shape=(3, 3)).toarray()
    >>> W
    array([[1, 0, 2],
           [0, 0, 3],
           [4, 5, 6]])
    >>> type_of_laplacian = 'norm1'
    >>> verbose = False
    >>> Laplacian = graph_laplacian(W, type_of_laplacian, verbose)
    >>> print(np.array_str(np.array(Laplacian),
    ...       precision=5, suppress_small=True))
    [[ 0.66667  0.      -0.29814]
     [ 0.       1.      -0.44721]
     [-0.59628 -0.74536  0.6    ]]

    """

    if type_of_laplacian is 'basic':
        if verbose:
            print("Calculate unnormalized Laplacian")
        Laplacian = diagonal_degree_matrix(W) - W

    elif type_of_laplacian is 'norm1':
        if verbose:
            print("Normalize the Laplacian")
        ddmi_sq = diagonal_degree_matrix(W, inverse=True, square_root=True)
        Laplacian = ddmi_sq * (diagonal_degree_matrix(W, inverse=False, square_root=False) - W) * ddmi_sq

    elif type_of_laplacian is 'norm2':
        if verbose:
            print("Normalize the Laplacian")
        ddmi_sq = diagonal_degree_matrix(W, inverse=True, square_root=True)
        Laplacian = ddmi_sq * W * ddmi_sq

    elif type_of_laplacian is 'norm3':
        if verbose:
            print("Normalize the Laplacian")
        ddmi = diagonal_degree_matrix(W, inverse=True, square_root=False)
        Laplacian = ddmi * (diagonal_degree_matrix(W, inverse=False, square_root=False) - W)

    elif type_of_laplacian is 'random_walk':
        if verbose:
            print("Compute Random Walk Laplacian")
        ddmi = diagonal_degree_matrix(W, inverse=True, square_root=False)
        Laplacian = ddmi * W

    else:
        if verbose:
            print('Option is not available')
        Laplacian = 0

    return Laplacian


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()