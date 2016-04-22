#!/usr/bin/env python
"""
Kernels.

Authors:
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)
    - Arno Klein, 2016  (arno@mindboggle.info)

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def rbf_kernel(x1, x2, sigma):
    """
    Compute normalized and unnormalized graph Laplacians.

    Parameters
    ----------
    x1 : Nx1 numpy array
    x2 : Nx1 numpy array
    sigma : float

    Returns
    -------
    rbf : float

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.kernels import rbf_kernel
    >>> x1 = np.array([0.1,0.2,0.4,0])
    >>> x2 = np.array([0.1,0.3,0.5,0])
    >>> sigma = 0.5
    >>> rbf = rbf_kernel(x1, x2, sigma)
    >>> print('{0:0.5f}'.format(rbf))
    0.96079

    """
    import numpy as np

    return np.exp(-np.linalg.norm(x1 - x2) ** 2 / (2 * sigma ** 2))


# def cotangent_kernel(Nodes, Meshes):
#     """
#     This function constructs weighted edges of a graph.
#
#     Parameters
#     ----------
#     Nodes : numpy array
#     Meshes : numpy array
#
#     Returns
#     -------
#     W : N x N matrix
#         weight matrix
#
#     Examples
#     --------
#     >>> import numpy as np
#     >>> from mindboggle.guts.kernels import cotangent_kernel
#     >>> Nodes = np.array([0,1,2,3,4])
#     >>> Meshes = np.array([[1,2,3],[0,1,2],[0,1,3],[0,1,4],[0,2,3],[0,3,4]])
#     >>> cotangent_kernel(Nodes, Meshes)
#     ValueError: 'axisa' out of bounds
#
#     """
#     import numpy as np
#     from scipy.sparse import lil_matrix
#
#     num_nodes = Nodes.shape[0]
#     W = lil_matrix((num_nodes, num_nodes))
#     #print('Constructing sparse affinity matrix...')
#     for c in Meshes:
#         # Obtain vertices which comprise face
#         v0, v1, v2 = Nodes[c[0]], Nodes[c[1]], Nodes[c[2]]
#
#         # Obtain cotangents of angles
#         cot0 = np.dot(v1-v0, v2-v0) / np.linalg.norm(np.cross(v1-v0, v2-v0))
#         cot1 = np.dot(v2-v1, v0-v1) / np.linalg.norm(np.cross(v2-v1, v0-v1))
#         cot2 = np.dot(v0-v2, v1-v2) / np.linalg.norm(np.cross(v0-v2, v1-v2))
#
#         # Update weight matrix accordingly
#         W[c[1], c[2]] += cot0
#         W[c[2], c[1]] += cot0
#         W[c[0], c[2]] += cot1
#         W[c[2], c[0]] += cot1
#         W[c[0], c[1]] += cot2
#         W[c[1], c[0]] += cot2
#
#     return W


def inverse_distance(x1, x2, epsilon):
    """
    This function constructs weighted edges of a graph,
    where the weight is the inverse of the distance between two nodes.

    Parameters
    ----------
    x1 : Nx1 numpy array
    x2 : Nx1 numpy array
    epsilon : float

    Returns
    -------
    d : float

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.kernels import inverse_distance
    >>> x1 = np.array([0.1,0.2,0.4,0])
    >>> x2 = np.array([0.1,0.3,0.5,0])
    >>> epsilon = 0.05
    >>> d = inverse_distance(x1, x2, epsilon)
    >>> print('{0:0.5f}'.format(d))
    5.22408

    """
    import numpy as np

    return 1.0/(np.linalg.norm(x1 - x2) + epsilon)


# ============================================================================
# Doctests
# ============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules