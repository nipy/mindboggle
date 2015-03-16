#!/usr/bin/env python
"""
Kernels.

Authors:
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def rbf_kernel(x1, x2, sigma):
    import numpy as np

    return np.exp(-np.linalg.norm(x1 - x2) ** 2 / (2 * sigma ** 2))


def cotangent_kernel(Nodes, Meshes):
    import numpy as np
    from scipy.sparse import lil_matrix

    num_nodes = Nodes.shape[0]
    W = lil_matrix((num_nodes, num_nodes))
    print('Constructing sparse affinity matrix...')
    for c in Meshes:
        # Obtain vertices which comprise face
        v0, v1, v2 = Nodes[c[0]], Nodes[c[1]], Nodes[c[2]]

        # Obtain cotangents of angles
        cot0 = np.dot(v1-v0, v2-v0) / np.linalg.norm(np.cross(v1-v0, v2-v0))
        cot1 = np.dot(v2-v1, v0-v1) / np.linalg.norm(np.cross(v2-v1, v0-v1))
        cot2 = np.dot(v0-v2, v1-v2) / np.linalg.norm(np.cross(v0-v2, v1-v2))

        # Update weight matrix accordingly
        W[c[1], c[2]] += cot0
        W[c[2], c[1]] += cot0
        W[c[0], c[2]] += cot1
        W[c[2], c[0]] += cot1
        W[c[0], c[1]] += cot2
        W[c[1], c[0]] += cot2

    return W


def inverse_distance(x1, x2, epsilon):
    import numpy as np

    return 1.0/(np.linalg.norm(x1 - x2) + epsilon)

