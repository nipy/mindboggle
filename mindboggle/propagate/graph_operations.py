import numpy as np
import networkx as nx
from scipy.sparse import lil_matrix

###############################################################################
# -----------------------------------------------------------------------------
#     Kernels
# -----------------------------------------------------------------------------
###############################################################################

def rbf_kernel(x1, x2, sigma):
    return np.exp(-np.linalg.norm(x1 - x2) ** 2 / (2 * sigma ** 2))

def cotangent_kernel(Nodes, Meshes):
    num_nodes = Nodes.shape[0]
    W = lil_matrix((num_nodes, num_nodes))
    print 'Constructing sparse affinity matrix...'
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
	return 1.0/(np.linalg.norm(x1 - x2) + epsilon)

###############################################################################
# -----------------------------------------------------------------------------
#     Diagonal degree matrix
# -----------------------------------------------------------------------------
###############################################################################
def compute_diagonal_degree_matrix(W, inverse=False, square_root=False):
    """
    Compute diagonal degree matrix.

    Input
    =====
    W: N x N sparse matrix in csr format
                    Affinity matrix
                inverse: boolean
                    Compute inverse of diagonal degree matrix?
                square_root: boolean
                    Compute square root of diagonal degree matrix?

    Return
    ======
    ddm: N x N sparse matrix in csr format
    Diagonal matrix

    """
    ddm = lil_matrix((W.shape[0], W.shape[0]))

    if inverse:
        if not square_root:
            ddm.setdiag(1 / W.sum(axis=1))
        else:
            ddm.setdiag(np.sqrt(1 / W.sum(axis=1)))

    else:
        ddm.setdiag(W.sum(axis=1))

    return ddm.tocsr()

###############################################################################
# -----------------------------------------------------------------------------
#     Graph Laplacian
# -----------------------------------------------------------------------------
###############################################################################
def graph_laplacian(W, which='norm1'):
    """
    Compute normalized and unnormalized graph Laplacians.

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

    """

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

###############################################################################
# -----------------------------------------------------------------------------
#     Matrix weights and affinity matrix
# -----------------------------------------------------------------------------
###############################################################################

def compute_weights(Nodes, Meshes, kernel=rbf_kernel, add_to_graph=True,
                    G=nx.Graph(), sigma=20):
    """
    Construct weighted edges of a graph and compute an affinity matrix.

    Input
    =====
    Nodes: numpy array
    Meshes: numpy array
    kernel: function which determines weights of edges
        rbf_kernel       - Gaussian kernel, with parameter sigma
        cotangent_kernel - weight calculation for Laplace_Beltrami_Operator
        inverse_distance - additional kernel where the weight is the inverse
                           of the distance between two nodes
    add_to_graph:  boolean if weights should be added to graph
    G:  networkx graph
    sigma:  float (parameter for rbf_kernel)

    Output
    ======
    G:  networkx graph
    affinity_matrix:  numpy array (sparse affinity matrix)

    """

    if kernel is rbf_kernel or kernel is inverse_distance:
        print('Computing weights using {} kernel with parameter = {}'.format(
              kernel, sigma))

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
            edge_mat = np.hstack((edges[0].T[:, np.newaxis], edges[1].T[:, np.newaxis]))
            weighted_edges = np.asarray([[edge_mat[i,0], edge_mat[i,1],affinity_matrix[edge_mat[i]]] for i in xrange(affinity_matrix.shape[0])])
            print('Adding weighted edges to the graph...')
            G.add_weighted_edges_from(weighted_edges)

    if add_to_graph:
        return G, affinity_matrix.tocsr()
    else:
        return affinity_matrix.tocsr()
