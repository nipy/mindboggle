#!/usr/bin/python
"""
Computing the Laplace-Beltrami Spectrum of a given structure. 

1. Geometric Laplacians 
2. FEM Laplacians (to be implemented)

Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net
    - Eliezer Stavsky  (eli.stavsky@gmail.com)

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

We follow the definitions and steps given in Reuter et al.'s 
Discrete Laplace-Beltrami Operators for Shape Analysis and Segmentation

"""

def gen_V(Meshes, W, Neighbor):
    """Computer V as in Reuter et al.'s paper
    
    Parameters
    -------------
    Meshes : 2-D numpy array
        Meshes[i] is a 3-element array containing the indexes of nodes
    W: 2-D numpy array
        W[i,j] is w_{ij} in Eq. (3) of Reuter's paper
    v: a 1-D list
        v_i = \sum_{j\in N(i)} w_{ij} and N(i) is the set of neighbors of node i.
    Neighbor: a 2-D list
        Neighbor[i] gives the list of neighbors of node i on the mesh, in node indexes. 
        
    Returns 
    ----------
    V : a sparse diagnonal matrix
        As described in Reuter's paper
    
    """
    from scipy.sparse import lil_matrix
    num_nodes = W.shape[0]
    V=lil_matrix((num_nodes, num_nodes))
    v = [sum([W[i,j] for j in Neighbor[i]]) for i in range(num_nodes)]
    V.setdiag(v)
    return V

def area(Nodes, Meshes):
    """Compute the areas all triangles on the mesh

    Parameters
    -------------
    
    Nodes : 2-D numpy array 
        Nodes[i] is the 3-D coordinates of nodes on a mesh 
    Meshes : 2-D numpy array
        Meshes[i] is a 3-element array containing the indexes of nodes 
    
    Returns
    --------
    Area: A 1-D numpy array
        Area[i] is the area of the i-th triangle 
    
    Notes
    ------
    Using Eliezer's compute_face_measures() in lbopy.py to do so.  
     
    """
    import numpy as np
    Area = np.zeros(Meshes.shape[0])
    i = 0
    for Triangle in Meshes: # Shoot, I cannot use enumerate() for numpy array
        a = np.linalg.norm(Nodes[Triangle[0]] - Nodes[Triangle[1]])
        b = np.linalg.norm(Nodes[Triangle[1]] - Nodes[Triangle[2]])
        c = np.linalg.norm(Nodes[Triangle[2]] - Nodes[Triangle[0]])
        s = (a+b+c)/2.0

        Area[i] = np.sqrt(s*(s-a)*(s-b)*(s-c))
        i += 1
    return Area

def masses(Nodes, Areas, Faces_at_Nodes):
    """Computer the mass matrix D = diag(d_1, ..., d_n) where 
       d_i = a(i)/3 and a(i) is the area of all triangles at node i 
       (Here we adopt Eq. (4) of Reuter's paper)    

    Parameters
    -----------
    
    Nodes : 2-D numpy array 
        Nodes[i] is the 3-D coordinates of nodes on a mesh 
    Meshes : 2-D numpy array
        Meshes[i] is a 3-element array containing the indexes of nodes 
    Area: A 1-D numpy array
        Area[i] is the area of the i-th triangle 
    Faces_at_Nodes: a 2-D list
        Faces_at_Nodes[i] is a list of IDs of faces at node i.
    d: a list of floats
        The sequence d_i, ..., d_n in Eq.(4)
    
    Returns
    ---------
    D: a sparse diagonal matrix
        The mass matrix
    
    """
#    import numpy as np
    from scipy.sparse import lil_matrix
    num_nodes = Nodes.shape[0]
    
    d = [sum([Areas[j] for j in Faces_at_Nodes[i]]) for i in range(num_nodes)]
    D = lil_matrix((num_nodes, num_nodes))
    D.setdiag(d)
    D /= 3 
    return D

def geometric_laplacian(Nodes, Faces):
    """The portal function to compute geometric laplacian
    
    Parameters
    ----------
    Nodes : 2-D numpy array 
        Nodes[i] is the 3-D coordinates of nodes on a mesh 
    Faces : 2-D numpy array
        Faces[i] is a 3-element array containing the indexes of nodes 

    Returns
    -------
    Spectrum : a list of floats
        The Laplacian-Beltrami Spectrum 
    
    Notes
    ------
    
    This algorithm is described in Section 2.1.1 Discrete geometric Laplacians
    Steps:
    1. Compute W (can directly use Eliezer's cotangent kernel)
    2. Compute V = diag(v_1,...v_n) where v_i = \sum_{j\in N(i)} w_{ij} 
       and N(i) is the set of neighbors of node i.
    3. Compute stiffness matrix A = V - W 
    4. Compute the mass matrix D = diag(d_1, ..., d_n) where 
       d_i = a(i)/3 and a(i) is the area of all triangles at node i 
       (Here we adopt Eq. (4) of Reuter's paper)
    5. L = inv(D)*A
    
    """
    import mindboggle.utils.kernels
    W = mindboggle.utils.kernels.cotangent_kernel(Nodes, Faces)
    W /= 2
    
    import mindboggle.utils.mesh_operations
    Neighbor = mindboggle.utils.mesh_operations.find_neighbors(Faces, len(Nodes))
     
    V = gen_V(Faces, W, Neighbor)
    A = V - W # the stiffness matrix
    Area = area(Nodes, Faces)
    
    Faces_at_Nodes = mindboggle.utils.mesh_operations.find_faces_at_vertexes(Faces, len(Nodes))
    D = masses(Nodes, Area, Faces_at_Nodes)
    
    import numpy
    L = numpy.linalg.inv(D)*A
    Spectrum = numpy.linalg.eig(L)
    
    return Spectrum