#!/usr/bin/python
"""
Computing the Laplace-Beltrami Spectrum of a given structure using linear FEM method

1. old_fem_laplacian (Forrest's old version. Has bugs. Kept temporarily.)
2. linear_fem_laplacians (Forrest's new version. Requires computeAB.py in the same directory)

We follow the definitions and steps given in Reuter et al.'s
Discrete Laplace-Beltrami Operators for Shape Analysis and Segmentation (2009)

The information about using SciPy to solve generalized eigenvalue problem is
at: http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html

DEPENDENCE: Scipy 0.10 or later to solve generalized eigenvalue problem. 

I just realize that PEP 8 says Capitalized_Letters_with_Underscores look UGLY!
This is contrary to what I remembered before. 
So I am replacing variable names gradually.  

Known things that you need to pay attention (not to be considered as bugs):
1. In ``nodes``, do NOT provide coordinates of vertices that do NOT appear on 
the 3-D structure whose LBS is to be calculated. 
For example, do not use coordinates of all POINTS from a VTK file as ``nodes`` 
and some of the faces (e.g., sulcus fold) as ``faces``. 
This will cause singular matrix error when inverting matrixes because some rows 
are all zeros. 

Acknowledgments:
    - Dr. Martin Reuter, MIT, http://reuter.mit.edu/ (Who provides his MATLAB code which is of great help and explains his paper in very details.)
    - Dr. Eric You Xu, Google, http://www.youxu.info/ (Who explains to Forrest how eigenvalue problems are solved numerically.)

Authors:
    - Forrest Sheng Bao, 2012-2013  (forrest.bao@gmail.com)  http://fsbao.net
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License
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
    return V/3

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

def old_fem_laplacian(Nodes, Faces):
    """The portal function to compute geometric laplacian
    
    Parameters
    ----------
    Nodes : 2-D numpy array 
        Nodes[i] is the 3-D coordinates of nodes on a mesh 
    Faces : 2-D numpy array
        Faces[i] is a 3-element array containing the indexes of nodes 

    Returns
    -------
    eigenvalues : a list of floats
        The Laplacian-Beltrami Spectrum 
    
    Notes
    ------
    
    This is how Forrest got the steps from the paper:
    1. The FEM Laplacian problem is given as 
       A_{cot}\mathbf{f} = - \lambda B \mathbf{f} in the paper (the next equation after Eq. 6)
       We denote this equation in the docstring as Eq.A.
    2. Let L' = - B^{-1} A_{cot}
    3. Then Eq.A can be rewritten as 
        L' \mathbf{f} = \lambda \mathbf{f}
    4. Similarly to geometric Laplacian, the FEM Laplacian spectrum is then the 
       eigenvalues of L' .
 
    Steps:
    We could heavily reuse the code for geometric Laplacian.
    
    1. Compute W (can directly use Eliezer's cotangent kernel)
    2. Compute V = diag(v_1,...v_n) where v_i = \sum_{j\in N(i)} w_{ij} 
       and N(i) is the set of neighbors of node i.
    3. Compute stiffness matrix A = V - W. 
       Note that W and V are two cases for A_{cot}.
       A is -A_{cot}   (A_{cot} should be W - V)
    4. Compute the mass matrix B according to the paper.    
       B = P + Q where P[i,j] = (Area[x] + Area[y])/2 for x and y 
       are the two faces sharing the edge (i,j) (0's, otherwise), and
       Q[i,j] = (\sum_{k\in N(i)} Area[k] )/6 for i=j (0's, otherwise)
        
       I assume by the notation \sum_{k\in N(i)} |t_k| in the paper, 
       the authors mean total area of all triangles centered at node i.
       There is some ambiguity here because N(i) is the set of neighbor nodes 
       of node i (defined earlier in the paper) whereas t_k is a triangle.
       This is my best guess.  
        
    5. L = inv(B)*A
    """
    
    def gen_P(edges, faces_at_edges, Area, num_nodes):
        """Generate the P mentioned in pseudocode above
        """
       
        from scipy.sparse import lil_matrix
        P = lil_matrix((num_nodes, num_nodes))
#        P = numpy.zeros((num_nodes, num_nodes))
        for [i,j] in edges:
            P[i,j] = sum([Area[face] for face in faces_at_edges[(i,j)]]) # this line replaces the block commented below
            # =-------------------
#            facing_edges = faces_at_edges[(i,j)]
#            if len(facing_edges) == 1:  
#                [t1]= facing_edges
#                P[i,j] = Area[t1] 
#            else:
#                [t1,t2]= facing_edges
#                P[i,j] = Area[t1] + Area[t2]
            #=---------------------------------
    
        return P/12
        
    def gen_Q(edges, faces_at_edges, Area, num_nodes, Neighbor, Faces_at_Nodes):
        """Generate the Q mentioned in pseudocode above
        """
        from scipy.sparse import lil_matrix
        Q = lil_matrix((num_nodes, num_nodes))
        q = [sum([Area[k] for k in Faces_at_Nodes[i]]) for i in range(num_nodes)]
        Q.setdiag(q)        

        return Q/6
    
    import numpy
          
    num_nodes = len(Nodes)
    
    if num_nodes < 5: # too small
        print "The input size is too small. Skipped."
        return numpy.array([-1,-1,-1, -1, -1])
    
    import mindboggle.utils.kernels
    W = mindboggle.utils.kernels.cotangent_kernel(Nodes, Faces)
    W /= 2
    
    import mindboggle.utils.mesh
    Neighbor = mindboggle.utils.mesh.find_neighbors(Faces, num_nodes)
     
    V = gen_V(Faces, W, Neighbor)
    A = V - W # the stiffness matrix
    
    Area = area(Nodes, Faces)
    Faces_at_Nodes = mindboggle.utils.mesh.find_faces_at_vertices(Faces, num_nodes)
    # up to this point, the computation is the same as in geometric Laplacian
    
    faces_at_edges = mindboggle.utils.mesh.find_faces_at_edges(Faces)
    edges = mindboggle.utils.mesh.find_edges(Faces.tolist())
    
    P = gen_P(edges, faces_at_edges, Area, num_nodes)
    Q = gen_Q(edges, faces_at_edges, Area, num_nodes, Neighbor, Faces_at_Nodes)
    B = P + Q

    from scipy.sparse.linalg import eigsh, eigs 
    # note eigs is for nonsymmetric matrixes while eigsh is for  real-symmetric or complex-hermitian matrices
    
    eigenvalues, eigenvectors = eigsh(A, k=3, M=B, which="SM")
    
    return eigenvalues

def linear_fem_laplacian(Nodes, Faces):
    """New Linear FEM laplacian code after studying Martin Reuter's code in MATLAB
    """

    num_nodes = len(Nodes)
    
    if num_nodes < 5: # too small
        print "The input size is too small. Skipped."
        return numpy.array([-1,-1,-1, -1, -1])

    import computeAB
    
    A, B = computeAB.computeAB(Nodes, Faces)    

    from scipy.sparse.linalg import eigsh, eigs 
    # note eigs is for nonsymmetric matrixes while eigsh is for  real-symmetric or complex-hermitian matrices
    
    eigenvalues, eigenvectors = eigsh(A, k=3, M=B, which="SM")

    return eigenvalues

if __name__ == "__main__":
    import numpy as np
    # You should get different output if you only change the coordinates of nodes.
    # If you do NOT see the changes, you are computing graph laplacian.
    

    # Use some vertices on a cube. First, define a cube.  
    nodes = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    nodes = np.array(nodes)
    # Then, pick some faces. 
    faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]] # note, all points must be on faces. O/w, you get singular matrix error when inverting D
    faces = np.array(faces)
    
#    print "the FEM LBS is:", list(old_fem_laplacian(nodes, faces))
    print "the linear FEM LBS is:", list(linear_fem_laplacian(nodes, faces))
    
