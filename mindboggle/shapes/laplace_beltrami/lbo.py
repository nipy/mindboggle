#!/usr/bin/python
"""
Computing the Laplace-Beltrami Spectrum of a given structure using linear FEM method

1. old_fem_laplacian (Forrest's old version. Has bugs. Kept temporarily.)
2. fem_laplacians (Forrest's new version that gives same output as Martin's.)

We follow the definitions and steps given in Martin Reuter et al.'s paper:
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
    - Dr. Martin Reuter, MIT, who provides his MATLAB code which is of great help 
      and explains his paper in very details               http://reuter.mit.edu/ 
    - Dr. Eric You Xu, Google, who explains to Forrest how eigenvalue problems are 
      solved numerically                                   http://www.youxu.info/ 

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

def computeAB(V, T):
    """Compute the matrices A and B for LBO

    Inputs
    ==========
    V: list of lists of 3 floats
        Vertices. Each element is the X, Y and Z coordinate of a vertex on the structure

    T: list of lists of 3 integers
        Triangles. Each element represents the 3 indices of vertices that form a triangle on the mesh 

    Outputs
    =============

    A : a numpy matrix
    B : a numpy matrix

    """
    import numpy as np
    from scipy import sparse

    [v, t] = map(np.array, [V, T])

    [vnum, tnum] = map(len, [V, T]) # numbers of vertices and triangles

    # linear local matrices on unit triangle:
    tB = (np.ones((3,3)) + np.eye(3) )/24.0;

    tA00 = np.array([[ 0.5,-0.5, 0.0], 
                     [-0.5, 0.5, 0.0],
                     [ 0.0, 0.0, 0.0]])

    tA11 = np.array([[ 0.5, 0.0,-0.5], 
                     [ 0.0, 0.0, 0.0],
                     [-0.5, 0.0, 0.5]])

    tA0110 = np.array([[ 1.0,-0.5,-0.5], 
                       [-0.5, 0.0, 0.5],
                       [-0.5, 0.5, 0.0]])

# replicate into third dimension for each triangle
    tB     = np.array([np.tile(tB    ,(1, 1)) for i in xrange(tnum)]) # 1st index is the 3rd index in MATLAB
    tA00   = np.array([np.tile(tA00  ,(1, 1)) for i in xrange(tnum)])
    tA11   = np.array([np.tile(tA11  ,(1, 1)) for i in xrange(tnum)])
    tA0110 = np.array([np.tile(tA0110,(1, 1)) for i in xrange(tnum)])

    # compute vertex coord and difference vector for each triangle
    v1 = v[t[:,0],:]
    v2 = v[t[:,1],:]
    v3 = v[t[:,2],:]
    v2mv1 = v2 - v1
    v3mv1 = v3 - v1

    def reshape_and_repmat(A):
        """For a given 1-D array A, does the two MATLAB code below 

        M = reshape(M,1,1,tnum);
        M = repmat(M,3,3);

        Please note that the a0 is a 3-D matrix. But the 3rd index in NumPy is the 1st index in MATLAB    

        Luckily, tnum is the size of A. 

        """
        return np.array([np.ones((3,3))*x for x in A])# may require an old enough version of numpy/python to support this semantics

    # compute length^2 of v3mv1 for each triangle
    a0 = np.sum(v3mv1 * v3mv1, axis=1)
    a0 = reshape_and_repmat(a0)

    # compute lenght^2 of v2mv1 for each triangle
    a1 = np.sum(v2mv1 * v2mv1, axis=1)
    a1 = reshape_and_repmat(a1)

    # compute dot product (v2mv1*v3mv1) for each triangle
    a0110 = np.sum(v2mv1 * v3mv1, axis=1)
    a0110 = reshape_and_repmat(a0110)

    # compute cross product and 2*vol for each triangle
    cr  = np.cross(v2mv1,v3mv1)
    vol = np.sqrt(np.sum(cr*cr, axis=1))
    # zero vol will cause division by zero below, so set to small value:
    vol_mean = np.mean(vol)
    vol = [vol_mean if x == 0 else x for x in vol]
    vol = reshape_and_repmat(vol)

    # construct all local A and B matrices (guess: for each triangle)
    localB = vol * tB
    localA = (1.0/vol) * (a0*tA00 + a1*tA11 - a0110 * tA0110)

    # construct row and col indices
    J = np.array([np.tile(x, (3,1))  for x in t]) # note, J in numpy is I in MATLAB after flattening, because numpy is row-major while MATLAB is column-major
    I = np.array([np.transpose(np.tile(x, (3,1)))  for x in t])

    # flatten arrays
    J_new = I.flatten() # note, swap I and J here.
    I_new = J.flatten()
    localA = localA.flatten()
    localB = localB.flatten()

    # construct sparse matrix
    A = sparse.csr_matrix((localA, (I_new, J_new)))
    B = sparse.csr_matrix((localB, (I_new, J_new)))

#    return [[I_new[i]+1, J_new[i]+1, localA[i], localB[i]] for i in xrange(54)]
    return A, B

def print_sparse_matrix(M):
    print "\nThe sparse matrix:\n"
    for Row in M.toarray().tolist():
        for E in Row:
            if E == 0:
                print "0\t".rjust(6),
            else:
                print '{0:2.4f}\t'.format(E),
        print ""

    print "The End.\n"

def test_computeAB():
    import numpy as np
    nodes = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    nodes = np.array(nodes)
    # Then, pick some faces.
    faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]] # note, all points must be on faces. O/w, you get singular matrix error when inverting D
    faces = np.array(faces)

    A, B = computeAB(nodes, faces)

    print_sparse_matrix(A)
    print_sparse_matrix(B)

def fem_laplacian(Nodes, Faces):
    """New Linear FEM laplacian code after studying Martin Reuter's code in MATLAB
    """

    num_nodes = len(Nodes)
    
    if num_nodes < 5: # too small
        print "The input size is too small. Skipped."
        return numpy.array([-1,-1,-1, -1, -1])
  
    A, B = computeAB(Nodes, Faces)    

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
    print "the linear FEM LBS is:", list(fem_laplacian(nodes, faces))
    
