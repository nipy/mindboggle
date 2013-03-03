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


if __name__ == "__main__":
    import numpy as np
    nodes = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    nodes = np.array(nodes)
    # Then, pick some faces.
    faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]] # note, all points must be on faces. O/w, you get singular matrix error when inverting D
    faces = np.array(faces)

    A, B = computeAB(nodes, faces)

    print_sparse_matrix(A)
    print_sparse_matrix(B)
