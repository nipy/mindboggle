#!/usr/bin/python
"""
Compute the Laplace-Beltrami Spectrum (LBS) using a linear finite element method.

We follow the definitions and steps given in Martin Reuter et al.'s 2009 paper:
    ``Discrete Laplace-Beltrami Operators for Shape Analysis and Segmentation''


Dependency:
    Scipy 0.10 or later to solve the generalized eigenvalue problem.
    Information about using Scipy to solve a generalized eigenvalue problem:
    http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html

NOTE ::
    For ``points``, only include coordinates of vertices in the 3-D structure
    whose LBS is to be calculated. For example, do not use coordinates of all
    POINTS from a VTK file as ``points`` and only corresponding faces as ``faces``.
    Otherwise this will cause a singular matrix error when inverting matrices
    because some rows are all zeros.

Acknowledgments:
    - Dr. Martin Reuter, MIT, helped us to better understand his articles about
      the Laplace-Beltrami operator and us with provided his MATLAB code.
    - Dr. Eric You Xu, Google (http://www.youxu.info/),
      explained how eigenvalue problems are solved numerically.

Authors:
    - Forrest Sheng Bao, 2012-2013  (forrest.bao@gmail.com)  http://fsbao.net
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)
    - Martin Reuter, 2009, http://reuter.mit.edu/ (original MATLAB code)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def computeAB(points, faces):
    """
    Compute matrices for the Laplace-Beltrami operator.

    The matrices correspond to A and B from Reuter's 2009 article.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex

    faces : list of lists of 3 integers
        each list contains indices to vertices that form a triangle on the mesh

    Returns
    -------
    A : numpy matrix
    B : numpy matrix

    Examples
    --------
    >>> # Define a cube, then compute A and B on a selection of faces.
    >>> import numpy as np
    >>> from mindboggle.shapes.laplace_beltrami import computeAB, print_sparse_matrix
    >>> points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    >>> points = np.array(points)
    >>> faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]] # note, all points must be on faces. O/w, you get singular matrix error when inverting D
    >>> faces = np.array(faces)
    >>> #
    >>> A, B = computeAB(points, faces)
    >>> #
    >>> print_sparse_matrix(A)
        The sparse matrix:

        1.5000	-1.0000	-0.5000	    0	    0	    0	    0	    0
        -1.0000	2.0000	    0	    0	-0.5000	    0	    0	-0.5000
        -0.5000	    0	2.0000	-0.5000	-1.0000	    0	    0	    0
            0	    0	-0.5000	2.5607	-0.3536	-1.2071	-0.5000	    0
            0	-0.5000	-1.0000	-0.3536	1.8536	    0	    0	    0
            0	    0	    0	-1.2071	    0	1.2071	    0	    0
            0	    0	    0	-0.5000	    0	    0	0.5000	    0
            0	-0.5000	    0	    0	    0	    0	    0	0.5000
    >>> print_sparse_matrix(B)
        The sparse matrix:

        0.2500	0.0833	0.0417	    0	0.0833	    0	    0	0.0417
        0.0833	0.1667	    0	    0	0.0417	    0	    0	0.0417
        0.0417	    0	0.1667	0.0417	0.0833	    0	    0	    0
            0	    0	0.0417	0.2845	0.1006	0.1006	0.0417	    0
        0.0833	0.0417	0.0833	0.1006	0.3679	0.0589	    0	    0
            0	    0	    0	0.1006	0.0589	0.2012	0.0417	    0
            0	    0	    0	0.0417	    0	0.0417	0.0833	    0
        0.0417	0.0417	    0	    0	    0	    0	    0	0.0833

    """
    import numpy as np
    from scipy import sparse

    points = np.array(points)
    faces = np.array(faces)
    nfaces = faces.shape[0]

    # Linear local matrices on unit triangle:
    tB = (np.ones((3,3)) + np.eye(3)) / 24.0

    tA00 = np.array([[ 0.5,-0.5, 0.0],
                     [-0.5, 0.5, 0.0],
                     [ 0.0, 0.0, 0.0]])

    tA11 = np.array([[ 0.5, 0.0,-0.5],
                     [ 0.0, 0.0, 0.0],
                     [-0.5, 0.0, 0.5]])

    tA0110 = np.array([[ 1.0,-0.5,-0.5],
                       [-0.5, 0.0, 0.5],
                       [-0.5, 0.5, 0.0]])

    # Replicate into third dimension for each triangle
    # (For tB, 1st index is the 3rd index in MATLAB.):
    tB = np.array([np.tile(tB, (1, 1)) for i in xrange(nfaces)])
    tA00 = np.array([np.tile(tA00, (1, 1)) for i in xrange(nfaces)])
    tA11 = np.array([np.tile(tA11, (1, 1)) for i in xrange(nfaces)])
    tA0110 = np.array([np.tile(tA0110,(1, 1)) for i in xrange(nfaces)])

    # Compute vertex coordinates and a difference vector for each triangle:
    v1 = points[faces[:, 0], :]
    v2 = points[faces[:, 1], :]
    v3 = points[faces[:, 2], :]
    v2mv1 = v2 - v1
    v3mv1 = v3 - v1

    def reshape_and_repeat(A):
        """
        For a given 1-D array A, run the MATLAB code below.

            M = reshape(M,1,1,nfaces);
            M = repmat(M,3,3);

        Please note that a0 is a 3-D matrix, but the 3rd index in NumPy
        is the 1st index in MATLAB.  Fortunately, nfaces is the size of A.

        """
        return np.array([np.ones((3,3))*x for x in A])

    # Compute length^2 of v3mv1 for each triangle:
    a0 = np.sum(v3mv1 * v3mv1, axis=1)
    a0 = reshape_and_repeat(a0)

    # Compute length^2 of v2mv1 for each triangle:
    a1 = np.sum(v2mv1 * v2mv1, axis=1)
    a1 = reshape_and_repeat(a1)

    # Compute dot product (v2mv1*v3mv1) for each triangle:
    a0110 = np.sum(v2mv1 * v3mv1, axis=1)
    a0110 = reshape_and_repeat(a0110)

    # Compute cross product and 2*vol for each triangle:
    cr  = np.cross(v2mv1,v3mv1)
    vol = np.sqrt(np.sum(cr*cr, axis=1))
    # zero vol will cause division by zero below, so set to small value:
    vol_mean = np.mean(vol)
    vol = [vol_mean if x == 0 else x for x in vol]
    vol = reshape_and_repeat(vol)

    # Construct all local A and B matrices (guess: for each triangle):
    localB = vol * tB
    localA = (1.0/vol) * (a0*tA00 + a1*tA11 - a0110*tA0110)

    # Construct row and col indices.
    # (Note: J in numpy is I in MATLAB after flattening,
    #  because numpy is row-major while MATLAB is column-major.)
    J = np.array([np.tile(x, (3,1)) for x in faces])
    I = np.array([np.transpose(np.tile(x, (3,1))) for x in faces])

    # Flatten arrays and swap I and J:
    J_new = I.flatten()
    I_new = J.flatten()
    localA = localA.flatten()
    localB = localB.flatten()

    # Construct sparse matrix:
    A = sparse.csr_matrix((localA, (I_new, J_new)))
    B = sparse.csr_matrix((localB, (I_new, J_new)))

    return A, B

def print_sparse_matrix(M):
    """
    Print sparse matrix.

    Note ::
        print() command not suitable for Python 3.0.
    """
    print("\nThe sparse matrix:\n")
    for Row in M.toarray().tolist():
        for E in Row:
            if E == 0:
                print "0\t".rjust(6),
            else:
                print "{0:2.4f}\t".format(E),
        print("")

def compute_area(points, faces):
    """
    Compute the areas of all triangles on the mesh.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh

    Returns
    -------
    area: 1-D numpy array
        area[i] is the area of the i-th triangle

    """
    import numpy as np

    area = np.zeros(len(faces))

    points = np.array(points)

    for i, triangle in enumerate(faces):

        a = np.linalg.norm(points[triangle[0]] - points[triangle[1]])
        b = np.linalg.norm(points[triangle[1]] - points[triangle[2]])
        c = np.linalg.norm(points[triangle[2]] - points[triangle[0]])
        s = (a+b+c) / 2.0

        area[i] = np.sqrt(s*(s-a)*(s-b)*(s-c))

    return area

def area_normalize(points, faces, spectrum):
    """Normalize the spectrum for one shape using areas as suggested in Reuter et al. 2006

    Parameters
    ------------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    spectrum : list of floats
        LB spectrum of a given shape defined by _points_ and _faces_

    Returns
    -----------
    new_spectrum : list of floats
        LB spectrum normalized by area
    """

    area = compute_area(points, faces)
    total_area = sum(area) # the area of the entire shape

    new_spectrum = [x/total_area for x in spectrum]

    return new_spectrum

def fem_laplacian(points, faces, n_eigenvalues=200, normalization=None):
    """
    Linear FEM laplacian code after Martin Reuter's MATLAB code.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    n_eigenvalues : integer
        number of eigenvalues to return
    normalization : string
        the method used to normalize eigenvalues (default: None)
        if "area", use area of the 2D structure as mentioned in Reuter et al. 2006

    Returns
    -------
    spectrum : list
        first n_eigenvalues eigenvalues for Laplace-Beltrami spectrum

    Examples
    --------
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> # Define a cube:
    >>> points = [[0,0,0], [0,1,0], [1,1,0], [1,0,0],
    >>>           [0,0,1], [0,1,1], [1,1,1], [1,0,1]]
    >>> faces = [[0,1,2], [2,3,0], [4,5,6], [6,7,4], [0,4,7], [7,3,0],
    >>>          [0,4,5], [5,1,0], [1,5,6], [6,2,1], [3,7,6], [6,2,3]]
    >>> print("The un-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=3)))
        The un-normalized linear FEM Laplace-Beltrami Spectrum is:
        [-2.74238008172841e-16, 4.583592135001265, 4.800000000000001]
    >>> print("The area-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=3, normalization="area")))
        The area-normalized linear FEM Laplace-Beltrami Spectrum is:
        [-7.4014869016002383e-16, 0.76393202250021075, 0.80000000000000049]
    >>> # Spectrum for a single fold:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.utils.mesh import remove_faces, renumber_faces
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> faces, lines, indices, points, npoints, scalars, name, input_vtk = read_vtk(fold_file)
    >>> indices = [i for i,x in enumerate(scalars) if x != -1]
    >>> points = np.array(points)
    >>> points = points[indices]
    >>> faces = remove_faces(faces, indices)
    >>> faces = renumber_faces(faces, indices)
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=3, normalization="area")))

    """
    from scipy.sparse.linalg import eigsh

    from mindboggle.shapes.laplace_beltrami import computeAB

    min_n_eigenvalues = 3
    npoints = len(points)

    if npoints < min_n_eigenvalues:
        print "The input size {0} is smaller than n_eigenvalue {1}. Skipped.".\
            format(npoints, n_eigenvalues)
        return None
    elif npoints < n_eigenvalues:
        n_eigenvalues = npoints

    A, B = computeAB(points, faces)

    # Note: eigs is for nonsymmetric matrices while
    #       eigsh is for real-symmetric or complex-Hermitian matrices.
    eigenvalues, eigenvectors = eigsh(A, k=n_eigenvalues, M=B, sigma=0)

    spectrum = eigenvalues.tolist()

    if normalization == "area":
        spectrum = area_normalize(points, faces, spectrum)

    return spectrum

def fem_laplacian_from_labels(vtk_file, n_eigenvalues=200, normalization=None):
    """
    Compute linear FEM Laplace-Beltrami spectra from each labeled region in a file.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file containing index scalars (labels)
    n_eigenvalues : integer
        number of eigenvalues to return
    normalization : string
        the method used to normalize eigenvalues (default: None)
        if "area", use area of the 2D structure as mentioned in Reuter et al. 2006

    Returns
    -------
    spectrum_lists : list of lists
        first n_eigenvalues eigenvalues for each label's Laplace-Beltrami spectrum

    Examples
    --------
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'tests', 'cube.vtk')
    >>> n_eigenvalues = 3
    >>> print("The un-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian_from_labels(vtk_file, n_eigenvalues)))
        The un-normalized linear FEM Laplace-Beltrami Spectrum is:
        Reduced 12 to 0 triangular faces
        The input size 2 is smaller than n_eigenvalue 3. Skipped.
        Reduced 12 to 6 triangular faces
        [None, [3.2049378106476917e-17, 4.583592135001263, 4.800000000000005]]
    >>> print("The area-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian_from_labels(vtk_file, n_eigenvalues, normalization="area")))
        The area-normalized linear FEM Laplace-Beltrami Spectrum is:
        Reduced 12 to 0 triangular faces
        The input size 2 is smaller than n_eigenvalue 3. Skipped.
        Reduced 12 to 6 triangular faces
        [None, [-2.1366252070792216e-17, 1.5278640450004206, 1.6000000000000001]]
    """
    import numpy as np

    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.utils.mesh import remove_faces, renumber_faces
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    # Read VTK surface mesh file:
    faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(vtk_file)
    points = np.array(points)

    # Loop through labeled regions:
    spectrum_lists = []
    for label in np.unique(labels):

        # Extract points and renumber faces for the labeled region:
        indices = [i for i,x in enumerate(labels) if x == label]
        if indices:
            label_points = points[indices]
            label_faces = remove_faces(faces, indices)
            if label_faces:
                label_faces = renumber_faces(label_faces, indices)

                # Compute Laplace-Beltrami spectrum for the labeled region:
                spectrum = fem_laplacian(label_points, label_faces, n_eigenvalues,
                                         normalization)

                # Append to a list of lists of spectra:
                spectrum_lists.append(spectrum)
            else:
                spectrum_lists.append([])
        else:
            spectrum_lists.append([])

    return spectrum_lists


if __name__ == "__main__":

    import numpy as np
    # You should get different outputs if you change the coordinates of points.
    # If you do NOT see changes, you may be computing the graph Laplacian.

    # Define a cube:
    points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1],
              [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    # Pick some faces:
    faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]

    print("The un-normalized linear FEM Laplace-Beltrami Spectrum is:\n\t{0}\n".format(
        fem_laplacian(points, faces, n_eigenvalues=3)))

    print("The area-normalized linear FEM Laplace-Beltrami Spectrum is:\n\t{0}\n".format(
        fem_laplacian(points, faces, n_eigenvalues=3, normalization="area")))


    """
    def compute_V(W, neighbors):
        ""
        Compute V as in Martin Reuter's 2009 paper.

        Parameters
        ----------
        W: 2-D numpy array
            W[i,j] is w_{ij} in Eq. (3) of Reuter's 2009 paper
        neighbors: 2-D list
            neighbors[i] gives the list of neighbors of i on the mesh,
            in indices to vertices

        Returns
        -------
        V : sparse diagonal matrix
            described in Reuter's 2009 paper

        ""
        from scipy.sparse import lil_matrix

        npoints = W.shape[0]
        V = lil_matrix((npoints, npoints))

        # v: 1-D list
        # v_i = \sum_{j\in N(i)} w_{ij},
        #       where N(i) is the set of neighbors of vertex i
        v = [sum([W[i,j] for j in neighbors[i]]) for i in range(npoints)]

        V.setdiag(v)

        return V / 3

    def old_fem_laplacian(points, faces):
        ""The portal function to compute geometric laplacian

        Parameters
        ----------
        points : 2-D numpy array
            points[i] is the 3-D coordinates of points on a mesh
        faces : 2-D numpy array
            faces[i] is a 3-element array containing the indices of points

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
           B = P + Q where P[i,j] = (area[x] + area[y])/2 for x and y
           are the two faces sharing the edge (i,j) (0's, otherwise), and
           Q[i,j] = (\sum_{k\in N(i)} area[k] )/6 for i=j (0's, otherwise)

           I assume by the notation \sum_{k\in N(i)} |t_k| in the paper,
           the authors mean total area of all triangles centered at node i.
           There is some ambiguity here because N(i) is the set of neighbor points
           of node i (defined earlier in the paper) whereas t_k is a triangle.
           This is my best guess.

        5. L = inv(B)*A
        ""

        def gen_P(edges, faces_at_edges, area, npoints):
            ""Generate the P mentioned in pseudocode above
            ""

            from scipy.sparse import lil_matrix
            P = lil_matrix((npoints, npoints))
            #        P = numpy.zeros((npoints, npoints))
            for [i,j] in edges:
                P[i,j] = sum([area[face] for face in faces_at_edges[(i,j)]]) # this line replaces the block commented below
                # =-------------------
            #            facing_edges = faces_at_edges[(i,j)]
            #            if len(facing_edges) == 1:
            #                [t1]= facing_edges
            #                P[i,j] = area[t1]
            #            else:
            #                [t1,t2]= facing_edges
            #                P[i,j] = area[t1] + area[t2]
            #=---------------------------------

            return P/12

        def gen_Q(edges, faces_at_edges, area, npoints, neighbors, faces_at_points):
            ""Generate the Q mentioned in pseudocode above
            ""
            from scipy.sparse import lil_matrix
            Q = lil_matrix((npoints, npoints))
            q = [sum([area[k] for k in faces_at_points[i]]) for i in range(npoints)]
            Q.setdiag(q)

            return Q/6

        import numpy

        npoints = len(points)

        if npoints < 5: # too small
            print("The input size is too small. Skipped.")
            return numpy.array([-1,-1,-1, -1, -1])

        import mindboggle.utils.kernels
        W = mindboggle.utils.kernels.cotangent_kernel(points, faces)
        W /= 2

        import mindboggle.utils.mesh
        neighbors = mindboggle.utils.mesh.find_neighbors(faces, npoints)

        V = compute_V(faces, W, neighbors)
        A = V - W # the stiffness matrix

        area = compute_area(points, faces)
        faces_at_points = mindboggle.utils.mesh.find_faces_at_vertices(faces, npoints)
        # up to this point, the computation is the same as in geometric Laplacian

        faces_at_edges = mindboggle.utils.mesh.find_faces_at_edges(faces)
        edges = mindboggle.utils.mesh.find_edges(faces.tolist())

        P = gen_P(edges, faces_at_edges, area, npoints)
        Q = gen_Q(edges, faces_at_edges, area, npoints, neighbors, faces_at_points)
        B = P + Q

        from scipy.sparse.linalg import eigsh, eigs
        # note eigs is for nonsymmetric matrices while eigsh is for  real-symmetric or complex-hermitian matrices

        eigenvalues, eigenvectors = eigsh(A, k=3, M=B, which="SM")

        return eigenvalues

    """
