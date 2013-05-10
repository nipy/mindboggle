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

def wesd(EVAL1, EVAL2, Vol1, Vol2, show_error=False, N=3):
    """
    Weighted Spectral Distance. See Konukoglu et al. (2012)

    Parameters
    ---------------

    EVAL1 : numpy array of floats
        LB spectrum from the 1st shape

    EVAL2 : numpy array of floats
        LB spectrum from the 2nd shape

    Vol1: float
        Volume (area for 2D) of the 1st shape

    Vol2: float
        Volume (area for 2D) of the 2nd shape

    show_error: boolean
        Whether display the error of WESD using Eqs.(9) and (10) in Konukoglu et al. (2012).
        default: false

    N : integer
        The length of spetrum used (N>=3, default: 3)

    """
    # At present, algorithm doesn't return normalized result. It therefore doesn't require calculation of volume.

    d = 2.0 # " a surface is a 2d manifold. It doesn't matter that it is usually embedded in 3d Euclidean space. -Martin"
    Ball = 4.0/3*np.pi # For Three Dimensions
    p = 2.0

    Vol = np.amax((Vol1, Vol2))
    mu = np.amax(EVAL1[1], EVAL2[1])

    C = ((d+2)/(d*4*np.pi**2)*(Ball*Vol)**(2/d) - 1/mu)**p + ((d+2)/(d*4*np.pi**2)*(Ball*Vol/2)**(2/d) - 1/mu*(d/(d+4)))**p

    K = ((d+2)/(d*4*np.pi**2)*(Ball*Vol)**(2/d) - (1/mu)*(d/(d+2.64)))**p

    W = (C + K*(zeta(2*p/d,1) - 1 - .5**(2*p/d)))**(1/p) # the right-hand side of Eq.(8) or the equation right after Eq.(4)

    holder = 0
    for i in xrange(1, np.amin((len(EVAL1), len(EVAL2) )) ):
        holder += (np.abs(EVAL1[i] - EVAL2[i])/(EVAL1[i]*EVAL2[i]))**p
    WESD = holder ** (1/p)

    nWESD = WESD/W

    if show_error:
        WN = (C + K * (sum([n**(-1*2*p/d) for n in range(3,N+1)])))**(1/p)
        # the second term on the right-hand side of Eq.(9)
        print "Truncation error of WESD is: ", W - WN
        print "Truncation error of nWESD is: ", 1 -  WN/W

    return WESD

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
    >>> from mindboggle.utils.io_vtk import read_vtk, read_faces_points
    >>> from mindboggle.utils.mesh import reindex_faces_points
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> faces, points, npoints = read_faces_points(fold_file)
    >>> faces, points = reindex_faces_points(faces, points)
    >>> # Test LBO:
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=6, normalization="area")))
        [8.6598495496215578e-20, 4.214922171245502e-19, 1.7613177561957697e-05,
         3.9602772997696686e-05, 7.1740562223650042e-05, 8.5687655524969452e-05]
    """
    from scipy.sparse.linalg import eigsh

    from mindboggle.shapes.laplace_beltrami import computeAB

    min_n_eigenvalues = 10 * n_eigenvalues
    npoints = len(points)

    if npoints < min_n_eigenvalues:
        print "The input size {0} should be much larger than n_eigenvalue {1}. Skipped.".\
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

def fem_laplacian_from_labels(vtk_file, n_eigenvalues=3, normalization=None):
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
        first eigenvalues for each label's Laplace-Beltrami spectrum
    label_list : list of integers
        list of unique labels for which spectra are obtained

    Examples
    --------
    >>> # Spectrum for a single fold (one label):
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> print("The un-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian_from_labels(fold_file, n_eigenvalues=6)))
        ([[5.0385680900266565e-17, 2.452371976244268e-16,
           0.010247891019253531, 0.02304211720531372,
           0.04174087615603567, 0.04985572605660033]], [1.0])
    >>> print("The area-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian_from_labels(fold_file, n_eigenvalues=6,
    >>>                    normalization="area")))
        ([[8.6598495496215614e-20, 4.214922171245503e-19,
           1.7613177561957693e-05, 3.9602772997696571e-05,
           7.1740562223650029e-05, 8.5687655524969235e-05]], [1.0])
    >>> #
    >>> # Case of too few points:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'tests', 'cube.vtk')
    >>> n_eigenvalues = 6
    >>> print("{0}".format(fem_laplacian_from_labels(vtk_file, n_eigenvalues)))
        The input size 6 should be much larger than n_eigenvalue 6. Skipped.
        ([None], [0.0])
    >>> #
    >>> # Spectra for multiple sulci:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> print("The area-normalized linear FEM Laplace-Beltrami Spectra:\n")
    >>> print("{0}".format(fem_laplacian_from_labels(sulci_file, n_eigenvalues=6,
    >>>                    normalization="area")))

    """
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.utils.mesh import reindex_faces_points, remove_faces
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    min_n_eigenvalues = 10 * n_eigenvalues

    # Read VTK surface mesh file:
    faces, foo1, foo2, points, foo3, labels, foo4, foo5 = read_vtk(vtk_file)

    # Loop through labeled regions:
    ulabels = []
    [ulabels.append(int(x)) for x in labels if x not in ulabels if x != -1]
    label_list = []
    spectrum_lists = []
    for label in ulabels:

        # Extract points and renumber faces for the labeled region:
        indices = [i for i,x in enumerate(labels) if x == label]
        if len(indices) >= min_n_eigenvalues:
            label_faces = remove_faces(faces, indices)
            if label_faces:
                label_faces, label_points = reindex_faces_points(label_faces,
                                                                 points)

                # Compute Laplace-Beltrami spectrum for the labeled region:
                spectrum = fem_laplacian(label_points, label_faces,
                                         n_eigenvalues, normalization)

                # Append to a list of lists of spectra:
                spectrum_lists.append(spectrum)
                label_list.append(label)
        else:
            print "The input size {0} for label {1} is too small. Skipped.".\
                format(len(indices), label)

    return spectrum_lists, label_list


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
