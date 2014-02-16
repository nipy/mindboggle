#!/usr/bin/python
"""
Compute the Laplace-Beltrami spectrum using a linear finite element method.

We follow the definitions and steps given in Martin Reuter et al.'s 2009 paper:
    ``Discrete Laplace-Beltrami Operators for Shape Analysis and Segmentation''

Dependency:
    Scipy 0.10 or later to solve the generalized eigenvalue problem.
    Information about using Scipy to solve a generalized eigenvalue problem:
    http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html

NOTE ::

    For ``points``, only include coordinates of vertices in the 3-D structure
    whose spectrum is to be calculated. For example, do not use coordinates
    of all POINTS from a VTK file as ``points`` and use only corresponding
    faces as ``faces``. Otherwise this will cause a singular matrix error
    when inverting matrices because some rows are all zeros.

Acknowledgments:
    - Dr. Martin Reuter, MIT, provided his MATLAB code and assisted us in
      understanding his articles about the Laplace-Beltrami operator.
    - Dr. Eric You Xu, Google (http://www.youxu.info/),
      explained how eigenvalue problems are solved numerically.

Authors:
    - Martin Reuter, 2009, http://reuter.mit.edu/ (original MATLAB code)
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)
    - Forrest Sheng Bao, 2012-2013  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def computeAB(points, faces):
    """
    Compute matrices for the Laplace-Beltrami operator.

    The matrices correspond to A and B from Reuter's 2009 article.

    Note ::
        All points must be on faces. Otherwise, a singular matrix error
        is generated when inverting D.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex

    faces : list of lists of 3 integers
        each list contains indices to vertices that form a triangle on a mesh

    Returns
    -------
    A : numpy matrix
    B : numpy matrix

    Examples
    --------
    >>> # Define a cube, then compute A and B on a selection of faces.
    >>> import numpy as np
    >>> from mindboggle.shapes.laplace_beltrami import computeAB
    >>> points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    >>> points = np.array(points)
    >>> faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    >>> faces = np.array(faces)
    >>> #
    >>> A, B = computeAB(points, faces)
    >>> #
    >>> print(A.toarray())
        [[ 1.5  -1.   -0.5   0.    0.    0.    0.    0.  ]
         [-1.    2.    0.    0.   -0.5   0.    0.    -0.5  ]
         [-0.5   0.    2.   -0.5  -1.    0.    0.    0.  ]
         [ 0.    0.   -0.5   2.56066017 -0.35355339 -1.20710678 -0.5   0.  ]
         [ 0.   -0.5  -1.   -0.35355339  1.85355339  0.    0.    0.  ]
         [ 0.    0.    0.   -1.20710678  0.    1.20710678    0.    0.  ]
         [ 0.    0.    0.   -0.5   0.    0.    0.5    0.  ]
         [ 0.   -0.5   0.    0.    0.    0.    0.    0.5]]
    >>> print(B.toarray())
        [[ 0.25  0.08333333  0.04166667  0.    0.08333333  0.    0.    0.04166667]
         [ 0.08333333  0.16666667  0.    0.    0.04166667  0.    0.    0.04166667]
         [ 0.04166667  0.    0.16666667  0.04166667  0.08333333  0.    0.    0.  ]
         [ 0.    0.    0.04166667  0.2845178   0.10059223  0.10059223    0.04166667  0.  ]
         [ 0.08333333  0.04166667  0.08333333  0.10059223  0.36785113  0.05892557  0.  0.  ]
         [ 0.    0.    0.    0.10059223  0.05892557  0.20118446    0.04166667  0.  ]
         [ 0.    0.    0.    0.04166667  0.    0.04166667    0.08333333  0.  ]
         [ 0.04166667  0.04166667  0.    0.    0.    0.    0.     0.08333333]]

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
    # (for tB, 1st index is the 3rd index in MATLAB):
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


def area_normalize(points, faces, spectrum):
    """
    Normalize a spectrum using areas as suggested in Reuter et al. (2006)

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
    from mindboggle.utils.mesh import area_of_faces

    area = area_of_faces(points, faces)
    total_area = sum(area)

    new_spectrum = [x/total_area for x in spectrum]

    return new_spectrum


def wesd(EVAL1, EVAL2, Vol1, Vol2, show_error=False, N=3):
    """
    Weighted Spectral Distance. See Konukoglu et al. (2012)

    Note ::
        At present, algorithm doesn't return normalized result.
        It therefore doesn't require calculation of volume.

    Parameters
    ----------
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
        print("Truncation error of WESD is: {0}".format(W - WN))
        print("Truncation error of nWESD is: {1}".format(1 -  WN/W))

    return WESD


def fem_laplacian(points, faces, spectrum_size=10, normalization=None):
    """
    Compute linear finite-element method Laplace-Beltrami spectrum
    after Martin Reuter's MATLAB code.

    Note ::

        Compare fem_laplacian() with Martin Reuter's Matlab eigenvalues:

        fem_laplacian() results for Twins-2-1 left hemisphere (6 values):
        [4.829758648026221e-18,
         0.0001284173002467199,
         0.0002715181572272745,
         0.0003205150847159417,
         0.0004701628070486448,
         0.0005768904023010318]

        Martin Reuter's shapeDNA-tria Matlab code:
        {-4.7207711983791511358e-18 ;
         0.00012841730024672144738 ;
         0.00027151815722727096853 ;
         0.00032051508471592313632 ;
         0.0004701628070486902353  ;
         0.00057689040230097490998 }

        fem_laplacian() results for Twins-2-1 left postcentral (1022):
        [6.3469513010430304e-18,
         0.0005178862383467463,
         0.0017434911095630772,
         0.003667561767487686,
         0.005429017880363784,
         0.006309346984678924]

        Martin Reuter's Matlab code:
        {-2.1954862991027e-18 ;
         0.0005178862383468 ;
         0.0017434911095628 ;
         0.0036675617674875 ;
         0.0054290178803611 ;
         0.006309346984678 }

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    spectrum_size : integer
        number of eigenvalues to be computed (the length of the spectrum)
    normalization : string
        the method used to normalize eigenvalues ('area' or None)
        if "area", use area of the 2D structure as in Reuter et al. 2006

    Returns
    -------
    spectrum : list
        first spectrum_size eigenvalues for Laplace-Beltrami spectrum

    Examples
    --------
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> # Define a cube:
    >>> points = [[0,0,0], [0,1,0], [1,1,0], [1,0,0],
    >>>           [0,0,1], [0,1,1], [1,1,1], [1,0,1]]
    >>> faces = [[0,1,2], [2,3,0], [4,5,6], [6,7,4], [0,4,7], [7,3,0],
    >>>          [0,4,5], [5,1,0], [1,5,6], [6,2,1], [3,7,6], [6,2,3]]
    >>> fem_laplacian(points, faces, spectrum_size=3, normalization=None)
    [7.401486830834377e-17, 4.58359213500127, 4.799999999999998]
    >>> fem_laplacian(points, faces, spectrum_size=3, normalization="area")
    [1.2335811384723967e-17, 0.76393202250021175, 0.79999999999999949]
    >>> # Spectrum for entire left hemisphere of Twins-2-1:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_faces_points
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels',
    >>>                         'lh.labels.DKT25.manual.vtk')
    >>> faces, points, npoints = read_faces_points(vtk_file)
    >>> fem_laplacian(points, faces, spectrum_size=6, normalization=None)
    [4.829758648026222e-18,
     0.0001284173002467197,
     0.000271518157227274,
     0.00032051508471594065,
     0.0004701628070486444,
     0.0005768904023010318]
    >>> # Spectrum for Twins-2-1 left postcentral pial surface (22):
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.utils.mesh import remove_faces, reindex_faces_points
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> faces, u1,u2, points, u3, labels, u4,u5 = read_vtk(label_file)
    >>> I22 = [i for i,x in enumerate(labels) if x==22] # postcentral
    >>> faces = remove_faces(faces, I22)
    >>> faces, points, o1 = reindex_faces_points(faces, points)
    >>> #from mindboggle.utils.io_vtk import read_faces_points
    >>> #label_file = os.path.join(path, 'arno', 'labels', 'label22.vtk')
    >>> #faces, points, npoints = read_faces_points(label_file)
    >>> fem_laplacian(points, faces, spectrum_size=6, normalization=None)
    [6.3469513010430304e-18,
     0.0005178862383467462,
     0.0017434911095630795,
     0.0036675617674876856,
     0.005429017880363785,
     0.006309346984678933]
    >>> # Area-normalized spectrum for a single label (postcentral):
    >>> fem_laplacian(points, faces, spectrum_size=6, normalization="area")
    [1.1410192787181146e-21,
     9.3102680973672063e-08,
     3.1343504525679647e-07,
     6.5933366810380741e-07,
     9.7599836081654446e-07,
     1.1342589857996233e-06]

    """
    from scipy.sparse.linalg import eigsh, lobpcg
    import numpy as np

    from mindboggle.shapes.laplace_beltrami import computeAB

    #-----------------------------------------------------------------
    # Compute A and B matrices (from Reuter et al., 2009):
    #-----------------------------------------------------------------
    A, B = computeAB(points, faces)
    if A.shape[0] <= spectrum_size:
        print("The 3D shape has too few vertices ({0} <= {1}). Skip.".
              format(A.shape[0], spectrum_size))
        return None

    #-----------------------------------------------------------------
    # Use the eigsh eigensolver:
    #-----------------------------------------------------------------
    try :

        # eigs is for nonsymmetric matrices while
        # eigsh is for real-symmetric or complex-Hermitian matrices:
        eigenvalues, eigenvectors = eigsh(A, k=spectrum_size, M=B,
                                          sigma=0)
        spectrum = eigenvalues.tolist()

    #-----------------------------------------------------------------
    # Use the lobpcg eigensolver:
    #-----------------------------------------------------------------
    except RuntimeError:     
           
        print("eigsh() failed. Now try lobpcg.")
        print("Warning: lobpcg can produce different results from "
              "Reuter (2006) shapeDNA-tria software.")
        # Initial eigenvector values:
        init_eigenvecs = np.random.random((A.shape[0], spectrum_size))

        # maxiter = 40 forces lobpcg to use 20 iterations.
        # Strangely, largest=false finds largest eigenvalues
        # and largest=True gives the smallest eigenvalues:
        eigenvalues, eigenvectors =  lobpcg(A, init_eigenvecs, B=B,
                                            largest=True, maxiter=40)
        # Extract the real parts:
        spectrum = [value.real for value in eigenvalues]

        # For some reason, the eigenvalues from lobpcg are not sorted:
        spectrum.sort()

    #-----------------------------------------------------------------
    # Normalize by area:
    #-----------------------------------------------------------------
    if normalization == "area":
        spectrum = area_normalize(points, faces, spectrum)
        print("Compute area-normalized linear FEM Laplace-Beltrami spectrum")
    else:
        print("Compute linear FEM Laplace-Beltrami spectrum")

    return spectrum


def spectrum_of_largest(points, faces, spectrum_size=10, exclude_labels=[-1],
                        normalization=None, areas=None):
    """
    Compute Laplace-Beltrami spectrum on largest connected segment.

    In case a surface patch is fragmented, we select the largest fragment,
    remove extraneous triangular faces, and reindex indices.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    spectrum_size : integer
        number of eigenvalues to be computed (the length of the spectrum)
    exclude_labels : list of integers
        background values to exclude
    normalization : string
        the method used to normalize eigenvalues ('area' or None)
        if "area", use area of the 2D structure as in Reuter et al. 2006
    areas : numpy array or list of floats (or None)
        surface area scalar values for all vertices

    Returns
    -------
    spectrum : list
        first spectrum_size eigenvalues for Laplace-Beltrami spectrum

    Examples
    --------
    >>> # Spectrum for left postcentral + pars triangularis pial surfaces:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> from mindboggle.utils.mesh import remove_faces, reindex_faces_points
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_of_largest
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> spectrum_size = 6
    >>> exclude_labels = [-1]
    >>> normalization = None
    >>> faces, lines, indices, points, u1, labels, u2,u3 = read_vtk(label_file,
    >>>      return_first=True, return_array=True)
    >>> I20 = [i for i,x in enumerate(labels) if x==20] # pars triangularis
    >>> I22 = [i for i,x in enumerate(labels) if x==22] # postcentral
    >>> I22.extend(I20)
    >>> faces = remove_faces(faces, I22)
    >>> faces, points, o1 = reindex_faces_points(faces, points)
    >>> areas, u1 = read_scalars(area_file, True, True)
    >>> #
    >>> spectrum_of_largest(points, faces, spectrum_size, exclude_labels,
    >>>                     normalization, areas)
    [6.3469513010430304e-18,
     0.0005178862383467463,
     0.0017434911095630772,
     0.003667561767487686,
     0.005429017880363784,
     0.006309346984678924]
    >>> # View both segments:
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> scalars = np.zeros(np.shape(labels))
    >>> scalars[I22] = 1
    >>> vtk_file = 'test_two_labels.vtk'
    >>> write_vtk(vtk_file, points, indices, lines, faces,
    >>>           scalars, scalar_names='scalars', scalar_type='int')
    >>> plot_surfaces(vtk_file)

    """
    from scipy.sparse.linalg import eigsh, lobpcg
    import numpy as np

    from mindboggle.utils.segment import select_largest
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    if isinstance(areas, list):
        areas = np.array(areas)

    # Check to see if there are enough points:
    min_points_faces = spectrum_size
    npoints = len(points) 
    if npoints < min_points_faces or len(faces) < min_points_faces:
        print("The input size {0} ({1} faces) should be much larger "
              "than spectrum_size ({2})".
              format(npoints, len(faces), spectrum_size))
        return None
    else:

        #---------------------------------------------------------------------
        # Select the largest segment (connected set of indices):
        #---------------------------------------------------------------------
        points, faces = select_largest(points, faces, exclude_labels, areas,
                                       reindex=True)

        # Alert if the number of indices is small:
        if len(points) < min_points_faces:
            print("The input size {0} is too small.".format(len(points)))
            return None
        elif faces:

            #-----------------------------------------------------------------
            # Compute spectrum:
            #-----------------------------------------------------------------
            spectrum = fem_laplacian(points, faces, spectrum_size,
                                     normalization)
            return spectrum
        else:
            return None


def spectrum_from_file(vtk_file, spectrum_size=10, exclude_labels=[-1],
                       normalization=None, area_file=''):
    """
    Compute Laplace-Beltrami spectrum of a 3D shape in a VTK file.

    Parameters
    ----------
    vtk_file : string
        the input vtk file
    spectrum_size : integer
        number of eigenvalues to be computed (the length of the spectrum)
    exclude_labels : list of integers
        labels to be excluded
    normalization : string
        the method used to normalize eigenvalues ('area' or None)
        if "area", use area of the 2D structure as in Reuter et al. 2006
    area_file :  string
        name of VTK file with surface area scalar values

    Returns
    -------
    spectrum : list of floats
        first spectrum_size of Laplace-Beltrami spectrum

    Examples
    --------
    >>> # Spectrum for entire left hemisphere of Twins-2-1:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_from_file
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> spectrum_from_file(vtk_file, spectrum_size=6)
    [4.829758648026223e-18,
     0.00012841730024671977,
     0.0002715181572272744,
     0.00032051508471594173,
     0.000470162807048644,
     0.0005768904023010327]
    >>> # Spectrum for Twins-2-1 left postcentral pial surface (22)
    >>> # (after running explode_scalars() with reindex=True):
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_from_file
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'label22.vtk')
    >>> spectrum_from_file(vtk_file, spectrum_size=6)
    [6.3469513010430304e-18,
     0.0005178862383467463,
     0.0017434911095630772,
     0.003667561767487686,
     0.005429017880363784,
     0.006309346984678924]
    >>> # Loop thru all MB 101 brains
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_from_file
    >>> for hemidir in os.listdir(header):
    >>>     print hemidir
    >>>     sulci_file = os.path.join(header, hemidir, "sulci.vtk")
    >>>     spectrum = spectrum_from_file(sulci_file)

    """
    from mindboggle.utils.io_vtk import read_vtk, read_scalars
    from mindboggle.shapes.laplace_beltrami import spectrum_of_largest

    faces, u1, u2, points, u4, u5, u6, u7 = read_vtk(vtk_file)

    # Area file:
    if area_file:
        areas, u1 = read_scalars(area_file)
    else:
        areas = None

    spectrum = spectrum_of_largest(points, faces, spectrum_size,
                                   exclude_labels, normalization, areas)

    return spectrum


def spectrum_per_label(vtk_file, spectrum_size=10, exclude_labels=[-1],
                       normalization='area', area_file='',
                       largest_segment=True):
    """
    Compute Laplace-Beltrami spectrum per labeled region in a file.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file containing index scalars (labels)
    spectrum_size : integer
        number of eigenvalues to be computed (the length of the spectrum)
    exclude_labels : list of integers
        labels to be excluded
    normalization : string
        the method used to normalize eigenvalues ('area' or None)
        if "area", use area of the 2D structure as in Reuter et al. 2006
    area_file :  string
        name of VTK file with surface area scalar values
    largest_segment :  Boolean
        compute spectrum only for largest segment with a given label?

    Returns
    -------
    spectrum_lists : list of lists
        first eigenvalues for each label's Laplace-Beltrami spectrum
    label_list : list of integers
        list of unique labels for which spectra are obtained

    Examples
    --------
    >>> # Uncomment "if label==22:" below to run example:
    >>> # Spectrum for Twins-2-1 left postcentral (22) pial surface:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_per_label
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> spectrum_size = 6
    >>> exclude_labels = [0]  #[-1]
    >>> largest_segment = True
    >>> spectrum_per_label(vtk_file, spectrum_size, exclude_labels, None,
    >>>                    area_file, largest_segment)
    ([[6.3469513010430304e-18,
       0.0005178862383467463,
       0.0017434911095630772,
       0.003667561767487686,
       0.005429017880363784,
       0.006309346984678924]],
     [22])

    """
    from mindboggle.utils.io_vtk import read_vtk, read_scalars
    from mindboggle.utils.mesh import remove_faces, reindex_faces_points
    from mindboggle.shapes.laplace_beltrami import fem_laplacian,\
        spectrum_of_largest

    # Read VTK surface mesh file:
    faces, u1, u2, points, u4, labels, u5, u6 = read_vtk(vtk_file)

    # Area file:
    if area_file:
        areas, u1 = read_scalars(area_file)
    else:
        areas = None

    # Loop through labeled regions:
    ulabels = []
    [ulabels.append(int(x)) for x in labels if x not in ulabels
     if x not in exclude_labels]
    label_list = []
    spectrum_lists = []
    for label in ulabels:
      #if label == 22:
      #  print("DEBUG: COMPUTE FOR ONLY ONE LABEL")

        # Determine the indices per label:
        Ilabel = [i for i,x in enumerate(labels) if x == label]
        print('{0} vertices for label {1}'.format(len(Ilabel), label))

        # Remove background faces:
        pick_faces = remove_faces(faces, Ilabel)
        pick_faces, pick_points, o1 = reindex_faces_points(pick_faces, points)

        # Compute Laplace-Beltrami spectrum for the label:
        if largest_segment:
            exclude_labels_inner = [-1]
            spectrum = spectrum_of_largest(pick_points, pick_faces,
                                           spectrum_size,
                                           exclude_labels_inner,
                                           normalization, areas)
        else:
            spectrum = fem_laplacian(pick_points, pick_faces,
                                     spectrum_size, normalization)

        # Append to a list of lists of spectra:
        spectrum_lists.append(spectrum)
        label_list.append(label)

    return spectrum_lists, label_list


#if __name__ == "__main__":

    # import numpy as np
    # # You should get different outputs if you change the coordinates of points.
    # # If you do NOT see changes, you may be computing the graph Laplacian.
    #
    # # Define a cube:
    # points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1],
    #           [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    # # Pick some faces:
    # faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    #
    # print("Linear FEM Laplace-Beltrami spectrum\n\t{0}\n".format(
    #     fem_laplacian(points, faces, spectrum_size=5)))

