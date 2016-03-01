#!/usr/bin/python
"""
Compute the Laplace-Beltrami spectrum using a linear finite element method.

We follow the definitions and steps given in Martin Reuter et al.'s 2009 paper:
    ``Discrete Laplace-Beltrami Operators for Shape Analysis and Segmentation''

References (please cite when using for publication):
Martin Reuter et al.
    Discrete Laplace-Beltrami Operators for Shape Analysis and Segmentation
    Computers & Graphics 33(3):381-390, 2009
Martin Reuter et al.
    Laplace-Beltrami spectra as "Shape-DNA" of surfaces and solids
    Computer-Aided Design 38(4):342-366, 2006

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
    - Martin Reuter, 2009-2016, http://reuter.mit.edu/ (original MATLAB code)
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)
    - Forrest Sheng Bao, 2012-2013  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
    A : csr_matrix
    B : csr_matrix

    Examples
    --------
    >>> # Define a cube, then compute A and B on a selection of faces.
    >>> import numpy as np
    >>> from mindboggle.shapes.laplace_beltrami import computeAB
    >>> points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1],
    ...           [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    >>> points = np.array(points)
    >>> faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    >>> faces = np.array(faces)
    >>> A, B = computeAB(points, faces)
    >>> print(np.array_str(A.toarray(), precision=5, suppress_small=True))
    [[ 1.5     -1.      -0.5      0.       0.       0.       0.       0.     ]
     [-1.       2.       0.       0.      -0.5      0.       0.      -0.5    ]
     [-0.5      0.       2.      -0.5     -1.       0.       0.       0.     ]
     [ 0.       0.      -0.5      2.56066 -0.35355 -1.20711 -0.5      0.     ]
     [ 0.      -0.5     -1.      -0.35355  1.85355  0.       0.       0.     ]
     [ 0.       0.       0.      -1.20711  0.       1.20711  0.       0.     ]
     [ 0.       0.       0.      -0.5      0.       0.       0.5      0.     ]
     [ 0.      -0.5      0.       0.       0.       0.       0.       0.5    ]]
    >>> print(np.array_str(B.toarray(), precision=5, suppress_small=True))
    [[ 0.25     0.08333  0.04167  0.       0.08333  0.       0.       0.04167]
     [ 0.08333  0.16667  0.       0.       0.04167  0.       0.       0.04167]
     [ 0.04167  0.       0.16667  0.04167  0.08333  0.       0.       0.     ]
     [ 0.       0.       0.04167  0.28452  0.10059  0.10059  0.04167  0.     ]
     [ 0.08333  0.04167  0.08333  0.10059  0.36785  0.05893  0.       0.     ]
     [ 0.       0.       0.       0.10059  0.05893  0.20118  0.04167  0.     ]
     [ 0.       0.       0.       0.04167  0.       0.04167  0.08333  0.     ]
     [ 0.04167  0.04167  0.       0.       0.       0.       0.       0.08333]]

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
    vol_mean = 0.001*np.mean(vol)
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
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    spectrum : list of floats
        LB spectrum of a given shape defined by _points_ and _faces_

    Returns
    -------
    new_spectrum : list of floats
        LB spectrum normalized by area

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.shapes.laplace_beltrami import area_normalize
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> # Define a cube:
    >>> points = [[0,0,0], [0,1,0], [1,1,0], [1,0,0],
    ...           [0,0,1], [0,1,1], [1,1,1], [1,0,1]]
    >>> faces = [[0,1,2], [2,3,0], [4,5,6], [6,7,4], [0,4,7], [7,3,0],
    ...          [0,4,5], [5,1,0], [1,5,6], [6,2,1], [3,7,6], [6,2,3]]
    >>> spectrum = fem_laplacian(points, faces, spectrum_size=3,
    ...                          normalization=None)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...       precision=5, suppress_small=True))
    [ 4.58359  4.8    ]
    >>> new_spectrum = area_normalize(points, faces, spectrum)
    >>> print(np.array_str(np.array(new_spectrum[1::]),
    ...       precision=5, suppress_small=True))
    [ 27.50155  28.8    ]

    """
    from mindboggle.guts.mesh import area_of_faces

    area = area_of_faces(points, faces)
    total_area = sum(area)

    new_spectrum = [x*total_area for x in spectrum]

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
        display the error of WESD using Eqs.(9) and (10) in Konukoglu (2012)?
        default: false
    N : integer
        The length of spetrum used (N>=3, default: 3)

    """

    # Martin Reuter: "a surface is a 2d manifold.
    # It doesn't matter that it is usually embedded in 3d Euclidean space."
    d = 2.0
    Ball = 4.0 / 3 * np.pi  # For three dimensions
    p = 2.0

    Vol = np.amax((Vol1, Vol2))
    mu = np.amax(EVAL1[1], EVAL2[1])

    C = ((d+2)/(d*4*np.pi**2)*(Ball*Vol)**(2/d) - 1/mu)**p + \
        ((d+2)/(d*4*np.pi**2)*(Ball*Vol/2)**(2/d) - 1/mu*(d/(d+4)))**p

    K = ((d+2)/(d*4*np.pi**2)*(Ball*Vol)**(2/d) - (1/mu)*(d/(d+2.64)))**p

    # the right-hand side of Eq.(8) or the equation right after Eq.(4):
    W = (C + K*(zeta(2*p/d,1) - 1 - .5**(2*p/d)))**(1/p)

    holder = 0
    for i in xrange(1, np.amin((len(EVAL1), len(EVAL2) )) ):
        holder += (np.abs(EVAL1[i] - EVAL2[i])/(EVAL1[i]*EVAL2[i]))**p
    WESD = holder ** (1/p)

    #nWESD = WESD/W

    if show_error:
        WN = (C + K * (sum([n**(-1*2*p/d) for n in range(3,N+1)])))**(1/p)
        # the second term on the right-hand side of Eq.(9)
        #print("Truncation error of WESD is: {0}".format(W - WN))
        #print("Truncation error of nWESD is: {1}".format(1 -  WN/W))

    return WESD


def fem_laplacian(points, faces, spectrum_size=10, normalization=None,
                  verbose=False):
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

    Julien Lefevre, regarding comparison with Spongy results:
    "I have done some comparisons between my Matlab codes and yours
    on python and everything sounds perfect:
    The evaluation has been done only for one mesh (about 10000 vertices)."
    - L2 error between your A and B matrices and mine are about 1e-16.
    - I have also compared eigenvalues of the generalized problem;
      even if the error is slightly increasing, it remains on the order
      of machine precision.
    - Another good point: computation time for 1000 eigenvalues was 67s
      with python versus 63s in matlab.
      And it is quite the same behavior for lower orders.
    - Since the eigenvalues are increasing with order,
      it is also interesting to look at the relative error...
      high frequencies are not so much perturbed."

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
    verbose : Boolean
        print statements?

    Returns
    -------
    spectrum : list
        first spectrum_size eigenvalues for Laplace-Beltrami spectrum

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> # Define a cube:
    >>> points = [[0,0,0], [0,1,0], [1,1,0], [1,0,0],
    ...           [0,0,1], [0,1,1], [1,1,1], [1,0,1]]
    >>> faces = [[0,1,2], [2,3,0], [4,5,6], [6,7,4], [0,4,7], [7,3,0],
    ...          [0,4,5], [5,1,0], [1,5,6], [6,2,1], [3,7,6], [6,2,3]]
    >>> spectrum = fem_laplacian(points, faces, spectrum_size=3,
    ...                          normalization=None, verbose=False)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...                    precision=5, suppress_small=True))
    [ 4.58359  4.8    ]
    >>> spectrum = fem_laplacian(points, faces, spectrum_size=3,
    ...                          normalization="area", verbose=False)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...                    precision=5, suppress_small=True))
    [ 27.50155  28.8    ]
    >>> # Spectrum for entire left hemisphere of Twins-2-1:
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['left_freesurfer_labels'])
    >>> points, f1,f2, faces, labels, f3,f4,f5 = read_vtk(label_file)
    >>> spectrum = fem_laplacian(points, faces, spectrum_size=6,
    ...                          normalization=None, verbose=False)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...                    precision=5, suppress_small=True))
    [ 0.00013  0.00027  0.00032  0.00047  0.00058]
    >>> # Spectrum for Twins-2-1 left postcentral pial surface (22):
    >>> from mindboggle.guts.mesh import keep_faces, reindex_faces_points
    >>> I22 = [i for i,x in enumerate(labels) if x==1022] # postcentral
    >>> faces = keep_faces(faces, I22)
    >>> faces, points, o1 = reindex_faces_points(faces, points)
    >>> spectrum = fem_laplacian(points, faces, spectrum_size=6,
    ...                          normalization=None, verbose=False)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...                    precision=5, suppress_small=True))
    [ 0.00057  0.00189  0.00432  0.00691  0.00775]
    >>> # Area-normalized spectrum for a single label (postcentral):
    >>> spectrum = fem_laplacian(points, faces, spectrum_size=6,
    ...                          normalization="area", verbose=False)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...                    precision=5, suppress_small=True))
    [  2.69259   8.97865  20.44857  32.74477  36.739  ]

    """
    from scipy.sparse.linalg import eigsh, lobpcg
    import numpy as np

    from mindboggle.shapes.laplace_beltrami import computeAB

    #-----------------------------------------------------------------
    # Compute A and B matrices (from Reuter et al., 2009):
    #-----------------------------------------------------------------
    A, B = computeAB(points, faces)
    if A.shape[0] <= spectrum_size:
        if verbose:
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
                                          sigma=-0.01)
        spectrum = eigenvalues.tolist()

    #-----------------------------------------------------------------
    # Use the lobpcg eigensolver:
    #-----------------------------------------------------------------
    except RuntimeError:     
           
        if verbose:
            print("eigsh() failed. Now try lobpcg.")
            print("Warning: lobpcg can produce different results from "
                  "Reuter (2006) shapeDNA-tria software.")
        # Initial eigenvector values:
        init_eigenvecs = np.random.random((A.shape[0], spectrum_size))

        # maxiter = 40 forces lobpcg to use 20 iterations.
        # Strangely, largest=false finds largest eigenvalues
        # and largest=True gives the smallest eigenvalues:
        eigenvalues, eigenvectors = lobpcg(A, init_eigenvecs, B=B,
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
        if verbose:
            print("Compute area-normalized linear FEM Laplace-Beltrami "
                  "spectrum")
    else:
        if verbose:
            print("Compute linear FEM Laplace-Beltrami spectrum")

    return spectrum


def spectrum_of_largest(points, faces, spectrum_size=10, exclude_labels=[-1],
                        normalization=None, areas=None, verbose=False):
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
    verbose : Boolean
        print statements?

    Returns
    -------
    spectrum : list
        first spectrum_size eigenvalues for Laplace-Beltrami spectrum

    Examples
    --------
    >>> # Spectrum for left postcentral + pars triangularis pial surfaces:
    >>> import numpy as np
    >>> from mindboggle.mio.vtks import read_scalars, read_vtk, write_vtk
    >>> from mindboggle.guts.mesh import keep_faces, reindex_faces_points
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_of_largest
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['left_freesurfer_labels'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> spectrum_size = 6
    >>> exclude_labels = [-1]
    >>> normalization = None
    >>> points, indices, lines, faces, labels, f1, npoints, f2 = read_vtk(label_file,
    ...     return_first=True, return_array=True)
    >>> I20 = [i for i,x in enumerate(labels) if x==1020] # pars triangularis
    >>> I22 = [i for i,x in enumerate(labels) if x==1022] # postcentral
    >>> I22.extend(I20)
    >>> faces = keep_faces(faces, I22)
    >>> faces, points, o1 = reindex_faces_points(faces, points)
    >>> areas, u1 = read_scalars(area_file, True, True)
    >>> verbose = False
    >>> spectrum = spectrum_of_largest(points, faces, spectrum_size,
    ...     exclude_labels, normalization, areas, verbose)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...                    precision=5, suppress_small=True))
    [ 0.00057  0.00189  0.00432  0.00691  0.00775]

    View both segments (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces
    >>> scalars = np.zeros(np.shape(labels))
    >>> scalars[I22] = 1
    >>> vtk_file = 'test_two_labels.vtk'
    >>> write_vtk(vtk_file, points, indices, lines, faces,
    ...           scalars, scalar_names='scalars', scalar_type='int')
    >>> plot_surfaces(vtk_file) # doctest: +SKIP

    """
    import numpy as np
    #from scipy.sparse.linalg import eigsh, lobpcg

    from mindboggle.guts.segment import select_largest
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    if isinstance(areas, list):
        areas = np.array(areas)

    # Check to see if there are enough points:
    min_points_faces = spectrum_size
    npoints = len(points) 
    if npoints < min_points_faces or len(faces) < min_points_faces:
        raise IOError("The input size {0} ({1} faces) should be much larger "
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
            raise IOError("The input size {0} is too small.".
                          format(len(points)))
            return None
        elif faces:

            #-----------------------------------------------------------------
            # Compute spectrum:
            #-----------------------------------------------------------------
            spectrum = fem_laplacian(points, faces, spectrum_size,
                                     normalization, verbose)
            return spectrum
        else:
            return None


def spectrum_from_file(vtk_file, spectrum_size=10, exclude_labels=[-1],
                       normalization=None, area_file='', verbose=False):
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
    verbose : Boolean
        print statements?

    Returns
    -------
    spectrum : list of floats
        first spectrum_size of Laplace-Beltrami spectrum

    Examples
    --------
    >>> # Spectrum for entire left hemisphere of Twins-2-1:
    >>> import numpy as np
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_from_file
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_per_label
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> spectrum = spectrum_from_file(vtk_file, spectrum_size=6)
    >>> print(np.array_str(np.array(spectrum[1::]),
    ...                    precision=5, suppress_small=True))
    [ 0.00013  0.00027  0.00032  0.00047  0.00058]

    """
    from mindboggle.mio.vtks import read_vtk, read_scalars
    from mindboggle.shapes.laplace_beltrami import spectrum_of_largest

    points, indices, lines, faces, scalars, scalar_names, npoints, \
            input_vtk = read_vtk(vtk_file)

    # Area file:
    if area_file:
        areas, u1 = read_scalars(area_file)
    else:
        areas = None

    spectrum = spectrum_of_largest(points, faces, spectrum_size,
                                   exclude_labels, normalization, areas,
                                   verbose)

    return spectrum


def spectrum_per_label(vtk_file, spectrum_size=10, exclude_labels=[-1],
                       normalization='area', area_file='',
                       largest_segment=True, verbose=False):
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
    area_file :  string (optional)
        name of VTK file with surface area scalar values
    largest_segment :  Boolean
        compute spectrum only for largest segment with a given label?
    verbose : Boolean
        print statements?

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
    >>> import numpy as np
    >>> from mindboggle.shapes.laplace_beltrami import spectrum_per_label
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> spectrum_size = 6
    >>> exclude_labels = [0]  #[-1]
    >>> largest_segment = True
    >>> verbose = False
    >>> spectrum_lists, label_list = spectrum_per_label(vtk_file,
    ...     spectrum_size, exclude_labels, None, area_file, largest_segment,
    ...     verbose)
    >>> print(np.array_str(np.array(spectrum_lists[0][1::]),
    ...                    precision=5, suppress_small=True))
    [ 0.00054  0.00244  0.00291  0.00456  0.00575]
    >>> label_list[0:10]
    [1029, 1005, 1011, 1021, 1008, 1025, 999, 1013, 1007, 1022]

    """
    from mindboggle.mio.vtks import read_vtk, read_scalars
    from mindboggle.guts.mesh import keep_faces, reindex_faces_points
    from mindboggle.shapes.laplace_beltrami import fem_laplacian,\
        spectrum_of_largest

    # Read VTK surface mesh file:
    points, indices, lines, faces, labels, scalar_names, npoints, \
        input_vtk = read_vtk(vtk_file)

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
        if verbose:
          print('{0} vertices for label {1}'.format(len(Ilabel), label))

        # Remove background faces:
        pick_faces = keep_faces(faces, Ilabel)
        pick_faces, pick_points, o1 = reindex_faces_points(pick_faces, points)

        # Compute Laplace-Beltrami spectrum for the label:
        if largest_segment:
            exclude_labels_inner = [-1]
            spectrum = spectrum_of_largest(pick_points, pick_faces,
                                           spectrum_size,
                                           exclude_labels_inner,
                                           normalization, areas, verbose)
        else:
            spectrum = fem_laplacian(pick_points, pick_faces, spectrum_size,
                                     normalization, verbose)

        # Append to a list of lists of spectra:
        spectrum_lists.append(spectrum)
        label_list.append(label)

    return spectrum_lists, label_list


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()


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
