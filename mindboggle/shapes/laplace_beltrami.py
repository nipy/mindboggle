#!/usr/bin/python
"""
Compute the Laplace-Beltrami Spectrum using a linear finite element method.

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
        each list contains indices to vertices that form a triangle on a mesh

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

        1.5000-1.0000-0.5000    0    0    0    0    0
        -1.00002.0000    0    0-0.5000    0    0-0.5000
        -0.5000    02.0000-0.5000-1.0000    0    0    0
        0    0-0.50002.5607-0.3536-1.2071-0.5000    0
        0-0.5000-1.0000-0.35361.8536    0    0    0
        0    0    0-1.2071    01.2071    0    0
        0    0    0-0.5000    0    00.5000    0
        0-0.5000    0    0    0    0    00.5000
    >>> print_sparse_matrix(B)
        The sparse matrix:

        0.25000.08330.0417    00.0833    0    00.0417
        0.08330.1667    0    00.0417    0    00.0417
        0.0417    00.16670.04170.0833    0    0    0
        0    00.04170.28450.10060.10060.0417    0
        0.08330.04170.08330.10060.36790.0589    0    0
        0    0    00.10060.05890.20120.0417    0
        0    0    00.0417    00.04170.0833    0
        0.04170.0417    0    0    0    0    00.0833

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
        print("Truncation error of WESD is: {0}".format(W - WN))
        print("Truncation error of nWESD is: {1}".format(1 -  WN/W))

    return WESD

def fem_laplacian(points, faces, n_eigenvalues=200, normalization=None):
    """
    Linear FEM laplacian code after Martin Reuter's MATLAB code.

    Note ::
        Output for label 2 in Twins-2-1 left hemisphere:
        [4.829758648026223e-18, 0.00012841730024671977,
         0.0002715181572272744, 0.00032051508471594173,
         0.000470162807048644, 0.0005768904023010327]

        Martin Reuter's code:
         Creator: ./shapeDNA-tria
         File: /home/forrest/Dropbox/Share_Research/NYTX/label2.vtk
         User: Forrest Bao s.bao@ttu.edu Stony Brook
         Refine: 0
         Degree: 1
         Dimension: 2
         Elements: 290134
         DoF: 145069
         NumEW: 6

         Area: 110016
         Volume: 534346
         BLength: 0
         EulerChar: 2

         Time(pre) : 2
         Time(calcAB) : 0
         Time(calcEW) : 7
         Time(total ) : 9

        Eigenvalues:
        {-4.7207711983791511358e-18 ; 0.00012841730024672144738 ;
          0.00027151815722727096853 ; 0.00032051508471592313632 ;
          0.0004701628070486902353  ; 0.00057689040230097490998 }

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
    >>> print("Computing the un-normalized linear FEM Laplace-Beltrami Spectrum\n")
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=3)))
        The un-normalized linear FEM Laplace-Beltrami Spectrum is:
        [-2.74238008172841e-16, 4.583592135001265, 4.800000000000001]
    >>> print("Computing the area-normalized linear FEM Laplace-Beltrami Spectrum\n")
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
    from scipy.sparse.linalg import eigsh, lobpcg
    import scipy

    from mindboggle.shapes.laplace_beltrami import computeAB

    min_npoints = 10 * n_eigenvalues
    npoints = len(points)

    if npoints < min_npoints:
        print("The input size {0} should be much larger than n_eigenvalues {1}".\
            format(npoints, n_eigenvalues))
        return None
    elif npoints < n_eigenvalues:
        n_eigenvalues = npoints

    A, B = computeAB(points, faces)

    try :
        eigenvalues, eigenvectors = eigsh(A, k=n_eigenvalues, M=B, sigma=0)
        # Note: eigs is for nonsymmetric matrices while
        #       eigsh is for real-symmetric or complex-Hermitian matrices:
        spectrum = eigenvalues.tolist()
    except RuntimeError:
        # Initial eigenvector values:
        init_eigenvecs = scipy.rand(A.shape[0], n_eigenvalues)

        # maxiter = 40 forces lobpcg to use 20 iterations.
        # Strangely, largest=false finds largest eigenvalues
        # and largest=True gives the smallest eigenvalues:
        eigenvalues, eigenvectors =  lobpcg(A, init_eigenvecs, B=B,
                                            largest=True, maxiter = 40)
        # Extract the real parts:
        spectrum = [value.real for value in eigenvalues]
        # For some reason, the eigenvalues from lobpcg are not sorted:
        spectrum.sort()

    if normalization == "area":
        spectrum = area_normalize(points, faces, spectrum)

    return spectrum

def fem_laplacian_from_labels(vtk_file, n_eigenvalues=3, exclude_labels=[-1],
                              normalization=None, area_file=''):
    """
    Compute linear FEM Laplace-Beltrami spectra from each labeled region in a file.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file containing index scalars (labels)
    n_eigenvalues : integer
        number of eigenvalues to return
    exclude_labels : list of integers
        labels to be excluded
    normalization : string
        the method used to normalize eigenvalues (default: None)
        if "area", use area of the 2D structure as mentioned in Reuter et al. 2006
    area_file :  string
        name of VTK file with surface area scalar values

    Returns
    -------
    spectrum_lists : list of lists
        first eigenvalues for each label's Laplace-Beltrami spectrum
    label_list : list of integers
        list of unique labels for which spectra are obtained

    Examples
    --------
    >>> # Spectrum for labels:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> n_eigenvalues = 6
    >>> exclude_labels = [0]  #[-1]
    >>> print("Computing the un-normalized linear FEM Laplace-Beltrami Spectrum\n")
    >>> print("{0}".format(fem_laplacian_from_labels(label_file, n_eigenvalues,
    >>>                    exclude_labels, normalization=None, area_file='')))
        Results for label 2:
        ([[6.6934294195879822e-05, 0.055666737899630753, 0.063101464519212155,
           0.067440354878078423, 0.084502865477046815, 0.09518313048174197]],
          [2])
    >>> print("Computing the area-normalized linear FEM Laplace-Beltrami Spectrum\n")
    >>> print("{0}".format(fem_laplacian_from_labels(label_file, n_eigenvalues,
    >>>       exclude_labels, normalization="area", area_file=area_file)))
        ([[6.4125855895383511e-05, 0.046637028413482488, 0.066348957202302425,
           0.069927320749993999, 0.090186784531755951, 0.097848820601276823]],
          [2])
        ([[-1.6433030975833292e-21, 5.272910722504152e-08, 2.0515299130701199e-07,
           3.2441092366550171e-07, 6.4210373686882463e-07, 8.0270533114561409e-07]],
          [2])
    >>> # Spectrum for one label (artificial composite), two pieces:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> n_eigenvalues = 6
    >>> exclude_labels = [0]  #[-1]
    >>> #
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_vtk, write_vtk
    >>> faces, lines, indices, points, foo1, labels, foo2, foo3 = read_vtk(label_file,
    >>>      return_first=True, return_array=True)
    >>> I2 = [i for i,x in enumerate(labels) if x==2]
    >>> I11 = [i for i,x in enumerate(labels) if x==11]
    >>> scalars = np.zeros(np.shape(labels))
    >>> scalars[I2] = 1
    >>> scalars[I11] = 1
    >>> vtk_file = 'test_two_labels.vtk'
    >>> write_vtk(vtk_file, points, indices, lines, faces,
    >>>           scalars, scalar_names='scalars')
    >>> print("Computing the un-normalized linear FEM Laplace-Beltrami Spectrum\n")
    >>> print("{0}".format(fem_laplacian_from_labels(vtk_file,
    >>>       n_eigenvalues, exclude_labels, normalization=None,
    >>>       area_file=area_file)))
        Load "scalars" scalars from test_two_labels.vtk
        Load "scalars" scalars from lh.pial.area.vtk
        16647 vertices for label 1
        2 segments for label 1
        Reduced 290134 to 14498 triangular faces
        ([[-8.764053090852845e-18, 0.00028121452203987146, 0.0010941205613292243,
           0.0017301461686759188, 0.0034244633555606295, 0.004280982704174599]],
         [1])

    """
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, read_scalars
    from mindboggle.utils.mesh import reindex_faces_points, remove_faces
    from mindboggle.utils.segment import segment
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    min_npoints = 0#10 * n_eigenvalues

    # Read VTK surface mesh file:
    faces, foo1, foo2, points, npoints, labels, foo4, foo5 = read_vtk(vtk_file)
    neighbor_lists = find_neighbors(faces, npoints)

    if area_file:
        areas, foo1 = read_scalars(area_file, True, True)

    # Loop through labeled regions:
    ulabels = []
    [ulabels.append(int(x)) for x in labels if x not in ulabels
     if x not in exclude_labels]
    label_list = []
    spectrum_lists = []
    for label in ulabels:
      if label==2:

        # Determine the indices per label:
        label_indices = [i for i,x in enumerate(labels) if x == label]
        print('{0} vertices for label {1}'.format(len(label_indices), label))

        """
        # Segment the label indices into connected sets of indices:
        segment_label = segment(label_indices, neighbor_lists, min_region_size=1,
            seed_lists=[], keep_seeding=False, spread_within_labels=False,
            labels=[], label_lists=[], values=[], max_steps='', verbose=False)

        # Select the largest segment (connected set of indices for the label):
        select_indices = []
        max_segment_area = 0
        unique_segments = [x for x in np.unique(segment_label)
                           if x not in exclude_labels]
        if len(unique_segments) > 1:
            print('{0} segments for label {1}'.
                  format(len(unique_segments), label))
        for segment_number in unique_segments:
            if segment_number != -1:
                segment_indices = [i for i,x in enumerate(segment_label)
                                   if x == segment_number]
                if area_file:
                    segment_area = np.sum(areas[segment_indices])
                else:
                    segment_area = len(segment_indices)
                if segment_area > max_segment_area:
                    select_indices = segment_indices
        """
        select_indices = label_indices


        # Extract points and renumber faces for the selected indices:
        if len(select_indices) >= min_npoints:
            select_faces = remove_faces(faces, select_indices)
            if select_faces:
                select_faces, label_points = reindex_faces_points(select_faces,
                                                                  points)

                # Compute Laplace-Beltrami spectrum for the labeled region:
                spectrum = fem_laplacian(label_points, select_faces,
                                         n_eigenvalues, normalization)

                # Append to a list of lists of spectra:
                spectrum_lists.append(spectrum)
                label_list.append(label)
        else:
            print("The input size {0} for label {1} is too small. Skipped.".\
                format(len(select_indices), label))

    return spectrum_lists, label_list

def fem_laplacian_from_single_vtk(vtk_file, n_eigenvalues=6):
    """
    Compute linear FEM Laplace-Beltrami spectrum of a 3D shape in a VTK file.

    Note ::
        Output for label 2 in Twins-2-1 left hemisphere:
        [4.829758648026223e-18, 0.00012841730024671977,
         0.0002715181572272744, 0.00032051508471594173,
         0.000470162807048644, 0.0005768904023010327]
        Martin Reuter's code:
        {-4.7207711983791511358e-18 ; 0.00012841730024672144738 ;
          0.00027151815722727096853 ; 0.00032051508471592313632 ;
          0.0004701628070486902353  ; 0.00057689040230097490998 }

    Parameters
    ----------
    vtk_file : string
        the input vtk file
    n_eigenvalues : integer
        number of eigenvalues to be computed (the length of the spectrum)

    >>> # Spectrum for one label (artificial composite), two pieces:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_single_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> fem_laplacian_from_single_vtk(vtk_file, n_eigenvalues=6)
        Results for label 2:
        [6.4507181039686386e-05, 0.043252103515002013, 0.054016629188302107, 0.065206572702442234, 0.07756974245771403, 0.10151186954732085]


    """
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    faces, lines, indices, points, npoints, scalars, name, input_vtk = read_vtk(vtk_file)

    try:
        print("{0}".format(fem_laplacian(points, faces, n_eigenvalues)))
    except RuntimeError:
        print("Failed to compute Laplace-Beltrami spectrum " \
              "due to RuntimeError (e.g., singularity error)")


if __name__ == "__main__":

    import os
    from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels

    path = os.environ['MINDBOGGLE_DATA']
    label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    n_eigenvalues = 6
    exclude_labels = [0]  #[-1]
    print("Computing the un-normalized linear FEM Laplace-Beltrami Spectrum\n")
    print("{0}".format(fem_laplacian_from_labels(label_file, n_eigenvalues,
                       exclude_labels, normalization="area", area_file=area_file)))


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
    # print("Computing the un-normalized linear FEM Laplace-Beltrami Spectrum\n\t{0}\n".format(
    #     fem_laplacian(points, faces, n_eigenvalues=5)))

