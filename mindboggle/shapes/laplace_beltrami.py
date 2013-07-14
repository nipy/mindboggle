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

def fem_laplacian(points, faces, n_eigenvalues=6, background=-1,
                  normalization=None, areas=None):
    """
    Linear FEM laplacian code after Martin Reuter's MATLAB code.

    Note ::

        Compare fem_laplacian with Martin Reuter's Matlab code output:

        fem_laplacian() results for Twins-2-1 left hemisphere:
        [4.829758648026223e-18,
         0.00012841730024671977,
         0.0002715181572272744,
         0.00032051508471594173,
         0.000470162807048644,
         0.0005768904023010327]

        Martin Reuter's Matlab code:
         Creator: ./shapeDNA-tria
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
        {-4.7207711983791511358e-18 ;
         0.00012841730024672144738 ;
         0.00027151815722727096853 ;
         0.00032051508471592313632 ;
         0.0004701628070486902353  ;
         0.00057689040230097490998 }

        fem_laplacian() results for Twins-2-1 left hemisphere postcentral:
        [6.346951301043029e-18,
         0.0005178862383467465,
         0.0017434911095630787,
         0.0036675617674876916,
         0.005429017880363785,
         0.006309346984678927]

        Martin Reuter's Matlab code:
         -2.1954862991027e-18
         0.0005178862383468
         0.0017434911095628
         0.0036675617674875
         0.0054290178803611
         0.006309346984678

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    n_eigenvalues : integer
        number of eigenvalues to return
    background : integer
        background value
    normalization : string
        the method used to normalize eigenvalues ('area' or None)
        if "area", use area of the 2D structure as in Reuter et al. 2006
    areas : numpy array or list of floats (or None)
        surface area scalar values for all vertices

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
    >>> print("Linear FEM Laplace-Beltrami spectrum")
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=3)))
        [7.401486830834377e-17, 4.58359213500127, 4.799999999999998]
    >>> print("Area-normalized linear FEM Laplace-Beltrami spectrum")
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=3, 
    >>>                                  normalization="area")))
        [1.2335811384723967e-17, 0.76393202250021175, 0.79999999999999949]
    >>> # Spectrum for entire left hemisphere of Twins-2-1:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_faces_points
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> fold_file = os.path.join(path, 'arno', 'labels',
    >>>                          'lh.labels.DKT25.manual.vtk')
    >>> faces, points, npoints = read_faces_points(fold_file)
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=6,
    >>>                    normalization=None)))
        [4.829758648026223e-18,
         0.00012841730024671904,
         0.00027151815722727406,
         0.00032051508471594146,
         0.0004701628070486449,
         0.0005768904023010303]
    >>> # Spectrum for Twins-2-1 left hemisphere postcentral (label 22):
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_faces_points
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'label22.vtk')
    >>> faces, points, npoints = read_faces_points(label_file)
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=6,
    >>>                    normalization=None)))
        [6.346951301043029e-18,
         0.0005178862383467465,
         0.0017434911095630787,
         0.0036675617674876916,
         0.005429017880363785,
         0.006309346984678927]
    >>> # Area-normalized spectrum for a single label (postcentral):
    >>> print("{0}".format(fem_laplacian(points, faces, n_eigenvalues=6,
    >>>                    normalization="area")))
        [1.1410192787181146e-21,
         9.310268097367214e-08,
         3.1343504525679715e-07,
         6.593336681038091e-07,
         9.759983608165455e-07,
         1.1342589857996225e-06]
    """
    from scipy.sparse.linalg import eigsh, lobpcg
    import numpy as np

    from mindboggle.utils.mesh import find_neighbors, remove_faces, \
        reindex_faces_points
    from mindboggle.utils.segment import segment
    from mindboggle.shapes.laplace_beltrami import computeAB

    # Check areas type:
    use_area = False
    if areas:
        if isinstance(areas, np.ndarray):
            use_area = True
        elif isinstance(areas, list):
            use_area = True
            areas = np.array(areas)

    # Check to see if there are enough points:
    min_npoints = n_eigenvalues
    npoints = len(points)
    if npoints < min_npoints or len(faces) < min_npoints:
        print("The input size {0} ({1} faces) should be much larger "
              "than n_eigenvalues {2}".\
              format(npoints, len(faces), n_eigenvalues))
        return None
    else:

        #---------------------------------------------------------------------
        # Segment the indices into connected sets of indices:
        #---------------------------------------------------------------------
        # Construct neighbor lists:
        neighbor_lists = find_neighbors(faces, npoints)

        # Determine the indices:
        indices = [x for sublst in faces for x in sublst]

        # Segment:
        segments = segment(indices, neighbor_lists, min_region_size=1,
            seed_lists=[], keep_seeding=False, spread_within_labels=False,
            labels=[], label_lists=[], values=[], max_steps='', verbose=False)

        #---------------------------------------------------------------------
        # Select the largest segment (connected set of indices):
        #---------------------------------------------------------------------
        unique_segments = [x for x in np.unique(segments) if x != background]
        if len(unique_segments) > 1:
            select_indices = []
            max_segment_area = 0
            print('{0} segments'.format(len(unique_segments)))
            for segment_number in unique_segments:
                segment_indices = [i for i,x in enumerate(segments)
                                   if x == segment_number]
                if use_area:
                    segment_area = np.sum(areas[segment_indices])
                else:
                    segment_area = len(segment_indices)
                if segment_area > max_segment_area:
                    select_indices = segment_indices

            #-----------------------------------------------------------------
            # Extract points and renumber faces for the selected indices:
            #-----------------------------------------------------------------
            faces = remove_faces(faces, select_indices)
        else:
            select_indices = indices

        # Alert if the number of indices is small:
        if len(select_indices) < min_npoints:
            print("The input size {0} is too small.".format(len(select_indices)))
            return None
        elif faces:

            #-----------------------------------------------------------------
            # Reindex indices in faces:
            #-----------------------------------------------------------------
            faces, points = reindex_faces_points(faces, points)

            #-----------------------------------------------------------------
            # Compute A and B matrices (from Reuter's et al., 2009):
            #-----------------------------------------------------------------
            A, B = computeAB(points, faces)

            #-----------------------------------------------------------------
            # Use the eigsh eigensolver:
            #-----------------------------------------------------------------
            try :

                # eigs is for nonsymmetric matrices while
                # eigsh is for real-symmetric or complex-Hermitian matrices:
                eigenvalues, eigenvectors = eigsh(A, k=n_eigenvalues, M=B,
                                                  sigma=0)
                spectrum = eigenvalues.tolist()

            #-----------------------------------------------------------------
            # Use the lobpcg eigensolver:
            #-----------------------------------------------------------------
            except RuntimeError:

                # Initial eigenvector values:
                init_eigenvecs = np.random.random((A.shape[0], n_eigenvalues))

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
                print("Area-normalized linear FEM Laplace-Beltrami spectrum:")
            else:
                print("Linear FEM Laplace-Beltrami spectrum:")

            return spectrum

        else:
            return None


def fem_laplacian_from_labels(vtk_file, n_eigenvalues=3, exclude_labels=[-1],
                              normalization='area', area_file=''):
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
        the method used to normalize eigenvalues ('area' or None)
        if "area", use area of the 2D structure as in Reuter et al. 2006
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
    >>> # Spectrum for label 22 (postcentral) in Twins-2-1:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> n_eigenvalues = 6
    >>> exclude_labels = [0]  #[-1]
    >>> fem_laplacian_from_labels(vtk_file, n_eigenvalues, exclude_labels,
    >>>                           normalization=None, area_file=area_file)
        Load "Labels" scalars from lh.labels.DKT25.manual.vtk
        Load "scalars" scalars from lh.pial.area.vtk
        7819 vertices for label 22
        Reduced 290134 to 15230 triangular faces
        Linear FEM Laplace-Beltrami spectrum:
        ([[6.3469513010430304e-18,
           0.0005178862383467463,
           0.0017434911095630772,
           0.003667561767487686,
           0.005429017880363784,
           0.006309346984678924]],
         [22])
    >>> # Spectrum for one label (artificial composite), two fragments:
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
    >>> I2 = [i for i,x in enumerate(labels) if x==2] # cingulate
    >>> I22 = [i for i,x in enumerate(labels) if x==22] # postcentral
    >>> scalars = np.zeros(np.shape(labels))
    >>> scalars[I2] = 1
    >>> scalars[I22] = 1
    >>> vtk_file = 'test_two_labels.vtk'
    >>> write_vtk(vtk_file, points, indices, lines, faces,
    >>>           scalars, scalar_names='scalars')
    >>> fem_laplacian_from_labels(vtk_file, n_eigenvalues, exclude_labels,
    >>>                           normalization=None, area_file=area_file)
        Load "scalars" scalars from test_two_labels.vtk
        Load "scalars" scalars from lh.pial.area.vtk
        15357 vertices for label 1
        Reduced 290134 to 29728 triangular faces
        2 segments
        Reduced 29728 to 14498 triangular faces
        Linear FEM Laplace-Beltrami spectrum:
        ([[-8.764053090852845e-18,
           0.00028121452203987146,
           0.0010941205613292243,
           0.0017301461686759188,
           0.0034244633555606295,
           0.004280982704174599]],
         [1])

    """
    from mindboggle.utils.io_vtk import read_vtk, read_scalars
    from mindboggle.utils.mesh import remove_faces
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

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
      #if label==22:

        # Determine the indices per label:
        label_indices = [i for i,x in enumerate(labels) if x == label]
        print('{0} vertices for label {1}'.format(len(label_indices), label))

        # Remove background faces:
        select_faces = remove_faces(faces, label_indices)

        # Compute Laplace-Beltrami spectrum for the label:
        spectrum = fem_laplacian(points, select_faces, n_eigenvalues,
            background=-1, normalization=normalization, areas=areas)

        # Append to a list of lists of spectra:
        spectrum_lists.append(spectrum)
        label_list.append(label)

    return spectrum_lists, label_list

def fem_laplacian_from_file(vtk_file, n_eigenvalues=6, normalization=None):
    """
    Compute linear FEM Laplace-Beltrami spectrum of a 3D shape in a VTK file.

    Parameters
    ----------
    vtk_file : string
        the input vtk file
    n_eigenvalues : integer
        number of eigenvalues to be computed (the length of the spectrum)
    normalization : string
        the method used to normalize eigenvalues ('area' or None)
        if "area", use area of the 2D structure as in Reuter et al. 2006

    Returns
    -------
    spectrum : list of floats
        first n_eigenvalues of Laplace-Beltrami spectrum

    Examples
    --------
    >>> # Spectrum for entire left hemisphere of Twins-2-1:
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_file
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> fem_laplacian_from_file(vtk_file, n_eigenvalues=6)
        Load "Labels" scalars from lh.labels.DKT25.manual.vtk
        Linear FEM Laplace-Beltrami spectrum:
        [4.829758648026221e-18,
         0.00012841730024672036,
         0.00027151815722727465,
         0.00032051508471594065,
         0.0004701628070486447,
         0.0005768904023010338]
    >>> # Spectrum for label 22 (postcentral) (after running explode_scalars()):
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'label22.vtk')
    >>> fem_laplacian_from_file(vtk_file, n_eigenvalues=6)
        Load "scalars" scalars from label22.vtk
        Linear FEM Laplace-Beltrami spectrum:
        [6.3469513010430304e-18,
         0.0005178862383467466,
         0.0017434911095630806,
         0.003667561767487689,
         0.005429017880363778,
         0.006309346984678918]

    """
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    faces, u1, u2, points, u4, u5, u6, u7 = read_vtk(vtk_file)

    spectrum = fem_laplacian(points, faces, n_eigenvalues)
    return spectrum


if __name__ == "__main__":

    import os
    from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels

    path = os.environ['MINDBOGGLE_DATA']
    label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    n_eigenvalues = 6
    exclude_labels = [0]  #[-1]
    print("Linear FEM Laplace-Beltrami spectrum")
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
    # print("Linear FEM Laplace-Beltrami spectrum\n\t{0}\n".format(
    #     fem_laplacian(points, faces, n_eigenvalues=5)))

