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
    from scipy.sparse.linalg import eigsh, lobpcg
    import scipy

    from mindboggle.shapes.laplace_beltrami import computeAB

    min_n_eigenvalues = 0# 10 * n_eigenvalues
    npoints = len(points)

    if npoints < min_n_eigenvalues:
        print "The input size {0} should be much larger than n_eigenvalue {1}. Skipped.".\
            format(npoints, n_eigenvalues)
        return None
    elif npoints < n_eigenvalues:
        n_eigenvalues = npoints

    A, B = computeAB(points, faces)

    try :
        eigenvalues, eigenvectors = eigsh(A, k=n_eigenvalues, M=B, sigma=0)
        # Note: eigs is for nonsymmetric matrices while
        #       eigsh is for real-symmetric or complex-Hermitian matrices.
        spectrum = eigenvalues.tolist()
    except RuntimeError:
        init_eigenvecs = scipy.rand(A.shape[0], n_eigenvalues) # initial eigenvector values 
            
        eigenvalues, eigenvectors =  lobpcg(A, init_eigenvecs, B=B, largest=True, maxiter = 40)
        # maxiter = 40 (forty) forces lobpcg to use 20 (twenty) iterations
        # Strangely, largest=false finds largest eigenvalues. 
        # I had to use largest=True to find smallest eigenvalues. 
        spectrum = [value.real for value in eigenvalues] # take the real parts out
        spectrum.sort() # for some reason, the eigenvalues from lobpcg is not sorted

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
    >>> print("The un-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian_from_labels(label_file, n_eigenvalues,
    >>>                    exclude_labels, normalization=None, area_file='')))
        ([[3.991471403036072e-17, 0.0008361843957687868, 0.0017871320208401006,
           0.0029994050521529613, 0.004667334518265755, 0.00513033437385627],
          [1.2391703345286005e-17, 0.0006134752952690859, 0.002609362442309575,
           0.004663168980214244, 0.00634652760526998, 0.00877724055345246],
          [1.0012068644990291e-17, 0.0022387207075428577, 0.006970152884917594,
           0.0091248802802332, 0.01564833626724935, 0.017910676892526053],
          [-2.1056485710445128e-17, 0.0009694328464937241, 0.00317818959971316,
           0.005054405066976481, 0.00789603110003099, 0.009229553258157068],
          [1.2368114019129475e-17, 0.003031104242333413, 0.009741447008073946,
           0.013007273963336891, 0.021752446526208086, 0.023014192373524967],
          [-4.974527972363915e-18, 0.0019016559845126855, 0.0025627314943418114,
           0.005959708974411091, 0.007055750503872705, 0.007881159600790508],
          [-6.864064879148695e-17, 0.0016543895024567043, 0.004288905710034005,
           0.006626629502359926, 0.0088839013872855, 0.013496308832233302],
          [6.346951301043029e-18, 0.0005178862383467466, 0.0017434911095630813,
           0.003667561767487689, 0.0054290178803637805, 0.006309346984678919],
          [-5.577706348301915e-17, 0.0013626966554686451, 0.0032593611935484747,
           0.004341403383623896, 0.008145164053283511, 0.01170208132644808],
          [9.328607502554145e-18, 0.0016340116861423826, 0.002377126729999986,
           0.004770252392447428, 0.005855207512818139, 0.007954700266141297],
          [3.1408972748645767e-17, 0.0032741466797619554, 0.008421528264955687,
           0.012592850550534159, 0.014793728473842147, 0.023426822648660153],
          [1.086715131229174e-17, 0.000396383812027545, 0.001584393280625915,
           0.002841274919883915, 0.004194276483443616, 0.005757938406357442],
          [3.605193776571674e-18, 0.0005807142594621487, 0.0020007773881150397,
           0.0042442307604486995, 0.005437927931201479, 0.006537436919971014],
          [2.8256444624793164e-17, 0.00045792420775603803, 0.00152199799445477,
           0.0024395035272770653, 0.004194506157056832, 0.006292423975821009],
          [-8.764053090852845e-18, 0.0002812145220398717, 0.0010941205613292206,
           0.0017301461686759227, 0.003424463355560629, 0.004280982704174603],
          [2.9604028746491564e-17, 0.0004528887382893767, 0.0015515885375304076,
           0.00331602414788001, 0.006619976369902932, 0.008205020906938507],
          [5.554957412505596e-18, 0.00035601986106019174, 0.001488072212933648,
           0.002028217741427394, 0.0021391981684121676, 0.0034678968678893123],
          [1.4259486522787315e-18, 0.009618295295697915, 0.025378077629970915,
           0.03248964988496929, 0.0508980016579518, 0.06719346053046973],
          [-2.4117396362925674e-17, 0.002615463222087654, 0.009472896293593986,
           0.012242050075865145, 0.01872759575825019, 0.020930902640074314],
          [-7.335755494644984e-18, 0.007912906429265811, 0.01652965473212952,
           0.023279981681191032, 0.037064614395027944, 0.048697661017445756],
          [9.499083277797273e-18, 0.0004338925675937878, 0.0016892631759901167,
           0.0031754688159043, 0.003511887898901262, 0.0047277285445355335],
          [1.9593383259665427e-17, 0.0014920483370167207, 0.0029582214145706227,
           0.004359091649425157, 0.005742293766092156, 0.009083586748755123],
          [1.7582836994280824e-17, 0.0015257857600227214, 0.004374858499594681,
           0.0054584077589882075, 0.00849735356103569, 0.011539236288697796],
          [4.852325972754763e-18, 0.0021908963589592333, 0.008674907947474586,
           0.014228894730428161, 0.017777570383697067, 0.018289655756614703]],
         [11, 29, 5, 25, 21, 8, 13, 22, 7, 31, 17, 30, 15, 24, 2, 9, 28,
          34, 35, 16, 3, 18, 12, 14])
    >>> print("The area-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
    >>> print("{0}".format(fem_laplacian_from_labels(label_file, n_eigenvalues,
    >>>       exclude_labels, normalization="area", area_file=area_file)))
        ([[5.9892025481818609e-21, 1.2546946246637232e-07, 2.6815914664981521e-07,
           4.5006070612756651e-07, 7.0033351031280842e-07, 7.6980663525618909e-07],
          [2.7442346880662513e-21, 1.3585865789706364e-07, 5.7786268186021908e-07,
           1.0326933848595218e-06, 1.4054856477642631e-06, 1.9437850730544671e-06],
          [4.8805153037603019e-21, 1.0912940233859737e-06, 3.3976932271045499e-06,
           4.4480435993556758e-06, 7.627988514532606e-06, 8.7307963792695699e-06],
          [-5.1999639167835664e-21, 2.3940442345569276e-07, 7.8486369788702351e-07,
           1.2482008788412443e-06, 1.9499491686589686e-06, 2.2792665675755041e-06],
          [8.4014236052977723e-21, 2.0589712135796223e-06, 6.6171788776208509e-06,
           8.8355927465686025e-06, 1.4776021423767969e-05, 1.5633101276778638e-05],
          [-1.1006988107353173e-21, 4.20773687917598e-07, 5.6704787343180347e-07,
           1.3186868416272294e-06, 1.5612046472757615e-06, 1.743840430287706e-06],
          [-2.0046182128355467e-20, 4.8315675713134583e-07, 1.2525549584453889e-06,
           1.9352763157144804e-06, 2.5945020677304866e-06, 3.9415342027628951e-06],
          [1.1410192787181146e-21, 9.3102680973672023e-08, 3.1343504525679774e-07,
           6.5933366810380837e-07, 9.7599836081654615e-07, 1.1342589857996225e-06],
          [-1.5130736911650156e-20, 3.6966099139584188e-07, 8.8417233966858203e-07,
           1.1776997267875037e-06, 2.209552219076263e-06, 3.174443092062819e-06],
          [1.8271159943703549e-21, 3.2004014381797145e-07, 4.6558784554276701e-07,
           9.3430926759838256e-07, 1.1468103137685179e-06, 1.5580203243313131e-06],
          [1.708810453193922e-20, 1.7813050163853839e-06, 4.5817466385126638e-06,
           6.8511615545240659e-06, 8.0485528960522633e-06, 1.2745402324205791e-05],
          [1.6992092565908249e-21, 6.197935624565444e-08, 2.4773886468983949e-07,
           4.442673618545527e-07, 6.5582536035063397e-07, 9.0032262897599203e-07],
          [6.6897355001613741e-22, 1.0775633815356515e-07, 3.7126080734338119e-07,
           7.8755215249623579e-07, 1.0090525442788043e-06, 1.2130755391790322e-06],
          [4.7563412979122466e-21, 7.708131187717483e-08, 2.5619436601068196e-07,
           4.1063592844973614e-07, 7.0605142027144236e-07, 1.0591890245777766e-06],
          [-1.6433030975833292e-21, 5.2729107225041559e-08, 2.0515299130701178e-07,
           3.2441092366550234e-07, 6.4210373686882643e-07, 8.0270533114561409e-07],
          [6.9244702442049931e-21, 1.0593202091090364e-07, 3.6292116695948432e-07,
           7.7562789638155606e-07, 1.5484321334529694e-06, 1.919178758057007e-06],
          [5.8885618235010395e-22, 3.7740072633636524e-08, 1.5774387763922336e-07,
           2.1500228849692985e-07, 2.2676682703375354e-07, 3.6761623155056537e-07],
          [2.2893255609858778e-21, 1.5441937013905384e-05, 4.074388072398815e-05,
           5.2161335424161561e-05, 8.1715492358327841e-05, 0.00010787745160225522],
          [-1.3823289456126841e-20, 1.4990965291904323e-06, 5.4295490891180001e-06,
           7.0167359356926422e-06, 1.0734035013065644e-05, 1.1996897236242909e-05],
          [-9.3212921887028401e-21, 1.0054658029823185e-05, 2.100366371930318e-05,
           2.9581071991349561e-05, 4.7096730648923597e-05, 6.1878442865411084e-05],
          [1.1610612220265047e-21, 5.3034152878328213e-08, 2.0647655253468472e-07,
           3.8813363311789536e-07, 4.2925372230907662e-07, 5.7786442341843429e-07],
          [4.204097029473209e-21, 3.2014460689876505e-07, 6.3473723229424336e-07,
           9.3531801076312435e-07, 1.2321077909034258e-06, 1.9490396100206598e-06],
          [4.9725555443470894e-21, 4.3150342819845607e-07, 1.2372421410133957e-06,
           1.5436778361815317e-06, 2.4031140467234829e-06, 3.2633808414170716e-06],
          [2.8658799060884637e-21, 1.2939868192529975e-06, 5.1235725945510507e-06,
           8.4038672840092365e-06, 1.0499785469438766e-05, 1.0802233241638491e-05]],
         [11, 29, 5, 25, 21, 8, 13, 22, 7, 31, 17, 30, 15, 24, 2, 9, 28,
          34, 35, 16, 3, 18, 12, 14])
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
    >>> print("The un-normalized linear FEM Laplace-Beltrami Spectrum is:\n")
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

    min_n_eigenvalues = 10 * n_eigenvalues

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

        # Determine the indices per label:
        label_indices = [i for i,x in enumerate(labels) if x == label]
        print('{0} vertices for label {1}'.
              format(len(label_indices), label))

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

        # Extract points and renumber faces for the selected indices:
        if len(select_indices) >= min_n_eigenvalues:
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
            print "The input size {0} for label {1} is too small. Skipped.".\
                format(len(select_indices), label)

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

    """
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.shapes.laplace_beltrami import fem_laplacian

    faces, lines, indices, points, npoints, scalars, name, input_vtk = read_vtk(vtk_file)

    try:
        print("{0}".format(fem_laplacian(points, faces, n_eigenvalues)))
    except RuntimeError:
        print "failed to compute LBS due to RuntimeError (e.g., singularity error)"


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
        fem_laplacian(points, faces, n_eigenvalues=5)))

    #print("The area-normalized linear FEM Laplace-Beltrami Spectrum is:\n\t{0}\n".format(
    #    fem_laplacian(points, faces, n_eigenvalues=3, normalization="area")))
