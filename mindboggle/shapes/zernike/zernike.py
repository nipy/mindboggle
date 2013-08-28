#!/usr/bin/python
"""
Compute the Zernike moments of a collection of points.


Authors:
    - Arthur Mikhno, 2013, Columbia University (original MATLAB code)
    - Brian Rossa, 2013, Tank Think Labs, LLC (port to Python)
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def zernike_moments(points, faces, order=10, scale_input=True,
                    decimate_fraction=0.75, decimate_smooth=25):
    """
    Compute the Zernike moments of a surface patch of points and faces.

    Optionally decimate the input mesh.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex
    faces : list of lists of 3 integers
        each list contains indices to vertices that form a triangle on a mesh
    order : integer
        order of the moments being calculated
    scale_input : Boolean
        translate and scale each object so it is bounded by a unit sphere?
        (this is the expected input to zernike_moments())
    decimate_fraction : float
        fraction of mesh faces to remove for decimation (1 for no decimation)
    decimate_smooth : integer
        number of smoothing steps for decimation

    Returns
    -------
    descriptors : list of floats
        Zernike descriptors

    Examples
    --------
    >>> # Example 1: simple cube (decimation doesn't have any effect):
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    >>> faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    >>> order = 3
    >>> decimate_fraction = 1
    >>> zernike_moments(points, faces, order, decimate_fraction)
    [0.0918881492369654,
     0.09357431096617608,
     0.04309029164656885,
     0.06466432586854755,
     0.03820155248327533,
     0.04138011726544602]
    >>> # Example 2: simple cube (with inner diagonal plane):
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> #path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = '/drop/cube.vtk'
    >>> faces, u1,u2, points, u3,u4,u5,u6 = read_vtk(vtk_file)
    >>> order = 3
    >>> decimate_fraction = 1
    >>> zernike_moments(points, faces, order, decimate_fraction)
    [0.0,
     1.5444366221695725e-21,
     0.0,
     2.081366518964347e-21,
     5.735003646768394e-05,
     2.433866250546253e-21]
    >>> # Example 3: Twins-2-1 left postcentral pial surface -- NO decimation:
    >>> # (zernike_moments took 142 seconds for order = 3 with no decimation)
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> faces, u1,u2, points, u3, labels, u4,u5 = read_vtk(label_file)
    >>> I22 = [i for i,x in enumerate(labels) if x==22] # postcentral
    >>> faces = remove_faces(faces, I22)
    >>> order = 3
    >>> scale_input = True
    >>> decimate_fraction = 1
    >>> zernike_moments(points, faces, order, scale_input, decimate_fraction)
     [0.005558794553842859,
     0.009838755429501177,
     0.003512500896236744,
     0.00899042745665395,
     0.001672289910738449,
     0.000919469614081582]
    >>> # Example 4: left postcentral -- with 0.75 decimation:
    >>> # Results closer to no decimation results than decimation with
    >>> # 0, 10, 20, 30, or 100 smoothing steps (21 seconds)
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> faces, u1,u2, points, u3, labels, u4,u5 = read_vtk(label_file)
    >>> I22 = [i for i,x in enumerate(labels) if x==22] # postcentral
    >>> faces = remove_faces(faces, I22)
    >>> order = 3
    >>> scale_input = True
    >>> decimate_fraction = 0.75
    >>> decimate_smooth = 25  # results closer to no decimation than 20 or 30
    >>> zernike_moments(points, faces, order, scale_input, decimate_fraction, decimate_smooth)
    [0.00528486237819844,
     0.009571754617699853,
     0.0033489494903015944,
     0.00875603468268444,
     0.0015879536633349918,
     0.0008080165707033097]
    >>> # Example 5: left postcentral + pars triangularis pial surfaces:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> faces, u1,u2, points, u3, labels, u4,u5 = read_vtk(label_file)
    >>> I20 = [i for i,x in enumerate(labels) if x==20] # pars triangularis
    >>> I22 = [i for i,x in enumerate(labels) if x==22] # postcentral
    >>> I22.extend(I20)
    >>> faces = remove_faces(faces, I22)
    >>> order = 3
    >>> scale_input = True
    >>> decimate_fraction = 0.75
    >>> decimate_smooth = 25
    >>> zernike_moments(points, faces, order, scale_input, decimate_fraction, decimate_smooth)
    [0.006591540793309832,
     0.010749937070447451,
     0.0034900573103799214,
     0.008972460360983864,
     0.0018220183025464518,
     0.0016893113500293917]

    """
    import numpy as np

    from mindboggle.utils.mesh import reindex_faces_0to1
    from mindboggle.shapes.zernike.multiproc import MultiprocPipeline
    from mindboggle.utils.mesh import decimate

    #-------------------------------------------------------------------------
    # Convert 0-indices (Python) to 1-indices (Matlab) for all face indices:
    #-------------------------------------------------------------------------
    index1 = True
    if index1:
        faces = reindex_faces_0to1(faces)

    # Convert lists to numpy arrays:
    if isinstance(points, list):
        points = np.array(points)
    if isinstance(faces, list):
        faces = np.array(faces)

    #-------------------------------------------------------------------------
    # Translate all points so that they are centered at their mean,
    # and scale them so that they are bounded by a unit sphere:
    #-------------------------------------------------------------------------
    if scale_input:
        center = np.mean(points, axis=0)
        points = points - center
        maxd = np.max(np.sqrt(np.sum(points**2, axis=1)))
        points = points / maxd

    #-------------------------------------------------------------------------
    # Decimate surface:
    #-------------------------------------------------------------------------
    if decimate_fraction < 1:
        points, faces, u1,u2 = decimate(points, faces,
            decimate_fraction, decimate_smooth, [], save_vtk=False)

        # Convert lists to numpy arrays:
        points = np.array(points)
        faces = np.array(faces)

    #-------------------------------------------------------------------------
    # Multiprocessor pipeline:
    #-------------------------------------------------------------------------
    pl = MultiprocPipeline()

    n_V = 3
    n_faces = pl.size(faces, 0)

    #-------------------------------------------------------------------------
    # Geometric moments:
    #-------------------------------------------------------------------------
    G = pl.geometric_moments_orig(points, faces, order, n_faces, n_V)

    #
    Z = pl.zernike(G, order)

    #-------------------------------------------------------------------------
    # Extract Zernike descriptors:
    #-------------------------------------------------------------------------
    descriptors = pl.feature_extraction(Z, order).tolist()

    return descriptors


def zernike_moments_per_label(vtk_file, order=10, exclude_labels=[-1],
                              scale_input=True,
                              decimate_fraction=0.75, decimate_smooth=25):
    """
    Compute the Zernike moments per labeled region in a file.

    Optionally decimate the input mesh.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file containing index scalars (labels)
    order : integer
        number of moments to compute
    exclude_labels : list of integers
        labels to be excluded
    scale_input : Boolean
        translate and scale each object so it is bounded by a unit sphere?
        (this is the expected input to zernike_moments())
    decimate_fraction : float
        fraction of mesh faces to remove for decimation (1 for no decimation)
    decimate_smooth : integer
        number of smoothing steps for decimation

    Returns
    -------
    descriptors_lists : list of lists of floats
        Zernike descriptors per label
    label_list : list of integers
        list of unique labels for which moments are computed

    Examples
    --------
    >>> # Uncomment "if label==22:" below to run example:
    >>> # Twins-2-1 left postcentral (22) pial surface:
    >>> import os
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments_per_label
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> order = 3
    >>> exclude_labels = [0]
    >>> scale_input = True
    >>> decimate_fraction = 0.75
    >>> decimate_smooth = 25
    >>> zernike_moments_per_label(vtk_file, order, exclude_labels, scale_input,
    >>>     decimate_fraction, decimate_smooth)
    ([[0.00528486237819844,
       0.009571754617699853,
       0.0033489494903015944,
       0.00875603468268444,
       0.0015879536633349918,
       0.0008080165707033097]],
     [22])

    """
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.utils.mesh import remove_faces
    from mindboggle.shapes.zernike.zernike import zernike_moments

    #-------------------------------------------------------------------------
    # Read VTK surface mesh file and area file:
    #-------------------------------------------------------------------------
    faces, u1,u2, points, u3, labels, u4,u5 = read_vtk(vtk_file)

    #-------------------------------------------------------------------------
    # Loop through labeled regions:
    #-------------------------------------------------------------------------
    ulabels = []
    [ulabels.append(int(x)) for x in labels if x not in ulabels
     if x not in exclude_labels]
    label_list = []
    descriptors_lists = []
    for label in ulabels:
      #if label==22:

        #---------------------------------------------------------------------
        # Determine the indices per label:
        #---------------------------------------------------------------------
        Ilabel = [i for i,x in enumerate(labels) if x == label]
        print('  {0} vertices for label {1}'.format(len(Ilabel), label))

        #---------------------------------------------------------------------
        # Remove background faces:
        #---------------------------------------------------------------------
        faces = remove_faces(faces, Ilabel)

        #---------------------------------------------------------------------
        # Compute Zernike moments for the label:
        #---------------------------------------------------------------------
        descriptors = zernike_moments(points, faces, order, scale_input,
                                      decimate_fraction, decimate_smooth)

        #---------------------------------------------------------------------
        # Append to a list of lists of spectra:
        #---------------------------------------------------------------------
        descriptors_lists.append(descriptors)
        label_list.append(label)

    return descriptors_lists, label_list
