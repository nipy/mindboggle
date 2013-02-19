#!/usr/bin/env python
"""
Functions related to reading and writing VTK format files.

1. Functions for reading VTK files
2. Functions for writing VTK files
3. Functions for converting FreeSurfer files to VTK format


Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
# Functions for reading VTK files
#=============================================================================
def read_vertices(Filename):
    """
    Load VERTICES segment from a VTK file (actually contains indices to vertices)

    Parameters
    ----------
    Filename : string
        The path/filename of a VTK format file.

    Returns
    -------
    indices : a list of integers
        Each element is an integer defined in VERTICES segment of the VTK file.
        The integer is an index referring to a point defined in POINTS segment of the VTK file.

    Notes ::

        We assume that VERTICES segment is organized as one line,
        the first column of which is the number of vertices.
        Vertices here are as vertices in VTK terminology. It may not be the vertices in your 3-D surface.

    """
    import vtk

    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(Filename)
    Reader.Update()

    Data = Reader.GetOutput()

    Vrts = Data.GetVerts()
    indices = [Vrts.GetData().GetValue(i) for i in xrange(1, Vrts.GetSize())]

    return indices

def read_lines(Filename):
    """
    Load LINES from a VTK file, along with the scalar values.

    The line that extracts vertices from a VTK
    iterates from 1 to Vrts.GetSize(), rather than from 0.

    Parameters
    ----------
    Filename : string
        The path/filename of a VTK format file.

    Returns
    -------
    lines : list of 2-tuple of integers
        Each element is a 2-tuple of IDs (i.e., indexes) of two points defined in POINTS segment of the VTK file
    scalars : list of floats
        Each element is a scalar value corresponding to a vertex

    """
    import vtk

    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(Filename)
    Reader.Update()

    Data = Reader.GetOutput()
    Lns = Data.GetLines()

    lines  = [[Lns.GetData().GetValue(j) for j in xrange(i*3+1, i*3+3) ]
              for i in xrange(Data.GetNumberOfLines())]

    PointData = Data.GetPointData()
    print "There are", Reader.GetNumberOfscalarsInFile(), "scalars in file", Filename
    print "Loading the scalar", Reader.GetScalarsNameInFile(0)
    ScalarsArray = PointData.GetArray(Reader.GetScalarsNameInFile(0))
    scalars = [ScalarsArray.GetValue(i) for i in xrange(0, ScalarsArray.GetSize())]

    return lines, scalars

def read_points(filename):
    """
    Load points of a VTK surface file.

    Parameters
    ----------
    filename : string
        path/filename of a VTK format file

    Returns
    -------
    points : list of lists of floats
        each element is a list of 3-D coordinates of a vertex on a surface mesh

    """
    import vtk

    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()

    Data = Reader.GetOutput()
    points = [list(Data.GetPoint(point_id))
              for point_id in xrange(Data.GetNumberOfPoints())]

    return points

def read_faces_points(filename):
    """
    Load points and faces of a VTK surface file.

    Parameters
    ----------
    filename : string
        path/filename of a VTK format file

    Returns
    -------
    faces : list of lists of integers
        each element is list of 3 indices of vertices that form a face
        on a surface mesh
    points : list of lists of floats
        each element is a list of 3-D coordinates of a vertex on a surface mesh
    npoints : integer
        number of points

    """
    import vtk

    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()

    Data = Reader.GetOutput()
    points = [list(Data.GetPoint(point_id))
              for point_id in xrange(Data.GetNumberOfPoints())]
    npoints = len(points)

    if Data.GetNumberOfPolys() > 0:
        faces = [[Data.GetPolys().GetData().GetValue(j)
                  for j in xrange(i*4 + 1, i*4 + 4)]
                  for i in xrange(Data.GetPolys().GetNumberOfCells())]
    else:
        faces = []

    return faces, points, npoints

def read_scalars(filename, return_first=True, return_array=False):
    """
    Load all scalar lookup tables from a VTK file.

    Parameters
    ----------
    filename : string
        The path/filename of a VTK format file.
    return_first : Boolean
        Return only the first list of scalar values?
    return_array : Boolean (only if return_first)
        Return first list of scalars as a numpy array?

    Returns
    -------
    scalars : list of lists of integers or floats
        each element is a list of scalar values for the vertices of a mesh
    scalar_name(s) : list of strings
        each element is the name of a lookup table

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> depths, name = read_scalars(depth_file)

    """
    import os
    import vtk
    if return_first and return_array:
        import numpy as np

    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()
    Data = Reader.GetOutput()
    PointData = Data.GetPointData()

    scalars = []
    scalar_names = []
    if Reader.GetNumberOfScalarsInFile() > 0:
        for scalar_index in range(Reader.GetNumberOfScalarsInFile()):
            scalar_name = Reader.GetScalarsNameInFile(scalar_index)

            n_scalars = scalar_index + 1
            if n_scalars == 1:
                print("Load \"{0}\" scalars from {1}".
                      format(scalar_name, os.path.basename(filename)))
            else:
                print("Load \"{0}\" (of {1} scalars) from {2}".
                      format(scalar_name, n_scalars, os.path.basename(filename)))

            scalar_array = PointData.GetArray(scalar_name)
            if scalar_array:
                scalar = [scalar_array.GetValue(i)
                          for i in xrange(scalar_array.GetSize())]
            else:
                print "Empty scalar -- Please check the source VTK"
                exit()
            scalars.append(scalar)
            scalar_names.append(scalar_name)

    if return_first:
        if len(scalars):
            scalars = scalars[0]
        if return_array:
            scalars = np.array(scalars)
        if len(scalar_names):
            scalar_names = scalar_names[0]
        else:
            scalar_names = ''

    return scalars, scalar_names


def read_vtk(filename, return_first=True, return_array=False):
    """
    Load faces, lines, indices, points, #points,
    and all scalar lookup tables from a VTK file.

    Note ::

        1. This supports copying lines, vertices (indices of points),
           and triangular faces from one surface to another.
        2. We assume that all vertices are written in one line in VERTICES segment.

    Parameters
    ----------
    filename : string
        The path/filename of a VTK format file.
    return_first : Boolean
        Return only the first list of scalar values?
    return_array : Boolean (only if return_first)
        Return first list of scalars as a numpy array?

    Returns
    -------
    faces : list of lists of integers
        each element is list of 3 indices of vertices that form a face
        on a surface mesh
    lines : list of 2-tuples of integers
        each element is an edge on the mesh, consisting of 2 integers
        representing the 2 vertices of the edge
    indices : list of integers
        indices of vertices
    points :  list of 3-tuples of floats
        each element has 3 numbers representing the coordinates of the points
    npoints : int
        number of vertices in the mesh
    scalars : array or list or list of lists of integers or floats
        array or list(s) of scalar values for the vertices of a mesh
    scalar_names : string or list of strings
        name(s) of lookup table(s)

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file)

    """
    import os
    import vtk
    if return_first and return_array:
        import numpy as np

    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()

    Data = Reader.GetOutput()
    PointData = Data.GetPointData()
    points = [list(Data.GetPoint(point_id))
              for point_id in xrange(0, Data.GetNumberOfPoints())]
    npoints = len(points)

    if Data.GetNumberOfPolys() > 0:
        faces = [[Data.GetPolys().GetData().GetValue(j)
                  for j in xrange(i*4 + 1, i*4 + 4)]
                  for i in xrange(Data.GetPolys().GetNumberOfCells())]
    else:
        faces = []

    if Data.GetNumberOfLines() > 0:
        lines  = [[Data.GetLines().GetData().GetValue(j)
                   for j in xrange(i*3+1, i*3+3) ]
                   for i in xrange(Data.GetNumberOfLines())]
    else:
        lines = []

    if Data.GetNumberOfVerts() > 0:
        indices = [Data.GetVerts().GetData().GetValue(i)
                    for i in xrange(1, Data.GetVerts().GetSize() )]
        # The reason the reading starts from 1 is because we need to avoid the
    else:
        indices = []

    scalars = []
    scalar_names = []

    if Reader.GetNumberOfScalarsInFile() > 0:
        for scalar_index in range(Reader.GetNumberOfScalarsInFile()):
            scalar_name = Reader.GetScalarsNameInFile(scalar_index)

            n_scalars = scalar_index + 1
            if n_scalars == 1:
                print("Load \"{0}\" scalars from {1}".
                      format(scalar_name, os.path.basename(filename)))
            else:
                print("Load \"{0}\" (of {1} scalars) from {2}".
                      format(scalar_name, n_scalars, os.path.basename(filename)))

            scalar_array = PointData.GetArray(scalar_name)
            if scalar_array:
                scalar = [scalar_array.GetValue(i)
                          for i in xrange(scalar_array.GetSize())]
                scalars.append(scalar)
                scalar_names.append(scalar_name)

    if return_first:
        if len(scalars):
            scalars = scalars[0]
        if return_array:
            scalars = np.array(scalars)
        if len(scalar_names):
            scalar_names = scalar_names[0]
        else:
            scalar_names = ''

    return faces, lines, indices, points, npoints, scalars, scalar_names

#=============================================================================
# Functions for writing VTK elements/files
#=============================================================================

def write_header(Fp, Header='# vtk DataFile Version 2.0',
                     Title='Generated by Mindboggle (www.mindboggle.info)',
                     fileType='ASCII', dataType='POLYDATA'):
    """
    Write header information for a VTK-format file::

        # vtk DataFile Version 2.0
        Generated by Mindboggle (www.mindboggle.info)
        ASCII
        DATASET POLYDATA

    This part matches three things in the VTK 4.2 File Formats doc:
      - Part 1: Header
      - Part 2: Title (256 characters maximum, terminated with newline character)
      - Part 3: Data type, either ASCII or BINARY
      - Part 4: Geometry/topology. dataType is one of:
          - STRUCTURED_POINTS
          - STRUCTURED_GRID
          - UNSTRUCTURED_GRID
          - POLYDATA
          - RECTILINEAR_GRID
          - FIELD

    """

    Fp.write('{0}\n{1}\n{2}\nDATASET {3}\n'.format(Header, Title, fileType, dataType))

def write_points(Fp, points, dataType="float"):
    """
    Write coordinates of points, the POINTS section in DATASET POLYDATA section::

        POINTS 150991 float
        -7.62268877029 -81.2403945923 -1.44539153576
        ...

    Indices are 0-offset. Thus the first point is point id 0::

        POINTS n dataType
        p0x p0y p0z
        ...
        p(n-1)x p(n-1)y p(n-1)z

    """
    import numpy as np

    Fp.write('POINTS {0} {1}\n'.format(len(points), dataType))

    n = np.shape(points)[1]
    for point in points:
        if n == 3:
            [R, A, S] = point
            Fp.write('{0} {1} {2}\n'.format(R, A, S))
        elif n == 2:
            [R, A] = point
            Fp.write('{0} {1}\n'.format(R, A))
        else:
            print('ERROR: Unrecognized number of coordinates per point')

def write_faces(Fp, faces):
    """
    Write indices to vertices forming triangular meshes or lines,
    the POLYGONS section in DATASET POLYDATA section:

        POLYGONS 301978 1207912
        3 0 1 4
        ...

    """
    import numpy as np

    n = np.shape(faces)[1]
    if n == 3:
        face_name = 'POLYGONS '
        Fp.write('{0} {1} {2}\n'.format(face_name, len(faces),
                 len(faces) * (n + 1)))
    elif n == 2:
        face_name = 'LINES '
        Fp.write('{0} {1} {2}\n'.format(face_name, len(faces),
                 len(faces) * (n + 1)))
    else:
        print('ERROR: Unrecognized number of vertices per face')

    for face in faces:
        if n == 3:
            [V0, V1, V2] = face
            Fp.write('{0} {1} {2} {3}\n'.format(n, V0, V1, V2))
        elif n == 2:
            [V0, V1] = face
            Fp.write('{0} {1} {2}\n'.format(n, V0, V1))

def write_lines(Fp, lines):
    """
    Save connected line segments to a VTK file.

    Parameters
    ----------
    Fp: pointer to a file
        pointer to the file to write lines
    lines : list of 2-tuples of integers
        each element is an edge on the mesh, consisting of 2 integers
        representing the 2 vertices of the edge
    """

    write_faces(Fp, lines)

def write_vertices(Fp, indices):
    """
    Write indices to vertices, the VERTICES section
    in the DATASET POLYDATA section::

        VERTICES 140200 420600
        3 130239 2779 10523
        ...

    Indices are 0-offset. Thus the first point is point id 0::

        VERTICES n size
        numPoints0 i0 j0 k0
        ...
        numPoints_[n-1] i_[n-1] j_[n-1] k_[n-1]

    Notes
    -------

    Currently we write all vertices in one line.

    """

    Fp.write('VERTICES {0} {1}\n{2} '.format(
             1, len(indices) + 1, len(indices)))
    [Fp.write('{0} '.format(i)) for i in indices]
    Fp.write('\n')

def write_scalars(Fp, scalars, scalar_name, begin_scalars=True):
    """
    Write per-VERTEX values as a scalar lookup table into a VTK file::

        POINT_DATA 150991
        SCALARS Max_(majority_labels) int 1
        LOOKUP_TABLE default
        11 11 11 11 11 11 11 11 11 11 ...

    Parameters
    ----------
    scalars :  list of floats
    begin_scalars : [Boolean] True if the first vertex lookup table in a VTK file

    """

    if begin_scalars:
        Fp.write('POINT_DATA {0}\n'.format(len(scalars)))
    Fp.write('SCALARS {0} float\n'.format(scalar_name))
    Fp.write('LOOKUP_TABLE {0}\n'.format(scalar_name))
    for Value in scalars:
        Fp.write('{0}\n'.format(Value))
    Fp.write('\n')


def write_vtk(output_vtk, points, indices=[], lines=[], faces=[],
              scalars=[], scalar_names=['scalars']):
    """
    Save lists of scalars into the lookup table of a VTK-format file.

    Scalar definition includes specification of a lookup table.
    The definition of a lookup table is optional. If not specified,
    the default VTK table will be used (and tableName should be "default").

    SCALARS dataName dataType numComp
    LOOKUP_TABLE tableName

    Parameters
    ----------
    output_vtk : string
        path of the output VTK file
    points :  list of 3-tuples of floats
        each element has 3 numbers representing the coordinates of the points
    indices : list of integers
        indices of vertices, default=[]
    lines : list of 2-tuples of integers
        Each element is an edge on the mesh, consisting of 2 integers
        representing the 2 vertices of the edge
        default=[]
    faces : list of 3-tuples of integers
        indices to the three vertices of a face on the mesh, default=[]
    scalars : list of lists of floats (or single list or array of floats)
        each list (lookup table) contains values assigned to the vertices, default=[]
    scalar_names : string or list of strings
        each element is the name of a lookup table, default=['scalars']
        if only one string is given for this field, the program will convert
        it into a list of only this string.

    Notes
    --------

    If you do not have all 7 parameters, it's safer to use syntax like
    ...``indices=indices, faces=faces``... (as in Toy example)
    than syntax like ...``indices, faces``... alone to
    ensure that your variables are aligned with their orders
    in parameter definition.

    Examples
    --------
    >>> # Toy example
    >>> import random, os
    >>> from mindboggle.utils.io_vtk import write_vtk
    >>> points = [[random.random() for i in [1,2,3]] for j in xrange(4)]
    >>> indices = [1,2,3,0]
    >>> lines = [[1,2],[3,4]]
    >>> faces = [[1,2,3],[0,1,3]]
    >>> scalar_names = ['curv','depth']
    >>> scalars = [[random.random() for i in xrange(4)] for j in [1,2]]
    >>> write_vtk('test_write_vtk.vtk', points,
    >>>          indices, lines, faces, scalars, scalar_names)
    >>> os.system('mayavi2 -m Surface -d test_write_vtk.vtk &')
    >>>
    >>> # Write vtk file with depth values on sulci
    >>> import os
    >>> from mindboggle.utils.mesh_operations import inside_faces
    >>> from mindboggle.utils.io_vtk import read_vtk, write_vtk
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file)
    >>> # Write to vtk file and view with mayavi2:
    >>> write_vtk('test_write_vtk.vtk', points, [], [], faces, depths, 'depths')
    >>> os.system('mayavi2 -m Surface -d test_write_vtk.vtk &')

    """
    import os
    from mindboggle.utils.io_vtk import write_header, write_points, \
         write_vertices, write_faces, write_scalars, \
         scalars_checker

    output_vtk = os.path.join(os.getcwd(), output_vtk)

    Fp = open(output_vtk,'w')
    write_header(Fp)
    write_points(Fp, points)

    if len(indices):
        write_vertices(Fp, indices)
    if len(lines):
        for i in xrange(0,len(lines)):
            lines[i] = [lines[i][0], lines[i][1]]
        write_faces(Fp, lines) # write_faces can write either lines or faces
    if len(faces):
        write_faces(Fp, faces)
    if len(scalars) > 0:
        scalars, scalar_names = scalars_checker(scalars, scalar_names)

        for i, scalar_list in enumerate(scalars):
            if i == 0:
                scalar_name = scalar_names[i]
                write_scalars(Fp, scalar_list, scalar_name)
            else:
                if len(scalar_names) < i + 1:
                    scalar_name = scalar_names[0]
                else:
                    scalar_name  = scalar_names[i]
                write_scalars(Fp, scalar_list, scalar_name, begin_scalars=False)
    else:
        print('Error: scalars is empty')
        exit()

    Fp.close()

    return output_vtk

def rewrite_scalars(input_vtk, output_vtk, new_scalars,
                    new_scalar_names=['scalars'], filter_scalars=[]):
    """
    Load VTK format file and save a subset of scalars into a new file.

    Parameters
    ----------
    input_vtk : string
        input VTK file name
    output_vtk : string
        output VTK file name
    new_scalars : list of lists of floats (or single list or array of floats)
        each list (lookup table) contains new values to assign to the vertices
    new_scalar_names : string or list of strings
        each element is the new name for a lookup table
    filter_scalars : list or numpy array (optional)
        scalar values used to filter faces (values > -1 retained)

    Returns
    -------
    output_vtk : string
        output VTK file name

    Examples
    --------
    >>> # Write vtk file with depth values on sulci
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> depths, name = read_scalars(depth_file)
    >>> sulci_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'features', 'lh.sulci.vtk')
    >>> sulci, name = read_scalars(sulci_file)
    >>> # Write to vtk file and view with mayavi2:
    >>> rewrite_scalars(depth_file, 'test_rewrite_scalars.vtk',
    >>>                      depths, 'depths', sulci)
    >>> os.system('mayavi2 -m Surface -d test_rewrite_scalars.vtk &')

    """
    import os
    from mindboggle.utils.mesh_operations import inside_faces
    from mindboggle.utils.io_vtk import write_header, write_points, \
         write_vertices, write_faces, write_scalars, read_vtk, \
         scalars_checker

    # Output VTK file to current working directory
    output_vtk = os.path.join(os.getcwd(), output_vtk)

    # Load VTK file
    faces, lines, indices, points, npoints, scalars, name = read_vtk(input_vtk)

    # Find indices to nonzero values
    indices = range(npoints)
    if len(filter_scalars):
        indices_filter = [i for i,x in enumerate(filter_scalars) if x > -1]
        indices_remove = [i for i,x in enumerate(filter_scalars) if x == -1]
        # Remove surface mesh faces whose three vertices are not all in indices
        faces = inside_faces(faces, indices_filter)

    # Write VTK file
    Fp = open(output_vtk,'w')
    write_header(Fp)
    write_points(Fp, points)
    if len(indices):
        write_vertices(Fp, indices)
    if len(faces):
        write_faces(Fp, faces)
    if len(new_scalars) > 0:
        new_scalars, new_scalar_names = scalars_checker(new_scalars, new_scalar_names)

        for i, new_scalar_list in enumerate(new_scalars):
            for iremove in indices_remove:
                new_scalar_list[iremove] = -1
            if i == 0:
                new_scalar_name = new_scalar_names[0]
                write_scalars(Fp, new_scalar_list, new_scalar_name)
            else:
                if len(new_scalar_names) < i + 1:
                    new_scalar_name = new_scalar_names[0]
                else:
                    new_scalar_name  = new_scalar_names[i]
                write_scalars(Fp, new_scalar_list, new_scalar_name, begin_scalars=False)
    else:
        print('Error: new_scalars is empty')
        exit()

    Fp.close()

    return output_vtk

def copy_scalars(output_vtk, points, faces, lines, indices, scalars,
                 scalar_names=['scalars'], header="Created by MindBoggle"):
    """
    Copy the scalars of one surface to another surface (same number of vertices).

    Scalar definition includes specification of a lookup table.
    The definition of a lookup table is optional. If not specified,
    the default VTK table will be used (and tableName should be "default").

    SCALARS dataName dataType numComp
    LOOKUP_TABLE tableName

    Parameters
    ----------
    output_vtk : string
        path of the output VTK file
    points :  list of 3-tuples of floats
        each element has 3 numbers representing the coordinates of the points
    indices : list of integers
        indices of vertices. they are called VERTICES in VTK format though
    faces : list of 3-tuples of integers
        indices to the three vertices of a face on the mesh
    lines : list of 2-tuples of integers
        Each element contains two indices  of the two mesh vertexes forming a line on the mesh
    scalars : list of lists of floats
        each list contains values assigned to the vertices
    scalar_names : list of strings
        each element is the name of a lookup table of scalars values
    header : string
        title for the target VTK file

    Returns
    -------
    output_vtk : string
        output VTK file name

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.mesh_operations import inside_faces
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> file1 = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                 'measures', 'lh.pial.depth.vtk')
    >>> # Don't have this yet:
    >>> file2 = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                 'surf', 'lh.inflated.vtk')
    >>> file3 = 'test_copy_scalars.vtk'
    >>> faces, lines, indices, points, npoints, scalars, \
    >>>     scalar_names = read_vtk(file1, return_first=False, return_array=False)
    >>> output_vtk = copy_scalars(file3, points, faces, lines, indices,
    >>>                           scalars, scalar_names, header)
    >>> # View with mayavi2:
    >>> os.system('mayavi2 -m Surface -d test_copy_scalars.vtk &')

    """
    import os
    from mindboggle.utils.io_vtk import write_header, write_points,\
        write_faces, write_lines, write_vertices, write_scalars,\
        scalars_checker

    output_vtk = os.path.join(os.getcwd(), output_vtk)

    Fp = open(output_vtk,'w')
    write_header(Fp, header)
    write_points(Fp, points)
    if len(faces):
        write_faces(Fp, faces)
    if len(lines):
        write_lines(Fp, lines)
    if len(indices):
        write_vertices(Fp, indices)
    if len(scalars) > 0:

#        # Make sure that new_scalars is a list
#        if type(scalars) != list:
#            if type(scalars[0]) == np.ndarray:
#                scalars = scalars.tolist()
#            scalars = [scalars]
#        # Make sure that new_scalars is a list of lists
#        if type(scalars[0]) != list:
#            scalars = [scalars]
#        # Make sure that scalar_names is a list
#        if type(scalar_names) != list:
#            scalar_names = [scalar_names]

        scalars, scalar_names = scalars_checker(scalars, scalar_names)

        for i, scalar_list in enumerate(scalars):
            if i == 0:
                scalar_name = scalar_names[i]
                write_scalars(Fp, scalar_list, scalar_name)
            else:
                if len(scalar_names) < i + 1:
                    scalar_name = scalar_names[0]
                else:
                    scalar_name  = scalar_names[i]
                write_scalars(Fp, scalar_list, scalar_name, begin_scalars=False)
    else:
        print('Error: scalars is empty')
        exit()

    Fp.close()

    return output_vtk

def explode_scalars(input_vtk, output_stem, exclude_values=[-1],
                    background_value=-1, output_scalar_name='scalars'):
    """
    Write out a separate VTK file for each integer (>-1)
    in (the first) scalar list of an input VTK file.

    Parameters
    ----------
    input_vtk : string
        path of the input VTK file
    output_stem : string
        path and stem of the output VTK file
    exclude_values : list or array
        values to exclude
    background_value : integer or float
        background value in output VTK files
    scalar_name : string
        name of a lookup table of scalars values

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import explode_scalars
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> sulci_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'features', 'lh.sulci.vtk')
    >>> output_stem = 'sulcus'
    >>> explode_scalars(sulci_file, output_stem)
    >>> example_vtk = os.path.join(os.getcwd(), output_stem + '0.vtk')
    >>> os.system('mayavi2 -m Surface -d ' + example_vtk + ' &')

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, write_vtk

    # Load VTK file
    faces, lines, indices, points, npoints, scalars, \
        scalar_names = read_vtk(input_vtk, return_first=False, return_array=False)
    print("Explode the scalar list in {0}".format(os.path.basename(input_vtk)))

    # Use first scalar list as a numpy array
    scalars = np.array(scalars[0])

    # Loop through unique (non-excluded) scalar values
    unique_scalars = [int(x) for x in np.unique(scalars) if x not in exclude_values]
    for scalar in unique_scalars:

        # Create array and indices for scalar value
        new_scalars = np.copy(scalars)
        new_scalars[scalars != scalar] = background_value
        indices = [i for i,x in enumerate(new_scalars) if x == scalar]
        print("  Scalar {0}: {1} vertices".format(scalar, len(indices)))

        # Write VTK file with scalar value
        output_vtk = os.path.join(os.getcwd(), output_stem + str(scalar) + '.vtk')
        write_vtk(output_vtk, points, indices, lines, faces,
                  [new_scalars.tolist()], [output_scalar_name])

def scalars_checker(scalars, scalar_names):
    """
    Check whether input scalars and scalar_names are in acceptable format. If not, reformat.

    Parameters
    ----------
    scalars : list of lists of floats (or single list or 1-/2-D array of floats)
    scalar_names : string or list of strings

    Examples
    --------
    >>> from mindboggle.utils.io_vtk import scalars_checker
    >>> scalars_checker([[1,2],[3,4]], "")
    >>> ([[1, 2], [3, 4]], [''])
    >>> scalars_checker([1,2,3,4], ["123"])
    >>> ([[1, 2, 3, 4]], ['123'])
    >>> scalars_checker(1, ["123"])
         Error: scalars is neither a list nor a numpy array.
    >>> % You will be kicked out from Python shell after the command above
    >>> import numpy as np
    >>> scalars_checker(np.array([1,2,3]), ["123"])
         Warning: the new_scalars is a 1-D numpy array. Conversion done but may have problems in final VTK.
    >>> ([[1, 2, 3]], ['123'])
    >>> scalars_checker(np.array([[1,2,3]]), ["123"])
    >>> ([[1, 2, 3]], ['123'])
    >>> scalars_checker(np.array([[1,2,3],[4,5,6]]), ["123"])
    >>> ([[1, 2, 3], [4, 5, 6]], ['123'])
    >>> scalars_checker(np.array([[[1,2,3]]]), ["123"])
        Error: Dimension of new_scalars is too high.

    Notes
    -----
    This function does not check all possible cases of scalars and scalar_names,
    but only those that are likely to happen when using Mindboggle.

    """
    import numpy as np
    if type(scalars) != list:
        if type(scalars) == np.ndarray:
            if len(scalars.shape) < 2: # this is a 1-D array
                scalars = np.array([scalars]) # increase dimension by 1
                print "Warning: scalars is a 1-D numpy array. Conversion done but may have problems in final VTK. "
                scalars = scalars.tolist()
            elif len(scalars.shape) == 2: # 2-D numpy array
                scalars = scalars.tolist()
            else:
                print "Error: Dimension of new_scalars is too high."
                exit()
        else:
            print "Error: scalars is neither a list nor a numpy array. "
            exit()
    else:  # a list, but may be 1-D
        if type(scalars[0]) == int or type(scalars[0]) == float: # this is an acceptable 1-D list
            scalars = [scalars]
        elif type(scalars[0]) == list:
            pass
        else:
            print "io_vtk.py: Error: scalars is a 1-D list containing unacceptable elements. "
            print "io_vtk.py: scalars type is:", type(scalars)
            print "io_vtk,py: scalar length is:", len(scalars)
            print "io_vtk.py: scalars[0] type is:", type(scalars[0])
            exit()

    if type(scalar_names) == str:
        scalar_names = [scalar_names]
    elif  type(scalar_names) == list:
        pass
    else:
        print "Error: scalar_names is neither a list nor a string"
        exit()

    return scalars, scalar_names

#=============================================================================
# Functions for writing tables
#=============================================================================

def write_mean_shapes_table(table_file, column_names, labels, depth_file,
                            mean_curvature_file, gauss_curvature_file,
                            max_curvature_file, min_curvature_file,
                            thickness_file='', convexity_file='',
                            norm_vtk_file='', exclude_labels=[]):
    """
    Make a table of mean values per label per measure.

    Parameters
    ----------
    filename :  output filename (without path)
    column_names :  names of columns [list of strings]
    labels :  name of label file or list of labels (same length as values)
    *shape_files :  arbitrary number of vtk files with scalar values
    norm_vtk_file :  name of file containing per-vertex normalization values
                     (e.g., surface areas)
    exclude_labels : list of integer labels to be excluded

    Returns
    -------
    means_file : table file name for mean shape values
    norm_means_file : table file name for mean shape values normalized by area

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import write_mean_shapes_table
    >>> table_file = 'test_write_mean_shapes_table.txt'
    >>> column_names = ['labels', 'area', 'depth', 'mean_curvature',
    >>>                 'gauss_curvature', 'max_curvature', 'min_curvature']
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> subject = 'MMRR-21-1'
    >>> labels_file = os.path.join(data_path, 'subjects', subject,
    >>>                            'labels', 'lh.labels.DKT25.manual.vtk')
    >>> exclude_values = [-1]
    >>> area_file = os.path.join(data_path, 'subjects', subject,
    >>>                                     'measures', 'lh.pial.area.vtk')
    >>> depth_file = os.path.join(data_path, 'subjects', subject,
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> mean_curvature_file = os.path.join(data_path, 'subjects', subject,
    >>>                                    'measures', 'lh.pial.curv.avg.vtk')
    >>> gauss_curvature_file = os.path.join(data_path, 'subjects', subject,
    >>>                                     'measures', 'lh.pial.curv.gauss.vtk')
    >>> max_curvature_file = os.path.join(data_path, 'subjects', subject,
    >>>                                   'measures', 'lh.pial.curv.max.vtk')
    >>> min_curvature_vector_file = os.path.join(data_path, 'subjects',
    >>>     subject, 'measures', 'lh.pial.curv.min.dir.txt')
    >>> write_mean_shapes_table(table_file, column_names, labels_file,
    >>>                         depth_file, mean_curvature_file,
    >>>                         gauss_curvature_file,
    >>>                         max_curvature_file, min_curvature_file,
    >>>                         thickness_file='', convexity_file='',
    >>>                         area_file, exclude_values)

    """
    import os
    import numpy as np
    from mindboggle.measure.measure_functions import mean_value_per_label
    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.io_file import write_table

    # Load per-vertex labels and normalization vtk files
    if type(labels) == str:
        labels, name = read_scalars(labels, return_first=True, return_array=True)
    if len(norm_vtk_file):
        norms, name = read_scalars(norm_vtk_file)
    else:
        norms = [1 for x in labels]

    # List files
    vtk_files = [depth_file, mean_curvature_file, gauss_curvature_file,
                 max_curvature_file, min_curvature_file,
                 thickness_file, convexity_file]

    # Append columns of values to table
    columns = []
    norm_columns = []
    for i, vtk_file in enumerate(vtk_files):
        if len(vtk_file):
            values, name = read_scalars(vtk_file)
            mean_values, norm_mean_values, norm_values, \
                label_list = mean_value_per_label(values, norms, labels,
                                                  exclude_labels)

            columns.append(mean_values)
            norm_columns.append(norm_mean_values)
        else:
            del(column_names[i])

    # Prepend with column of normalization values
    columns.insert(0, norm_values)
    norm_columns.insert(0, norm_values)
    column_names.insert(0, 'area')

    # Prepend with column of labels and write tables
    column_names.insert(0, 'label')

    means_file = os.path.join(os.getcwd(), table_file)
    write_table(label_list, columns, column_names, means_file)

    norm_means_file = os.path.join(os.getcwd(), 'norm_' + table_file)
    write_table(label_list, norm_columns, column_names, norm_means_file)

    return means_file, norm_means_file

def write_vertex_shapes_table(table_file, column_names,
                              labels_file, sulci_file, fundi_file,
                              area_file, depth_file,
                              mean_curvature_file, gauss_curvature_file,
                              max_curvature_file, min_curvature_file,
                              thickness_file='', convexity_file=''):
    """
    Make a table of shape values per vertex.

    Parameters
    ----------
    table_file : output filename (without path)
    column_names : names of columns [list of strings]
    *vtk_files : arbitrary number of vtk files with per-vertex scalar values
                 (set each missing file to an empty string)

    Returns
    -------
    shape_table : table file name for vertex shape values

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import write_vertex_shape_table
    >>> from mindboggle.utils.mesh_operations import find_neighbors, detect_boundaries
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> table_file = 'test_write_vertex_shape_table.txt'
    >>> column_names = ['sulci', 'area', 'depth', 'mean_curvature',
    >>>                 'gauss_curvature', 'max_curvature', 'min_curvature']
    >>> subject = 'MMRR-21-1'
    >>> labels_file = ''
    >>> fundi_file = ''
    >>> sulci_file = os.path.join(data_path, 'subjects', subject,
    >>>                                      'features', 'lh.sulci.vtk')
    >>> area_file = os.path.join(data_path, 'subjects', subject,
    >>>                                     'measures', 'lh.pial.area.vtk')
    >>> depth_file = os.path.join(data_path, 'subjects', subject,
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> mean_curvature_file = os.path.join(data_path, 'subjects', subject,
    >>>                                    'measures', 'lh.pial.curv.avg.vtk')
    >>> gauss_curvature_file = os.path.join(data_path, 'subjects', subject,
    >>>                                     'measures', 'lh.pial.curv.gauss.vtk')
    >>> max_curvature_file = os.path.join(data_path, 'subjects', subject,
    >>>                                   'measures', 'lh.pial.curv.max.vtk')
    >>> min_curvature_vector_file = os.path.join(data_path, 'subjects',
    >>>     subject, 'measures', 'lh.pial.curv.min.dir.txt')
    >>> write_vertex_shape_table(table_file, column_names,
    >>>     labels_file, sulci_file, fundi_file, area_file, depth_file,
    >>>     mean_curvature_file, gauss_curvature_file,
    >>>     max_curvature_file, min_curvature_file, '', '')

    """
    import os
    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.io_file import write_table

    # List files
    vtk_files = [labels_file, sulci_file, fundi_file, area_file, depth_file,
                 mean_curvature_file, gauss_curvature_file,
                 max_curvature_file, min_curvature_file,
                 thickness_file, convexity_file]

    # Append columns of values to table
    columns = []
    for i, vtk_file in enumerate(vtk_files):
        if len(vtk_file):
            values, name = read_scalars(vtk_file)
            if len(columns) == 0:
                indices = range(len(values))
            columns.append(values)
        else:
            del(column_names[i])

    # Prepend with column of indices and write table
    column_names.insert(0, 'index')
    shape_table = os.path.join(os.getcwd(), table_file)
    write_table(indices, columns, column_names, shape_table)

    return shape_table

#=============================================================================
# Functions for converting FreeSurfer files to VTK format
#=============================================================================

def freesurface_to_vtk(surface_file):
    """
    Convert FreeSurfer surface file to VTK format.

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import freesurface_to_vtk
    >>> subjects_path = os.environ['SUBJECTS_DIR']
    >>> surface_file = os.path.join(subjects_path, 'MMRR-21-2', 'surf', 'lh.pial')
    >>> freesurface_to_vtk(surface_file)

    """
    import os
    from mindboggle.utils.io_free import read_surface
    from mindboggle.utils.io_vtk import write_header, write_points, write_faces

    points, faces = read_surface(surface_file)

    output_vtk = os.path.join(os.getcwd(),
                            os.path.basename(surface_file + '.vtk'))
    Fp = open(output_vtk, 'w')
    write_header(Fp, Title='vtk output from ' + surface_file)
    write_points(Fp, points)
    write_faces(Fp, faces)
    Fp.close()

    return output_vtk

def freecurvature_to_vtk(file_string, surface_file, hemi, subject, subjects_path):
    """
    Convert FreeSurfer curvature, thickness, or convexity file to VTK format.

    Parameters
    ----------
    file_string : string
        string for FreeSurfer file: 'curv', 'thickness', 'sulc.pial'
    surface_file : string  (name of VTK surface file)
    hemi : string indicating left or right hemisphere
    subject : string
        name of subject directory
    subjects_path: string
        path to subject directory

    Returns
    -------
    output_vtk : string
        name of output VTK file, where each vertex is assigned
        the corresponding shape value.

    """
    import os
    from mindboggle.utils.io_free import read_curvature
    from mindboggle.utils.io_vtk import read_vtk, write_vtk

    filename = os.path.join(subjects_path, subject, 'surf',
                            hemi + '.' + file_string)
    output_vtk = os.path.join(os.getcwd(), file_string + '.vtk')

    curvature_values = read_curvature(filename)

    # Load VTK surface
    faces, lines, indices, points, npoints, scalars, name = read_vtk(surface_file)

    scalars = [curvature_values]
    scalar_names = [file_string]
    write_vtk(output_vtk, points, indices, lines, faces, scalars, scalar_names)

    return output_vtk

def freeannot_to_vtk(surface_file, hemi, subject, subjects_path, annot_name):
    """
    Load a FreeSurfer .annot file and save as a VTK format file.

    Parameters
    ----------
    surface_file : string  (name of VTK surface file)
    annot_file : strings  (name of FreeSurfer .annot file)

    Returns
    -------
    labels : list of integers (one label per vertex)
    output_vtk : output VTK file

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import freeannot_to_vtk
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                        'measures', 'lh.pial.depth.vtk')
    >>> subjects_path = os.environ['SUBJECTS_DIR']
    >>> annot_name = "aparcNMMjt" #labels.DKT31.manual"
    >>> freeannot_to_vtk(surface_file, "lh", MMRR-21-1, subjects_path, annot_name)

    """
    import os
    import nibabel as nb
    from mindboggle.utils.io_vtk import read_vtk, write_vtk

    annot_file = os.path.join(subjects_path, subject, 'label',
                              hemi + '.' + annot_name + '.annot')

    labels, colortable, names = nb.freesurfer.read_annot(annot_file)

    # Load FreeSurfer surface
    #from utils.io_file import read_surface
    #points, faces = read_surface(surface_file)

    # Load VTK surface
    faces, lines, indices, points, npoints, scalars, name = read_vtk(surface_file)

    output_stem = os.path.join(os.getcwd(),
                  os.path.basename(surface_file.strip('.vtk')))
    output_vtk = output_stem + '.' + annot_name.strip('.annot') + '.vtk'

    scalars = [labels.tolist()]
    scalar_names = ['Labels']
    write_vtk(output_vtk, points, indices, lines, faces, scalars, scalar_names)

    return labels, output_vtk

def vtk_to_freelabels(hemi, surface_file, label_numbers, label_names,
                      RGBs, scalar_name):
    """
    Write FreeSurfer .label files from a labeled VTK surface mesh.

    From https://surfer.nmr.mgh.harvard.edu/fswiki/LabelsClutsAnnotationFiles:

        "A label file is a text file capturing a list of vertices belonging to a region,
        including their spatial positions(using R,A,S coordinates). A label file
        corresponds only to a single label, thus contains only a single list of vertices"::

            1806
            7  -22.796  -66.405  -29.582 0.000000
            89  -22.273  -43.118  -24.069 0.000000
            138  -14.142  -81.495  -30.903 0.000000
            [...]

    Parameters
    ----------
    hemi :  hemisphere [string]
    surface_file :  vtk surface mesh file with labels [string]
    label_numbers :  label numbers [list of strings]
    label_names :  label names [list of strings]
    RGBs :  list of label RGB values for later conversion to a .annot file
    scalar_name :  name of scalar values in vtk file [string]

    Returns
    -------
    label_files :  list of .label file names (order must match label list)
    colortable :  file with list of labels and RGB values
                 NOTE: labels are identified by the colortable's RGB values

    """
    import os
    import numpy as np
    import vtk

    def string_vs_list_check(var):
        """
        Check type to make sure it is a string.

        (if a list, return the first element)
        """

        # Check type:
        if type(var) == str:
            return var
        elif type(var) == list:
            return var[0]
        else:
            os.error("Check format of " + var)

    # Check type to make sure the filename is a string
    # (if a list, return the first element)
    surface_file = string_vs_list_check(surface_file)

    # Initialize list of label files and output colortable file
    label_files = []
    #relabel_file = os.path.join(os.getcwd(), 'relabel_annot.txt')
    #f_relabel = open(relabel_file, 'w')
    colortable = os.path.join(os.getcwd(), 'colortable.ctab')
    f_rgb = open(colortable, 'w')

    # Loop through labels
    irgb = 0
    for ilabel, label_number in enumerate(label_numbers):

        # Check type to make sure the number is an int
        label_number = int(label_number)
        label_name = label_names[ilabel]

        # Load surface
        reader = vtk.vtkDataSetReader()
        reader.SetFileName(surface_file)
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        d = data.GetPointData()
        labels = d.GetArray(scalar_name)

        # Write vertex index, coordinates, and 0
        count = 0
        npoints = data.GetNumberOfPoints()
        L = np.zeros((npoints,5))
        for i in range(npoints):
            label = labels.GetValue(i)
            if label == label_number:
                L[count,0] = i
                L[count,1:4] = data.GetPoint(i)
                count += 1

        # Save the label file
        if count > 0:
            irgb += 1

            # Write to relabel_file
            #if irgb != label_number:
            #    f_relabel.writelines('{0} {1}\n'.format(irgb, label_number))

            # Write to colortable
            f_rgb.writelines('{0} {1} {2}\n'.format(
                             irgb, label_name, RGBs[ilabel]))

            # Store in list of .label files
            label_file = hemi + '.' + label_name + '.label'
            label_file = os.path.join(os.getcwd(), label_file)
            label_files.append(label_file)

            # Write to .label file
            f = open(label_file, 'w')
            f.writelines('#!ascii label\n' + str(count) + '\n')
            for i in range(npoints):
                if any(L[i,:]):
                    pr = '{0} {1} {2} {3} 0\n'.format(
                         np.int(L[i,0]), L[i,1], L[i,2], L[i,3])
                    f.writelines(pr)
                else:
                    break
            f.close()
    f_rgb.close()
    #f_relabel.close()

    return label_files, colortable  #relabel_file


#if __name__ == "__main__" :
