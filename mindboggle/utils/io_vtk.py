#!/usr/bin/env python
"""
Functions related to reading and writing VTK format files.

1. Functions for reading VTK files
2. Functions for writing VTK files
3. Functions for converting FreeSurfer files to VTK format


Authors:
    - Forrest Sheng Bao, 2012-2013  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

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
    print("There are {0} scalars in file {1}".format(
        Reader.GetNumberOfscalarsInFile(), Filename))
    print("Loading the scalar {0}".format(Reader.GetScalarsNameInFile(0)))
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
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
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
            scalar = [scalar_array.GetValue(i) for i in range(scalar_array.GetSize())]
            scalars.append(scalar)
            scalar_names.append(scalar_name)

    if return_first:
        if scalars:
            scalars = scalars[0]
        if return_array:
            scalars = np.array(scalars)
        if scalar_names:
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
    scalars : list or list of lists of floats or integers
        scalar values for the vertices of a mesh
    scalar_names : string or list of strings
        name(s) of lookup table(s)

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
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
        if scalars:
            scalars = scalars[0]
        if return_array:
            scalars = np.array(scalars)
        if scalar_names:
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
    Fp : pointer to a file
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
    scalars : list of lists of floats (or single list of floats)
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
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> points = [[random.random() for i in [1,2,3]] for j in xrange(4)]
    >>> indices = [1,2,3,0]
    >>> lines = [[1,2],[3,4]]
    >>> faces = [[1,2,3],[0,1,3]]
    >>> scalar_names = ['curv','depth']
    >>> scalars = [[random.random() for i in xrange(4)] for j in [1,2]]
    >>> #
    >>> write_vtk('test_write_vtk.vtk', points,
    >>>          indices, lines, faces, scalars, scalar_names)
    >>> #
    >>> # View:
    >>> plot_vtk('test_write_vtk.vtk')
    >>> #
    >>> # Write vtk file with depth values on sulci and view:
    >>> from mindboggle.utils.io_vtk import read_vtk, write_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file)
    >>> write_vtk('test_write_vtk.vtk', points, [], [], faces, depths, 'depths')
    >>> plot_vtk('test_write_vtk.vtk')

    """
    import os
    from mindboggle.utils.io_vtk import write_header, write_points, \
         write_vertices, write_faces, write_scalars, \
         scalars_checker

    output_vtk = os.path.join(os.getcwd(), output_vtk)

    Fp = open(output_vtk,'w')
    write_header(Fp)
    write_points(Fp, points)

    if indices:
        write_vertices(Fp, indices)
    if lines:
        for i in xrange(0,len(lines)):
            lines[i] = [lines[i][0], lines[i][1]]
        write_faces(Fp, lines) # write_faces can write either lines or faces
    if faces:
        write_faces(Fp, faces)
    if scalars:
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
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> depths, name = read_scalars(depth_file, True,True)
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> sulci, name = read_scalars(sulci_file)
    >>> #
    >>> rewrite_scalars(depth_file, 'test_rewrite_scalars.vtk',
    >>>                 [depths, sulci], ['depths', 'sulci'], sulci)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_rewrite_scalars.vtk')

    """
    import os
    import numpy as np

    from mindboggle.utils.mesh import remove_faces
    from mindboggle.utils.io_vtk import write_header, write_points, \
         write_vertices, write_faces, write_scalars, read_vtk, \
         scalars_checker

    # Convert numpy arrays to lists
    if isinstance(new_scalars, np.ndarray):
        new_scalars = new_scalars.tolist()
    if isinstance(filter_scalars, np.ndarray):
        filter_scalars = filter_scalars.tolist()

    # Output VTK file to current working directory
    output_vtk = os.path.join(os.getcwd(), output_vtk)

    # Load VTK file
    faces, lines, indices, points, npoints, scalars, name = read_vtk(input_vtk)

    # Find indices to nonzero values
    indices = range(npoints)
    if filter_scalars:
        indices_filter = [i for i,x in enumerate(filter_scalars) if x > -1]
        indices_remove = [i for i,x in enumerate(filter_scalars) if x == -1]
        # Remove surface mesh faces whose three vertices are not all in indices
        faces = remove_faces(faces, indices_filter)

    # Write VTK file
    Fp = open(output_vtk,'w')
    write_header(Fp)
    write_points(Fp, points)
    if indices:
        write_vertices(Fp, indices)
    if faces:
        write_faces(Fp, faces)
    if new_scalars:
        new_scalars, new_scalar_names = scalars_checker(new_scalars, new_scalar_names)

        for i, new_scalar_list in enumerate(new_scalars):
            if filter_scalars:
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
    >>> from mindboggle.utils.io_vtk import explode_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> output_stem = 'sulcus'
    >>> #
    >>> explode_scalars(sulci_file, output_stem)
    >>> #
    >>> # View:
    >>> example_vtk = os.path.join(os.getcwd(), output_stem + '0.vtk')
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk(example_vtk)

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
    Check whether input scalars and scalar_names are in acceptable format.
    If not, reformat.

    Parameters
    ----------
    scalars : list of lists of floats (or single list or 1-/2-D array of floats)
    scalar_names : string or list of strings

    Returns
    -------
    scalars : list of lists of floats
    scalar_names : list of strings

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import scalars_checker
    >>> scalars_checker([[1,2],[3,4]], ["list1", "list2"])
      ([[1, 2], [3, 4]], ['list1', 'list2'])
    >>> scalars_checker([[1,2],[3,4]], "")
      ([[1, 2], [3, 4]], ['', ''])
    >>> scalars_checker([[1,2],[3,4]], ["list1", "list2", "list3"])
      ([[1, 2], [3, 4]], ['list1', 'list2', 'list3'])
    >>> scalars_checker([[1,2],[3,4]], "list1")
      ([[1, 2], [3, 4]], ['list1', 'list1'])
    >>> scalars_checker([1,2,3,4], ["123"])
      ([[1, 2, 3, 4]], ['123'])
    >>> scalars_checker(1, ["123"])
      Error: scalars is neither a list nor a numpy array.
    >>> scalars_checker(np.array([1,2,3]), ["123"])
      ([[1, 2, 3]], ['123'])
    >>> scalars_checker(np.array([[1,2,3]]), ["123"])
      ([[1, 2, 3]], ['123'])
    >>> scalars_checker(np.array([[1,2,3],[4,5,6]]), ["123"])
      ([[1, 2, 3], [4, 5, 6]], ['123', '123'])
    >>> scalars_checker(np.array([[[1,2,3]]]), ["123"])
      Error: Dimension of new_scalars is too high.
    >>> scalars_checker(np.array([np.array([0,7,9]),[1,2,3]]), ["123"])
      ([[0, 7, 9], [1, 2, 3]], ['123', '123'])

    Notes
    -----
    This function does not check all possible cases of scalars and scalar_names,
    but only those that are likely to happen when using Mindboggle.

    """
    import sys
    import numpy as np

    # If not a list, convert to a list.
    if not isinstance(scalars, list):
        if isinstance(scalars, np.ndarray):
            if len(scalars.shape) < 2: # this is at most a 1-D array
                scalars = [scalars.tolist()]
            elif len(scalars.shape) == 2: # 2-D numpy array
                scalars = scalars.tolist()
            else:
                print("Error: Dimension of new_scalars is too high.")
                sys.exit()
        else:
            print("Error: scalars is neither a list nor a numpy array.")
            sys.exit()

    # If the list contains integers or floats, put in a list.
    if isinstance(scalars[0], int) or isinstance(scalars[0], float):
        scalars = [scalars]
    # If the list contains all lists, accept format.
    elif all([isinstance(x, list) for x in scalars]):
        pass
    # If the list contains arrays (and optionally lists), convert arrays to lists.
    elif all([isinstance(x, list) or isinstance(x, np.ndarray) for x in scalars]):
        scalars2 = []
        for x in scalars:
            if isinstance(x, list):
                scalars2.append(x)
            else:
                scalars2.append(x.tolist())
        scalars = scalars2
    else:
        print("Error: scalars is a 1-D list containing unacceptable elements.")
        print("scalars type is: {0}".format(type(scalars)))
        print("scalars length is: {0}".format(len(scalars)))
        print("scalars[0] type is: {0}".format(type(scalars[0])))
        sys.exit()

    # If scalar_names is a string, create a list containing
    # as many of this string as there are scalar lists.
    if isinstance(scalar_names, str):
        scalar_names = [scalar_names for x in scalars]
    elif isinstance(scalar_names, list):
        if len(scalar_names) < len(scalars):
            scalar_names = [scalar_names[0] for x in scalars]
        else:
            pass
    else:
        print("Error: scalar_names is neither a list nor a string")
        sys.exit()

    return scalars, scalar_names


#if __name__ == "__main__" :
