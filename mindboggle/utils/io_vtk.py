#!/usr/bin/python
"""
Functions related to reading and writing VTK format files.

1. Functions for writing basic VTK elements
2. Functions for loading and writing VTK files
3. Functions specific to Mindboggle features


Authors:
Forrest Sheng Bao  .  http://fsbao.net
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

#=========================================
# Functions for writing basic VTK elements
#=========================================

def write_header(Fp, Title='', Header='# vtk DataFile Version 2.0',
                 FileType='ASCII', DataType='POLYDATA'):
    """
    Write header information for a VTK-format file.

    This part matches three things in the VTK 4.2 File Formats doc:
    Part 1: Header
    Part 2: Title (256 characters maximum, terminated with newline \n character)
    Part 3: Data type, either ASCII or BINARY
    Part 4: Geometry/topology. Type is one of:
        STRUCTURED_POINTS
        STRUCTURED_GRID
        UNSTRUCTURED_GRID
        POLYDATA
        RECTILINEAR_GRID
        FIELD
    """
    Fp.write('{}\n{}\n{}\nDATASET {}\n'.format(Header, Title, FileType, DataType))

def write_points(Fp, point_list, Type="float"):
    """
    Write coordinates of points, the POINTS section in DATASET POLYDATA section.
    """
    Fp.write('POINTS {} {}\n'.format(len(point_list), Type))

    for point in point_list:
        [R, A, S] = point
        Fp.write('{} {} {}\n'.format(R, A, S))

def write_faces(Fp, face_list, vertices_per_face=3):
    """
    Write vertices forming triangular meshes,
    the POLYGONS section in DATASET POLYDATA section.
    """
    if vertices_per_face == 3:
        face_name = 'POLYGONS '
    elif vertices_per_face == 2:
        face_name = 'LINES '
    else:
        print('ERROR: Unrecognized number of vertices per face')

    Fp.write('{} {} {}\n'.format(face_name, len(face_list),
             len(face_list) * (vertices_per_face + 1)))

    for face in face_list:
        [V0, V1, V2] = face
        Fp.write('{} {} {} {}\n'.format(vertices_per_face, V0, V1, V2))

def write_vertices(Fp, vertex_list):
    """
    Write vertices, the VERTICES section in DATASET POLYDATA section.
    """
    Fp.write('VERTICES {} {}\n{} '.format(
             len(vertex_list), len(vertex_list) + 1, len(vertex_list)))
    [Fp.write('{} '.format(i)) for i in vertex_list]
    Fp.write('\n')

def write_vertex_LUT(Fp, LUT, LUTName, at_LUT_begin=True):
    """
    Write per-VERTEX values as a scalar LUT into a VTK file.

    Input
    =====
    LUT:  list of floats
    at_LUT_begin: [Boolean] True if the first vertex LUT in a VTK file

    """
    if at_LUT_begin:
        Fp.write('POINT_DATA {}\n'.format(len(LUT)))
    Fp.write('SCALARS {} float\n'.format(LUTName))
    Fp.write('LOOKUP_TABLE {}\n'.format(LUTName))
    for Value in LUT:
        Fp.write('{}\n'.format(Value))
    Fp.write('\n')

#============================================
# Functions for loading and writing VTK files
#============================================

def load_VTK_vertex(Filename):
    """
    Load VERTICES from a VTK file, along with the scalar values.

    Note:  vertex extraction iterates from 1 to Vrts.GetSize(), rather than 0

    Input
    =====
    Filename : string
        The path/filename of a VTK format file.

    Output
    ======
    Vertexes : list of integers
        Each element is an ID (i.e., index) of a point defined in POINTS segment of the VTK file

    Scalars : list of floats
        Each element is a scalar value corresponding to a vertex

    """
    import vtk
    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(Filename)
    Reader.Update()

    Data = Reader.GetOutput()
    Vrts = Data.GetVerts()
    Vertexes = [Vrts.GetData().GetValue(i) for i in xrange(1, Vrts.GetSize())]

    PointData = Data.GetPointData()
    print "There are", Reader.GetNumberOfScalarsInFile(), "scalars in file", Filename
    print "Loading the scalar", Reader.GetScalarsNameInFile(0)
    ScalarsArray = PointData.GetArray(Reader.GetScalarsNameInFile(0))
    Scalars = [ScalarsArray.GetValue(i) for i in xrange(0, ScalarsArray.GetSize())]

    return Vertexes, Scalars

def load_VTK_line(Filename):
    """
    Load VERTICES from a VTK file, along with the scalar values.

    The line that extracts vertices from a VTK
    iterates from 1 to Vrts.GetSize(), rather than from 0.

    Input
    =====
    Filename : string
        The path/filename of a VTK format file.
    Output
    ======
    Vertexes : list of integers
        Each element is an ID (i.e., index) of a point defined in POINTS segment of the VTK file
    Scalars : list of floats
        Each element is a scalar value corresponding to a vertex
    """

    import vtk
    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(Filename)
    Reader.Update()

    Data = Reader.GetOutput()
    Lns = Data.GetLines()

    Lines  = [[Lns.GetData().GetValue(j) for j in xrange(i*3+1, i*3+3) ] for i in xrange(Data.GetNumberOfLines())]

    PointData = Data.GetPointData()
    print "There are", Reader.GetNumberOfScalarsInFile(), "scalars in file", Filename
    print "Loading the scalar", Reader.GetScalarsNameInFile(0)
    ScalarsArray = PointData.GetArray(Reader.GetScalarsNameInFile(0))
    Scalars = [ScalarsArray.GetValue(i) for i in xrange(0, ScalarsArray.GetSize())]

    return Lines, Scalars

def write_scalars(vtk_file, Points, Vertices, Faces, LUTs=[], LUT_names=[]):
    """
    Save scalars into a VTK-format file.

    Inputs:
    ======
    vtk_file : string
        The path of the VTK file to save sulci
    Points :  list of 3-tuples of floats
        Each element has 3 numbers representing the coordinates of the points
    Vertices: list of integers
        IDs of vertices that are part of a sulcus
    Faces: list of 3-tuples of integers
        Each element is a face on the mesh, consisting of 3 integers
        representing the 3 vertices of the face
    LUTs: list of lists of integers
        Each element is a list of integers representing a scalar map for the mesh
    LUT_names: list of strings
        Each element is the name of a scalar map, e.g., curv, depth.

    Example:
    =======
    import random
    import io_vtk
    Points = [[random.random() for i in [1,2,3]] for j in xrange(0,4)]
    Vertices = [1,2,3,0]
    Faces = [[1,2,3],[0,1,3]]
    LUT_names = ['curv','depth']
    LUTs=[[random.random() for i in xrange(1,5)] for j in [1,2]]
    io_vtk.write_scalars('test.vtk',Points, Vertices, Faces, LUTs=LUTs,
                      LUT_names=LUT_names)

    """
    import os
    from utils import io_vtk

    vtk_file = os.path.join(os.getcwd(), vtk_file)

    Fp = open(vtk_file,'w')
    io_vtk.write_header(Fp)
    io_vtk.write_points(Fp, Points)
    io_vtk.write_vertices(Fp, Vertices)
    io_vtk.write_faces(Fp, Faces)
    if len(LUTs) > 0:
        # Make sure that LUTs is a list of lists
        if type(LUTs[0]) != list:
            LUTs = [LUTs]
        for i, LUT in enumerate(LUTs):
            if i == 0:
                io_vtk.write_vertex_LUT(Fp, LUT, LUT_names[i])
            else:
                io_vtk.write_vertex_LUT(Fp, LUT, LUT_names[i],
                                        at_LUT_begin=False)
    Fp.close()

    return vtk_file

def write_scalar_subset(input_vtk, output_vtk, new_scalars, filter_scalars=[]):
    """
    Load VTK format file and save a subset of (non-zero) scalars into a new file.

    Inputs:
    ======
    input_vtk:  input VTK file [string]
    output_vtk:  output VTK file [string]
    new_scalars:  new scalar values for VTK file
    filter_scalars:  (optional)
                     scalar values used to filter faces (non-zero are retained)

    """
    import os
    from utils import io_vtk

    # Load VTK file
    Points, Faces, Scalars = io_vtk.load_scalar(input_vtk, return_arrays=1)

    # Find indices to nonzero values
    indices = range(len(Scalars))
    if len(filter_scalars) > 0:
        indices_nonzero = [i for i,x in enumerate(filter_scalars) if int(x) > 0]
    else:
        indices_nonzero = indices
    # Remove surface mesh faces whose three vertices are not all in indices
    faces_indices = io_vtk.inside_faces(Faces, indices_nonzero)

    # Lookup lists for saving to VTK format files
    LUTs = [new_scalars]
    LUT_names = ['scalars']

    # Output VTK file to current working directory
    output_vtk = os.path.join(os.getcwd(), output_vtk)

    # Write VTK file
    Fp = open(output_vtk,'w')
    io_vtk.write_header(Fp)
    io_vtk.write_points(Fp, Points)
    io_vtk.write_vertices(Fp, indices)
    io_vtk.write_faces(Fp, faces_indices)
    if len(LUTs) > 0:
        for i, LUT in enumerate(LUTs):
            if i == 0:
                io_vtk.write_vertex_LUT(Fp, LUT, LUT_names[i])
            else:
                io_vtk.write_vertex_LUT(Fp, LUT, LUT_names[i],
                                        at_LUT_begin=False)
    Fp.close()

    return output_vtk

def write_mean_scalar_table(filename, column_names, labels, *shape_files):
    """
    Make a table of mean values per label per measure.

    Inputs:
    ======
    filename:  output filename (without path)
    column_names:  names of columns [list of strings]
    labels:  list (same length as values)
    *shape_files:  arbitrary number of vtk files with scalar values

    Output:
    ======
    filename:  table file name

    """
    import os
    from utils.io_vtk import load_scalar
    from utils.io_file import write_table
    from measure.measure_functions import mean_value_per_label

    columns = []
    for shape_file in shape_files:

        Points, Faces, Scalars = load_scalar(shape_file, return_arrays=1)

        mean_values, label_list = mean_value_per_label(Scalars, labels)

        columns.append(mean_values)

    filename = os.path.join(os.getcwd(), filename)
    write_table(label_list, columns, column_names, filename)

    return filename

def load_scalar(filename, return_arrays=1):
    """
    Load a VTK-format scalar map that contains only one SCALAR segment.

    Inputs
    =======

    filename : string
        The path/filename of a VTK format file.
    return_arrays: return numpy arrays instead of lists of lists below (1=yes, 0=no)

    Outputs
    =========
    Points : list of lists of floats (see return_arrays)
        Each element is a list of 3-D coordinates of a vertex on a surface mesh

    Faces : list of lists of integers (see return_arrays)
        Each element is list of 3 IDs of vertices that form a face
        on a surface mesh

    Scalars : list of floats (see return_arrays)
        Each element is a scalar value corresponding to a vertex

    """
    import vtk
    if return_arrays:
        import numpy as np

    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()

    Data = Reader.GetOutput()
    Points = [list(Data.GetPoint(point_id))
              for point_id in xrange(0, Data.GetNumberOfPoints())]

    CellArray = Data.GetPolys()
    Polygons = CellArray.GetData()
    Faces = [[Polygons.GetValue(j) for j in xrange(i*4 + 1, i*4 + 4)]
             for i in xrange(0, CellArray.GetNumberOfCells())]

    PointData = Data.GetPointData()
    print("Loading {} {} scalars in file {}...".
          format(Reader.GetNumberOfScalarsInFile,
                 Reader.GetScalarsNameInFile(0), filename))
    ScalarsArray = PointData.GetArray(Reader.GetScalarsNameInFile(0))
    if ScalarsArray:
        Scalars = [ScalarsArray.GetValue(i) for i in xrange(0, ScalarsArray.GetSize())]
    else:
        Scalars = []

    if return_arrays:
        return np.array(Points), np.array(Faces), np.array(Scalars)
    else:
        return Points, Faces, Scalars

def inside_faces(faces, indices):
    """
    Remove surface mesh faces whose three vertices are not all in "indices"

    Inputs:
    ======
    faces: triangular surface mesh vertex indices [#faces x 3]
    indices: vertex indices to mesh

    Return:
    ======
    faces: reduced array of faces

    """
    import numpy as np

    len_faces = len(faces)
    fs = frozenset(indices)
    faces = [lst for lst in faces if len(fs.intersection(lst)) == 3]
    faces = np.reshape(np.ravel(faces), (-1, 3))
    print('  Reduced {} to {} triangular faces.'.format(len_faces, len(faces)))

    return faces



"""
def write_mean_scalar_table(filename):

    import numpy as np

    Points, Faces, Scalars = load_scalar(filename, return_arrays=1)

    unique_scalars = np.unique(Scalars)

    Table = []

    for scalar in unique_scalars:
        Index = [i for i in xrange(len(Label)) if Label[i] == GyralLabel]
        Measure_of_Row = [ [Measure[i] for i in Index ] for Measure in Measures]
        Row = map(np.mean, Measure_of_Row)
        Row = [GyralLabel] + Row
        Table.append(Row)

    return Table
"""

#==========================================
# Functions specific to Mindboggle features
#==========================================

def write_feature_to_face(Fp, Vertex, Face, index_pair_list):
    """
    Load IDs of faces (in my output format) and original surface file
    to output a face from a feature in VTK format at STDIO.
    """

    from utils import io_vtk

    io_vtk.write_header(Fp, 'write_feature_to_face() output by Forrest Bao')
    io_vtk.write_points(Fp, Vertex)

    Fundi = []
    for i in xrange(0, len(index_pair_list)):
        Fundi.append(Face[index_pair_list[i]])

    io_vtk.write_faces(Fp, Fundi)

def write_vertices_to_fundi(Fp, Vertex, index_pair_list):
    """
    Load IDs of fundus vertices (in my output format) and original surface file
    to output fundi in VTK format.
    """

    from utils import io_vtk

    io_vtk.write_header(Fp, 'write_vertices_to_fundi() output by Forrest Bao')
    io_vtk.write_points(Fp, Vertex)

    io_vtk.write_vertices(Fp, index_pair_list)

def write_line_segments_to_fundi(Fp, Vertex, index_pair_list):
    """
    Load IDs of fundus curve segments (in my output format)
    and original surface file to output fundi in VTK format.
    """

    from utils import io_vtk

    io_vtk.write_header(Fp, 'Created by Mindboggle')
    io_vtk.write_points(Fp, Vertex)
    io_vtk.write_faces(Fp, index_pair_list, vertices_per_face=2)

def load_fundi_list(filename):
    """
    Load the file storing face/vertex IDs, which are fundi faces/vertices.
    """
    f = open(filename, 'r')
    lines = f.readlines()
    for i in xrange(0,len(lines)):
        lines[i] = int(lines[i][:-1])
    return lines

def load_segmented_fundi(filename):
    """
    Load the file storing fundi as curve segments.
    """

    Fp = open(filename, 'r')
    lines = Fp.readlines()
    Fp.close()

    Segs = []
    for line in lines:
        Segs.append(line)
        # I do NOT convert strings to integers because
        # later we write strings into VTK files.
        # Also no need to split. So break line is also included.

    return Segs

def write_fundi(vtk_file, Points, Vertices, Lines, LUTs=[], LUT_names=[]):
    """
    Save fundi into VTK files

    Parameters
    =============

    vtk_file : string
        The path of the VTK file to save fundi

    Points :  list of 3-tuples of floats
        Each element has 3 numbers representing the coordinates of the points

    Vertices: list of integers
        IDs of vertices that are part of a fundus

    Lines: list of 2-tuples of integers
        Each element is an edge on the mesh, consisting of 2 integers
        representing the 2 vertices of the edge

    LUTs: list of lists of integers
        Each element is a list of integers representing a map for the mesh

    LUT_names: list of strings
        Each element is the name of a map, e.g., curv, depth.

    Example
    ===========
    import random
    from utils import io_vtk
    Points = [[random.random() for i in range(3)] for j in xrange(5)]
    Vertices = [0,1,2,3,4]
    Faces = [[1,2,3],[0,3,4]]
    LUT_names = ['curv','depth']
    LUTs=[[random.random() for i in range(6)] for j in range(2)]
    io_vtk.write_scalars('test.vtk',Points, Vertices, Faces, LUTs=LUTs,
                      LUT_names=LUT_names)
    """

    import os
    from utils import io_vtk

    vtk_file = os.path.join(os.getcwd(), vtk_file)

    Fp = open(vtk_file,'w')
    io_vtk.write_header(Fp)
    io_vtk.write_points(Fp, Points)
    io_vtk.write_vertices(Fp, Vertices)
    for i in xrange(0,len(Lines)):
        Lines[i] = str(Lines[i][0]) + " " + str(Lines[i][1]) + "\n"
    io_vtk.write_faces(Fp, Lines, vertices_per_face=2)
    if len(LUTs) > 0:
        for i, LUT in enumerate(LUTs):
            if i == 0:
                io_vtk.write_vertex_LUT(Fp, LUT, LUT_names[i])
            else:
                io_vtk.write_vertex_LUT(Fp, LUT, LUT_names[i],
                                        at_LUT_begin=False)
    Fp.close()
