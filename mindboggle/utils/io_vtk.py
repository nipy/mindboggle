#!/usr/bin/python
"""
Functions related to reading and writing VTK format files.

This is not part of the official vtk module.
We found that VTK's Python binding has limited documentation.

Authors:
Forrest Sheng Bao  .  http://fsbao.net
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import io_file

#=============================================
# LEVEL 1: Basic VTK element writing functions
#=============================================

def write_header(Fp, Title='', Header='# vtk DataFile Version 2.0',
                 FileType='ASCII', DataType='POLYDATA'):
    """
    Write header information for a VTK-format file.

    This part matches three things in the VTK 4.2 File Formats doc

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

    Fp.write(Header)
    Fp.write("\n")
    Fp.write(Title)
    Fp.write("\n")
    Fp.write(FileType)
    Fp.write("\n")
    Fp.write("DATASET ")
    Fp.write(DataType)
    Fp.write("\n")

def write_points(Fp, point_list, Type="float"):
    """
    Write coordinates of points, the POINTS section in DATASET POLYDATA section.
    """
    Fp.write("POINTS " + str(len(point_list)) + " " + Type + "\n")
    for i in xrange(0, len(point_list)):
        [R, A, S] = point_list[i]
        Fp.write(str(R) + " " + str(A) + " " + str(S) + "\n")

def write_faces(Fp, face_list, vertices_per_face=3):
    """
    Write vertices forming triangular meshes,
    the POLYGONS section in DATASET POLYDATA section.
    """
    Fp.write("POLYGONS " + str(len(face_list)) + " " +
             str( (vertices_per_face + 1) * len(face_list)  )  + '\n' )
    for i in xrange(0, len(face_list)):
        [V0, V1, V2] = face_list[i]
        Fp.write( str(vertices_per_face) + " " + str(V0) + " " +
                  str(V1) + " " + str(V2) + "\n")

def write_vertices(Fp, vertex_list):
    """
    Write vertices, the VERTICES section in DATASET POLYDATA section.
    """
    # One possible solution
    Fp.write("VERTICES " + str(len(vertex_list)) + " " + str(len(vertex_list)+1) +
             "\n" + str(len(vertex_list)) + " ")
    [Fp.write(str(i)+" ") for i in vertex_list]
    Fp.write("\n")

def write_line_segments(Fp, index_pair_list):
    """
    Write line segments, each with two end points, into VTK format.

    Parameters
    ===========

    index_pair_list: list of strings (NOT integers)
       each element of index_pair_list is a string containing IDs of two vertices, like "1 3"

    """
    Fp.write("LINES " + str(len(index_pair_list)) + " " + str(len(index_pair_list)*3) + "\n")
    [Fp.write("2 " + Vrtx) for Vrtx in index_pair_list]

def write_vertex_LUT(Fp, LUT, LUTName, at_LUT_begin=True):
    """
    Write per-VERTEX values as a scalar LUT into a VTK file.

    This function is called by face_list_to_vtk

    Parameters
    ==========

    LUT    : list of floats

    at_LUT_begin: Boolean
        True, if this vertex LUT is the first vertex LUT in a VTK file.

    """
    if at_LUT_begin:
        Fp.write('POINT_DATA ' + str(len(LUT)) +'\n')
    Fp.write('SCALARS '+ LUTName +' float\n')
    Fp.write('LOOKUP_TABLE '+ LUTName +'\n')
    for Value in LUT:
        Fp.write(str(Value) + '\n')
    Fp.write('\n')


#========================================
# LEVEL 2: Saving features into VTK files
#========================================

def write_scalars(vtk_file, Points, Vertices, Faces, LUTs=[], LUT_names=[]):
    """
    Save scalars into a VTK-format file.

    Parameters
    =============

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

    Example
    ===========
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
    Fp = open(vtk_file,'w')
    write_header(Fp)
    write_points(Fp, Points)
    write_vertices(Fp, Vertices)
    write_faces(Fp, Faces)
    if len(LUTs) > 0:
        # Make sure that LUTs is a list of lists
        if type(LUTs[0]) != list:
            LUTs = [LUTs]
        for i, LUT in enumerate(LUTs):
            if i == 0:
                write_vertex_LUT(Fp, LUT, LUT_names[i])
            else:
                write_vertex_LUT(Fp, LUT, LUT_names[i], at_LUT_begin=False)
    Fp.close()

    return vtk_file

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
    import io_vtk
    Points = [[random.random() for i in [1,2,3]] for j in xrange(0,4)]
    Vertices = [1,2,3,0]
    Liness = [[1,2],[0,3]]
    LUT_names = ['curv','depth']
    LUTs=[[random.random() for i in xrange(1,5)] for j in [1,2]]
    io_vtk.write_scalars('test.vtk',Points, Vertices, Lines, LUTs=LUTs,
                      LUT_names=LUT_names)
    """
    Fp = open(vtk_file,'w')
    write_header(Fp)
    write_points(Fp, Points)
    write_vertices(Fp, Vertices)
    for i in xrange(0,len(Lines)):
        Lines[i] = str(Lines[i][0]) + " " + str(Lines[i][1]) + "\n"
    write_line_segments(Fp, Lines)
    if len(LUTs) > 0:
        for i, LUT in enumerate(LUTs):
            if i == 0:
                write_vertex_LUT(Fp, LUT, LUT_names[i])
            else:
                write_vertex_LUT(Fp, LUT, LUT_names[i], at_LUT_begin=False)
    Fp.close()

def surf_to_vtk(surface_file):
    Vertex, Face = io_file.readSurf(surface_file)
    Fp = open(surface_file + '.vtk', 'w')
    write_header(Fp, Title='vtk output from'+surface_file)
    write_points(Fp, Vertex)
    write_faces(Fp, Face)
    Fp.close()

def face_list_to_vtk(vtk_file, surface_file, index_pair_file, LUT=[], LUTname=[]):
    """
    Load a face list file and a surface file to map faces onto the surface
    and save the result into VTK format.

    This function is called by libbasin.getBasin()

    """
    Fp = open(vtk_file,'w')
    Vertex, Face = io_file.readSurf(surface_file)
    index_pair_list = load_fundi_list(index_pair_file)
    write_feature_to_face(Fp, Vertex, Face, index_pair_list)

    if LUT!=[] :
        for i in xrange(0, len(LUT)):
            if i == 0:
                write_vertex_LUT(Fp, LUT[i], LUTname[i])
            else:
                write_vertex_LUT(Fp, LUT[i], LUTname[i], at_LUT_begin=False)

    Fp.close()

def vertex_list_to_vtk(vtk_file, surface_file, index_pair_file, LUT=[], LUTname=[]):
    """
    Load a vertex list file and a surface file to map vertices onto the surface
    and save the result into VTK format.

    Parameters
    ============

    LUT    : list of LUTs
        LUT[i] is the i-th LUT, e.g., distance transformation values,
                                      curvature, convexity values.
        LUT[i] is a list of float numbers
        So more than one LUTs may be written to the final VTK file

    LUTname    : list of strings
        LUTname[i] is the name of the lookup table to be inserted into VTK file

    """
    Fp = open(vtk_file,'w')
    Vertex, Face = io_file.readSurf(surface_file)
    index_pair_list = load_fundi_list(index_pair_file)
    write_vertices_to_fundi(Fp, Vertex, index_pair_list)
    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            if i == 0:
                write_vertex_LUT(Fp, LUT[i], LUTname[i])
            else:
                write_vertex_LUT(Fp, LUT[i], LUTname[i], at_LUT_begin=False)

    Fp.close()

def line_segments_to_vtk(vtk_file, surface_file, index_pair_file, LUT=[], LUTname=[]):
    """
    Load a fundus curve segment list file and a surface file
    to map curve segments onto the surface and save the result into VTK format.

    Parameters
    ============

    LUT    : list of LUTs
        LUT[i] is the i-th LUT, e.g., distance transformation values,
                                      curvature, convexity values.
        LUT[i] is a list of float numbers
        So more than one LUTs may be written to the final VTK file

    LUTname    : list of strings
        LUTname[i] is the name of the lookup table to be inserted into VTK file

    """
    Fp = open(vtk_file,'w')
    Vertex, Face = io_file.readSurf(surface_file)
    index_pair_list = load_segmented_fundi(index_pair_file)
    write_line_segments_to_fundi(Fp, Vertex, index_pair_list)

    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            if i == 0:
                write_vertex_LUT(Fp, LUT[i], LUTname[i])
            else:
                write_vertex_LUT(Fp, LUT[i], LUTname[i], at_LUT_begin=False)

    Fp.close()


#=============================================================================
# LEVEL 3: Functions for loading and writing called by LEVEL 2 functions above
#=============================================================================

def write_line_segments(Fp, index_pair_list):
    """
    Write a line segment with two end points into VTK format.

    Parameters
    ===========

    index_pair_list: list of strings (NOT integers)
       each element of index_pair_list is a string containing IDs of two vertices,
       like "1 3"

    """
    Fp.write("LINES " + str(len(index_pair_list)) + " " + str(len(index_pair_list)*3) + "\n")
    [Fp.write("2 " + Vrtx) for Vrtx in index_pair_list]

def load_fundi_list(filename):
    """
    Load the file storing face/vertex IDs, which are fundi faces/vertices.
    """
    f = open(filename, 'r')
    lines = f.readlines()
    for i in xrange(0,len(lines)):
        lines[i] = int(lines[i][:-1])
    return lines

def write_feature_to_face(Fp, Vertex, Face, index_pair_list):
    """
    Load IDs of faces (in my output format) and original surface file
    to output a face from a feature in VTK format at STDIO.
    """
    write_header(Fp, 'write_feature_to_face() output by Forrest Bao')
    write_points(Fp, Vertex)

    Fundi = []
    for i in xrange(0, len(index_pair_list)):
        Fundi.append(Face[index_pair_list[i]])

    write_faces(Fp, Fundi)

def write_vertices_to_fundi(Fp, Vertex, index_pair_list):
    """
    Load IDs of fundus vertices (in my output format) and original surface file
    to output fundi in VTK format.
    """
    write_header(Fp, 'write_vertices_to_fundi() output by Forrest Bao')
    write_points(Fp, Vertex)

    write_vertices(Fp, index_pair_list)

def write_line_segments_to_fundi(Fp, Vertex, index_pair_list):
    """
    Load IDs of fundus curve segments (in my output format)
    and original surface file to output fundi in VTK format.
    """

    write_header(Fp, 'Created by Mindboggle')
    write_points(Fp, Vertex)
    write_line_segments(Fp, index_pair_list)

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


#----------------------------------------------------------
# Functions that use VTK's python binding to load VTK files
#----------------------------------------------------------

def load_scalar(filename):
    """
    Load a VTK-format scalar map that contains only one SCALAR segment.

    Inputs
    =======

    filename : string
        The path/filename of a VTK format file.


    Outputs
    =========
    Points : list of lists of floats
        Each element is a list of 3-D coordinates of a vertex on a surface mesh

    Faces : list of lists of integers
        Each element is a 3-tuple, the IDs of vertices that form a face
        on a surface mesh

    Scalars : list of floats
        Each element is a scalar value corresponding to a vertex

    """
    import vtk
    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()

    Data = Reader.GetOutput()
    Points = [list(Data.GetPoint(PointID))
              for point_id in xrange(0, Data.GetNumberOfPoints())]

    Polys = Data.GetPolys()
    CellArray = Polys.GetData()
    Faces = [[CellArray.GetValue(j) for j in xrange(i*4 + 1, i*4 + 4)]
             for i in xrange(0,  CellArray.GetNumberOfCells()) ]

    PointData = Data.GetPointData()
    print("Loading {} {} scalars in file {}...".
          format(Reader.GetNumberOfScalarsInFile,
                 Reader.GetScalarsNameInFile(0), filename))
    ScalarsArray = PointData.GetArray(Reader.GetScalarsNameInFile(0))
    Scalars = [ScalarsArray.GetValue(i) for i in xrange(0, ScalarsArray.GetSize())]

    return Points, Faces, Scalars



