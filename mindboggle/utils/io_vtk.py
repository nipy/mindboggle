#!/usr/bin/python
"""
Functions related to reading and writing VTK format files.

This is not the official vtk module. One reason that we do not use VTK's
python binding is the lack of documentation.

Authors:  Forrest Sheng Bao http://fsbao.net
Version:  0.2, last update on 2012-06-29

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import io_file


# Level 1: basic VTK element writing functions ----

def writeHeader(Fp, Title='', Header='# vtk DataFile Version 2.0',
                FileType='ASCII', DataType='POLYDATA'):
    """Write the all non-data information for a VTK-format file

    This part matches three things in VTK 4.2 File Formats doc

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

def writePoint(Fp, PointList, Type="float"):
    """
    Write coordinates of points, the POINTS section in DATASET POLYDATA section
    """
    Fp.write("POINTS " + str(len(PointList)) + " " + Type + "\n")
    for i in xrange(0, len(PointList)):
        [R, A, S] = PointList[i]
        Fp.write(str(R) + " " + str(A) + " " + str(S) + "\n")

def writeFace(Fp, FaceList, VertexPerFace=3):
    """
    Write vertexes forming triangular meshes,
    the POLYGONS section in DATASET POLYDATA section
    """
    Fp.write("POLYGONS " + str(len(FaceList)) + " " +
             str( (VertexPerFace + 1) * len(FaceList)  )  + '\n' )
    for i in xrange(0, len(FaceList)):
        [V0, V1, V2] = FaceList[i]
        Fp.write( str(VertexPerFace) + " " + str(V0) + " " +
                  str(V1) + " " + str(V2) + "\n")

def writeVrtx(Fp, VrtxList):
    """Write vertexes, the VERTICES section in DATASET POLYDATA section
    """
    # One possible solution
    Fp.write("VERTICES " + str(len(VrtxList)) + " " + str(len(VrtxList)+1) +
             "\n" + str(len(VrtxList)) + " ")
    [Fp.write(str(i)+" ") for i in VrtxList]
    Fp.write("\n")

def writeSeg(Fp, FundiList):
    """Write a line segment with two ending points into VTK format

    Parameters
    ===========

    FundiList: list of strings (NOT integers)
       each element of FundiList is a string containing IDs of two vertexes, like "1 3"

    """
    Fp.write("LINES " + str(len(FundiList)) + " " + str(len(FundiList)*3) + "\n")
    [Fp.write("2 " + Vrtx) for Vrtx in FundiList]

def wrtVrtxLUT(Fp, LUT, LUTName, AtLUTBegin=True):
    """write per-VERTEX values as a scalar LUT into a VTK file

    This function is called by fcLst2VTK

    Parameters
    ==========

    LUT    : list of floats

    AtLUTBegin: Boolean
        True, if this vertex LUT is the first vertex LUT in a VTK file.

    """
    if AtLUTBegin:
        Fp.write('POINT_DATA ' + str(len(LUT)) +'\n')
    Fp.write('SCALARS '+ LUTName +' float\n')
    Fp.write('LOOKUP_TABLE '+ LUTName +'\n')
    for Value in LUT:
        Fp.write(str(Value) + '\n')
    Fp.write('\n')

# end of Level 1 functions  ---

# Level 2: saving features into VTK files ---

def writeSulci(VTKFile, Points, Vertexes, Faces, LUTs=[], LUTNames=[]):
    """
    Save scalars into a VTK-format file.

    Parameters
    =============

    VTKFile : string
        The path of the VTK file to save sulci

    Points :  list of 3-tuples of floats
        Each element has 3 numbers representing the coordinates of the points

    Vertexes: list of integers
        IDs of vertices that are part of a sulcus

    Faces: list of 3-tuples of integers
        Each element is a face on the mesh, consisting of 3 integers
        representing the 3 vertices of the face

    LUTs: list of lists of integers
        Each element is a list of integers representing a scalar map for the mesh

    LUTNames: list of strings
        Each element is the name of a scalar map, e.g., curv, depth.

    Example
    ===========
    import random
    import io_vtk
    Points = [[random.random() for i in [1,2,3]] for j in xrange(0,4)]
    Vertexes = [1,2,3,0]
    Faces = [[1,2,3],[0,1,3]]
    LUTNames = ['curv','depth']
    LUTs=[[random.random() for i in xrange(1,5)] for j in [1,2]]
    io_vtk.writeSulci('test.vtk',Points, Vertexes, Faces, LUTs=LUTs,
                      LUTNames=LUTNames)

    """
    Fp = open(VTKFile,'w')
    writeHeader(Fp)
    writePoint(Fp, Points)
    writeVrtx(Fp, Vertexes)
    writeFace(Fp, Faces)
    if len(LUTs) > 0:
        for i, LUT in enumerate(LUTs):
            if i == 0:
                wrtVrtxLUT(Fp, LUT, LUTNames[i])
            else:
                wrtVrtxLUT(Fp, LUT, LUTNames[i], AtLUTBegin=False)
    Fp.close()

def writeFundi(VTKFile, Points, Vertexes, Lines, LUTs=[], LUTNames=[]):
    """
    Save fundi into VTK files

    Parameters
    =============

    VTKFile : string
        The path of the VTK file to save fundi

    Points :  list of 3-tuples of floats
        Each element has 3 numbers representing the coordinates of the points

    Vertexes: list of integers
        IDs of vertexes that are part of a fundus

    Lines: list of 2-tuples of integers
        Each element is an edge on the mesh, consisting of 2 integers
        representing the 2 vertexes of the edge

    LUTs: list of lists of integers
        Each element is a list of integers representing a map for the mesh

    LUTNames: list of strings
        Each element is the name of a map, e.g., curv, depth.

    Example
    ===========
    import random
    import io_vtk
    Points = [[random.random() for i in [1,2,3]] for j in xrange(0,4)]
    Vertexes = [1,2,3,0]
    Liness = [[1,2],[0,3]]
    LUTNames = ['curv','depth']
    LUTs=[[random.random() for i in xrange(1,5)] for j in [1,2]]
    io_vtk.writeSulci('test.vtk',Points, Vertexes, Lines, LUTs=LUTs,
                      LUTNames=LUTNames)
    """
    Fp = open(VTKFile,'w')
    writeHeader(Fp)
    writePoint(Fp, Points)
    writeVrtx(Fp, Vertexes)
    for i in xrange(0,len(Lines)):
        Lines[i] = str(Lines[i][0]) + " " + str(Lines[i][1]) + "\n"
    writeSeg(Fp, Lines)
    if len(LUTs) > 0:
        for i, LUT in enumerate(LUTs):
            if i == 0:
                wrtVrtxLUT(Fp, LUT, LUTNames[i])
            else:
                wrtVrtxLUT(Fp, LUT, LUTNames[i], AtLUTBegin=False)
    Fp.close()

def surf2VTK(SurfFile):
    Vertex, Face = io_file.readSurf(SurfFile)
    Fp = open(SurfFile + '.vtk', 'w')
    writeHeader(Fp, Title='vtk output from'+SurfFile)
    writePoint(Fp, Vertex)
    writeFace(Fp, Face)
    Fp.close()


# new version, activated 03/15/2011
def fcLst2VTK(VTKFile, SurfaceFile, FundiFile, LUT=[], LUTname=[]):
    """
    Load a face list file and a surface file to map faces onto the surface
    and save the result into VTK format

    This function is called libbasin.getBasin()

    """
    Fp = open(VTKFile,'w')
    Vertex, Face = io_file.readSurf(SurfaceFile)
    FundiList = loadFundiList(FundiFile)
    wrtFcFtr(Fp, Vertex, Face, FundiList)

# commented out Forrest 2011-05-30 13:53
#    if CurvFile!= '':
#        Curvature = fileio.readCurv(CurvFile)
#        wrtVrtxLUT(Fp, Curvature, LUTName = 'curvature')
# End of commented out Forrest 20110-05-30 13:53
# Replaced by:
    if LUT!=[] :
        for i in xrange(0, len(LUT)):
            if i == 0:
                wrtVrtxLUT(Fp, LUT[i], LUTname[i])
            else:
                wrtVrtxLUT(Fp, LUT[i], LUTname[i], AtLUTBegin=False)
# End of Replaced by

    Fp.close()

# new version, activated 03/15/2011
def vrtxLst2VTK(VTKFile, SurfaceFile, FundiFile, LUT=[], LUTname=[]):
    """
    Load a face list file and a surface file to map faces onto the surface
    and save the result into VTK format

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
    Fp = open(VTKFile,'w')
    Vertex, Face = io_file.readSurf(SurfaceFile)
    FundiList = loadFundiList(FundiFile)
    writeVrtxFundi(Fp, Vertex, FundiList)
    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            if i == 0:
                wrtVrtxLUT(Fp, LUT[i], LUTname[i])
            else:
                wrtVrtxLUT(Fp, LUT[i], LUTname[i], AtLUTBegin=False)

    Fp.close()

# new version, activated 03/15/2011
def seg2VTK(VTKFile, SurfaceFile, FundiFile, LUT=[], LUTname=[]):
    """
    Load a fundus curve segment list file and a surface file
    to map curve segments onto the surface and save the result into VTK format

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
    Fp = open(VTKFile,'w')
    Vertex, Face = io_file.readSurf(SurfaceFile)
    FundiList = loadSegFundi(FundiFile)
    writeSegFundi(Fp, Vertex, FundiList)

    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            if i == 0:
                wrtVrtxLUT(Fp, LUT[i], LUTname[i])
            else:
                wrtVrtxLUT(Fp, LUT[i], LUTname[i], AtLUTBegin=False)

    Fp.close()

# End of Level 2 functions ---

# Level 3 functions, called by Level 2 functions ---
def writeSeg(Fp, FundiList):
    """
    Write a line segment with two ending points into VTK format

    Parameters
    ===========

    FundiList: list of strings (NOT integers)
       each element of FundiList is a string containing IDs of two vertexes,
       like "1 3"

    """
    Fp.write("LINES " + str(len(FundiList)) + " " + str(len(FundiList)*3) + "\n")
    [Fp.write("2 " + Vrtx) for Vrtx in FundiList]

def loadFundiList(filename):
    """
    Load the file storing face/vertex IDs, which are fundi faces/vertexes.
    """
    f = open(filename, 'r')
    lines = f.readlines()
    for i in xrange(0,len(lines)):
        lines[i] = int(lines[i][:-1])
    return lines

def wrtFcFtr(Fp, Vertex, Face, FundiList):
    """
    Load IDs of faces (in my output format) and original surface file
    to output a face-from feature in VTK format at STDIO
    """
    writeHeader(Fp, 'patient 290 fundi')
    writePoint(Fp, Vertex)

    Fundi = []
    for i in xrange(0, len(FundiList)):
        Fundi.append(Face[FundiList[i]])

    writeFace(Fp, Fundi)

def writeVrtxFundi(Fp, Vertex, FundiList):
    """
    Load IDs of fundi faces (in my output format) and original surface file
    to output fundi in VTK format
    """
    writeHeader(Fp, 'vtk output by Forrest')
    writePoint(Fp, Vertex)

    writeVrtx(Fp, FundiList)

def writeSegFundi(Fp, Vertex, FundiList):
    """
    Load IDs of fundi curve segments (in my output format)
    and original surface file to output fundi in VTK format
    """

    writeHeader(Fp, 'Created by Mindboggle')
    writePoint(Fp, Vertex)
    writeSeg(Fp, FundiList)

def loadSegFundi(Filename):
    """
    Load the file storing fundi as curve segments
    """

    Fp = open(Filename, 'r')
    lines = Fp.readlines()
    Fp.close()

    Segs = []
    for line in lines:
        Segs.append(line)
        # I do NOT convert strings to integers because
        # later we write strings into VTK files.
        # Also no need to split. So break line is also included.

    return Segs


# functions that use VTK's python binding to load VTK files ---
def load_vtk_map(filename):
    """
    Load a VTK-format map that contains only one SCALAR segment

    Inputs
    =======

    filename : string
        The path/filename of a VTK format file.


    Outputs
    =========
    vertices : list of lists of floats
        Each element is a list of 3-D coordinates of a vertex on a surface mesh

    faces : list of lists of integers
        Each element is a 3-tuple, the IDs of vertexes that form a face
        on a surface mesh

    scalars : list of floats
        Each element is a scalar value corresponding to a vertex

    """
    import vtk
    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()

    Data = Reader.GetOutput()
    vertices = [list(Data.GetPoint(PointID))
                for point_id in xrange(0, Data.GetNumberOfPoints())]

    Polys = Data.GetPolys()
    CellArray = Polys.GetData()
    faces = [[CellArray.GetValue(j) for j in xrange(i*4 + 1, i*4 + 4)]
             for i in xrange(0,  CellArray.GetNumberOfCells()) ]

    PointData = Data.GetPointData()
    print("Loading {} {} scalars in file {}...".
          format(Reader.GetNumberOfScalarsInFile,
                 Reader.GetScalarsNameInFile(0), filename))
    ScalarsArray = PointData.GetArray(Reader.GetScalarsNameInFile(0))
    scalars = [ScalarsArray.GetValue(i) for i in xrange(0, ScalarsArray.GetSize())]

    return vertices, faces, scalars



