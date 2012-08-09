#!/usr/bin/python
"""
This Python library reads and writes different file types.
In particular, it has functions to read some FreeSurfer files,
including surface, curvature, and convexity files.
The function read_surface reads in surface files,
while the function read_curvature reads in both
curvature (.curv) and convexity (.sulc) files.

Authors:
Forrest Sheng Bao  .  http://fsbao.net
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def read_surface(filename):
    """
    Read in a FreeSurfer triangle surface mesh in binary format.

    Parameters
    ===========
    filename : string
        A binary FreeSurfer Triangle Surface file

    Outputs
    =======
    Vertex : list of 3-tuples of floats
        Each element is a 3-tuple (list) of floats, which are the X, Y and Z coordinates of a vertex, respectively.
        A 3-tuple's index in the list *Vertex* is the ID of a vertex.

    Face : list of 3-tuples of integers
        Each element is a 3-tuple (list) of integers, which are the IDs of 3 vertexes that form one face

    Example
    ========

    >>> import readFreeSurfer as rfs
    >>> Vrtx, Face = rfs.read_surface('lh.pial')
    >>> len(Vrtx)
    130412
    >>> len(Face)
    260820
    >>> Vrtx[10]
    [-7.902474880218506, -95.6839370727539, -21.856534957885742]
    >>> Face[10]
    [2, 39, 3]

    """

    import os
    import struct

    f = open(filename, "rb")
    f.seek(3)  # skip the first 3 Bytes "Magic" number

    s = f.read(50)   # the second field is a string of variable length
    End2 = s.find('\n\n',0)  # end of the second field is a '\n\n'

    f.seek(3+End2+2)  # jump to immediate Byte after the creating information

    s = f.read(8)
    VertexCount, FaceCount = struct.unpack(">ii", s)
    # print("This hemisphere has", VertexCount, "Vertexes and", FaceCount, "Faces")

    Vertex, Face = [], []

    for i in xrange(0, VertexCount):
        s = f.read(8)
        R, A = struct.unpack(">ff", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        A, S = struct.unpack(">ff", s)
        Vertex.append([R,A,S]) # R, A, S are the coordinates of vertexes

    for i in xrange(0, FaceCount):
        s = f.read(8)
        V0, V1 = struct.unpack(">ii", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        V1, V2 = struct.unpack(">ii", s)
        Face.append([V0, V1, V2])

    return Vertex, Face

def read_curvature(filename):
    """
    Read in a FreeSurfer curvature (per-vertex) file.

    Parameters
    ==========

    filename : string
        A binary FreeSurfer curvature (pre-vertex) file

    Outputs
    ========

    Curvature : list of floats
        Each element is the curvature value of a FreeSurfer mesh vertex.
        Elements are ordered by orders of vertexes in FreeSurfer surface file.

    Example
    ========

    >>> import readFreeSurfer as rfs
    >>> Curv = rfs.read_curvature('lh.curv')
    >>> len(Curv)
    130412
    >>> Curv[10]
    -0.37290969491004944

    """

    import struct

    f = open(filename, "rb")

    f.seek(3) # skip the first 3 Bytes "Magic" number

    s = f.read(8)  # get the VertexCount and FaceCount
    VertexCount, FaceCount = struct.unpack(">ii", s)
    # print("# of Vertexes:", VertexCount, ", # of Faces:", FaceCount)

    Curvature = [0.0]

    s = f.read(8)
    ValsPerVertex, Curvature[0] = struct.unpack(">if", s)

    VertexCount -= 1  # because the first curvature value has been loaded

    while VertexCount > 1:
        s = f.read(8)
        VertexVal1, VertexVal2  =  struct.unpack(">ff", s)
        Curvature += [VertexVal1, VertexVal2]
        VertexCount -= 2

    if VertexCount != 0:  # number of vertexes is even (NOT ODD!!!)
        f.seek(-4, os.SEEK_CUR)  # backward 4 Bytes from current position
        s = f.read(8)
        VertexVal1, VertexVal2 = struct.unpack(">ff", s)
        Curvature.append(VertexVal2)

    f.close()

    return Curvature

def write_list(filename, List):
    """
    Write a list in to a file, each line of which is a list element.
    """
    Fp = open(filename,'w')
    for Element in List:
        Fp.write(str(Element) + '\n')
    Fp.close()

def load_vertex_neighbor_list(filename):
    """
    Load neighbor list of vertexes from a file.

    Input
    ======
        filename: string
            the file from which neighbor list will be loaded

    """
    NbrLst = []
    Fp = open(filename, 'r')
    lines = Fp.readlines()
    for line in lines:
        NbrLst.append([int(i) for i in line.split()])
    Fp.close()
    return NbrLst

def load_faces_neighbor_list(filename):
    """
    Load neighbor list of faces from a file.

    Input
    ======
        filename: string
            the file from which neighbor list will be loaded

    """
    NbrLst = []
    Fp = open(filename, 'r')
    lines = Fp.readlines()
    for line in lines:
        six = [int(i) for i in line.split()]
        NbrLst.append([six[0:3], six[3:6]])
    Fp.close()
    return NbrLst

def write_line_segments(filename, line_paths):
    """
    Write a curve (e.g., fundus) as curve segments, each line contains a segment.
    consisting of the path from one (e.g., fundus) vertex to the nearest (e.g., fundus) vertex.
    """
    Fp = open(filename, 'w')
    for line_path in line_paths:
        if len(line_path) > 1:
            [Fp.write(str(line_path[i]) + '\t' + str(line_path[i+1]) + '\n') for i in xrange(0,len(line_path)-1)]

    Fp.close()

def write_distance_transform(filename, DT_map):
    """
    Write distance transform map into file, each line for a connected component.
    """

    Fp = open(filename, 'w')
    for DT_map in DT_map:
        [Fp.write(str(dist) + '\t') for dist in DT_map]
        Fp.write('\n')
    Fp.close()

def write_lists(filename, input_lists):
    """
    Output list of lists, each line in the file contains one element
    (also a list) of the top-level list.

    Parameters
    ==========

    input_lists : List of lists (2-D so far)
        Each element of input_lists is a 2-D list of equal or unequal size

    Notes
    ======

    2-D lists are seperated by a delimiter which is 4 dashes now: \n----\n

    """

    Fp = open(filename, 'w')
    for input_list in input_lists:
        for Row in input_list:
            for Element in Row:
                Fp.write(str(Element) + '\t')
            Fp.write('\n')
        Fp.write('----\n')
    Fp.close()

def read_list_strings(filename):
    """
    Read a list file.
    """

    import re

    Fp = open(filename, 'r')
    lines = Fp.readlines()
    List = []
    for line in lines:
        if len(re.findall(r'\S+', line)):
            List.append(re.findall(r'\S+', line)[0])
    Fp.close()

    return List

def read_list_2strings(filename):
    """
    Read a 2-column list file.
    """

    import re

    Fp = open(filename, 'r')
    lines = Fp.readlines()
    column1 = []
    column2 = []
    for line in lines:
        if len(re.findall(r'\S+', line)):
            column1.append(re.findall(r'\S+', line)[0])
            column2.append(re.findall(r'\S+', line)[1])
    Fp.close()

    return column1, column2

def read_lists(filename):
    """The reverse function of write_lists

    Assume all data are supposed to be integers. Change if floats are needed.

    """
    Fp = open(filename, 'r')
    Lists = [[]]  # initially, there is one list in lists
    while True:
        Line = Fp.readline()
        if len(Line) < 1 :
            Fp.close()
            break
        else:
            if Line == "----\n":
                Lists.append([])
            else:
                Row = [int(i) for i in Line.split()]
                Lists[-1].append(Row)
    Fp.close()

    return Lists[:-1] # because last one is an empty list

def read_float_lists(filename):
    """Read in float type lists

    """
    Fp = open(filename, 'r')
    Lists = [[]]  # initially, there is one list in lists
    while True:
        Line = Fp.readline()
        if len(Line) < 1 :
            Fp.close()
            break
        else:
            if Line == "----\n":
                Lists.append([])
            else:
                Row = [float(i) for i in Line.split()]
                Lists[-1].append(Row)
    Fp.close()

    return Lists[:-1] # because last one is an empty list

def np_loadtxt(filename):
    """
    Load numpy array from text file.
    """

    from numpy import loadtxt

    return loadtxt(filename)

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
        import os
        os.error("Check format of " + var)