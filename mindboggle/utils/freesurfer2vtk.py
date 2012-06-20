#!/usr/bin/python

#Copyright (C) 2011 by Forrest Sheng Bao http://fsbao.net

# This software is licensed under MIT license.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# version 0.2 Last update 2012-06-19

import sys  # this line imports pulic libraries

def writeHeader(Fp, Title='', Header='# vtk DataFile Version 2.0', FileType='ASCII', DataType='POLYDATA'):
    '''Write the all non-data information for a VTK-format file
    
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

    '''
    
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
    """Print coordinates of points, the POINTS section in DATASET POLYDATA section 
    """
    Fp.write("POINTS " + str(len(PointList)) + " " + Type + "\n")
    for i in xrange(0, len(PointList)):
        [R, A, S] = PointList[i]
        Fp.write(str(R) + " " + str(A) + " " + str(S) + "\n")
    
def writeFace(Fp, FaceList, VertexPerFace=3):
    """Print vertexes forming triangular meshes, the POLYGONS section in DATASET POLYDATA section 
    """
    Fp.write("POLYGONS " + str(len(FaceList)) + " " + str( (VertexPerFace + 1) * len(FaceList)  )  + '\n' )
    for i in xrange(0, len(FaceList)):
        [V0, V1, V2] = FaceList[i]
        Fp.write( str(VertexPerFace) + " " + str(V0) + " " + str(V1) + " " + str(V2) + "\n")

def readSurf(filename):
    import struct, os
    f = open(filename, "rb")
    f.seek(3)  # skip the first 3 Bytes "Magic" number
    
    s = f.read(50)   # the second field is string of creation information of variable length
    End2 = s.find('\n\n',0)  # end of the second field is a '\n\n'
    
    f.seek(3+End2+2)  # jump to immediate Byte after the creating information  
    
    s = f.read(8)
    VertexCount, FaceCount = struct.unpack(">ii", s)
#    print "This hemisphere has", VertexCount, "Vertexes and", FaceCount, "Faces"
        
    Vertex, Face = [], []
    
    for i in xrange(0, VertexCount):
        s = f.read(8)
        R, A = struct.unpack(">ff", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        A, S = struct.unpack(">ff", s)
        Vertex.append([R,A,S]) # R, A, S are the coordinates of vertexes
        #print i
        
    for i in xrange(0, FaceCount):
        s = f.read(8)
        V0, V1 = struct.unpack(">ii", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        V1, V2 = struct.unpack(">ii", s)
        Face.append([V0, V1, V2])      
        #print i, V0, V1, V2

    return Vertex, Face

if len(sys.argv) < 1:
    print "Usage: python fsSurf2vtk.py FreeSurfer_Surface VTK_Output"

InputFile = sys.argv[1]
OutputFile = sys.argv[2]

Vrtx, Face = readSurf(InputFile)

OutFP = open(OutputFile, 'w')

writeHeader(OutFP)
writePoint(OutFP, Vrtx)
writeFace(OutFP, Face)

OutFP.close()
