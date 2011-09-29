# Copyright (C) 2011 by Forrest Sheng Bao http://fsbao.net

# This software is a free software licensed in MIT license. 

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import fileio

# Level 1: basic VTK element writing functions ----

def writeHeader(Fp, Title='', Header='# vtk DataFile Version 2.0', FileType='ASCII'):
    '''Write the all non-data information for a VTK-format file
    
    This part matches three things in VTK 4.2 File Formats doc
    
    Part 1: Header
    Part 2: Title (256 characters maximum, terminated with newline \n character)
    Part 3: Data type, either ASCII or BINARY 

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

def writeVrtx(Fp, VrtxList):
    """Print vertexes, the VERTICES section in DATASET POLYDATA section 
    """
    # One possible solution
    Fp.write("VERTICES " + str(len(VrtxList)) + " " + str(len(VrtxList)+1) + "\n" + str(len(VrtxList)) + " ")
    [Fp.write(str(i)+" ") for i in VrtxList]
    
    # The other solution, but requires larger files
#    Fp.write("VERTICES " + str(len(VrtxList)) + " " + str(len(VrtxList)*2) + "\n")
#    [Fp.write("1 " + str(i)+"\n") for i in VrtxList]

# end of Level 1 functions  ---

# Level 2: converting my own fundi storage files into VTK representation ---

def surf2VTK(SurfFile):
    Vertex, Face = fileio.readSurf(SurfFile)
    Fp = open(SurfFile + '.vtk', 'w')
    writeHeader(Fp, Title='vtk output from'+SurfFile)
    writePoint(Fp, Vertex)
    writeFace(Fp, Face)
    Fp.close()
    

def fcLst2VTK(VTKFile, SurfaceFile, FundiFile, CurvFile='', LUT=[], LUTname=[]):  # new version, activated 03/15/2011
    '''Load a face list file and a surface file to map faces onto the surface and save the result into VTK format
    
    This function is called libbasin.getBasin()
         
    '''
    Fp = open(VTKFile,'w')
    Vertex, Face = fileio.readSurf(SurfaceFile)    
    FundiList = loadFundiList(FundiFile)
    wrtFcFtr(Fp, Vertex, Face, FundiList)
    
# commented out Forrest 2011-05-30 13:53
#    if CurvFile!= '':
#        Curvature = fileio.readCurv(CurvFile)      
#        wrtVrtxLUT(Fp, Curvature, LUTName = 'curvature')
# End of commented out Forrest 20110-05-30 13:53
# Replaced by:               
    if LUT!=[] and CurvFile != '':
        for i in xrange(0, len(LUT)):
            wrtVrtxLUT(Fp, LUT[i], LUTName=LUTname[i])
# End of Replaced by
    
    Fp.close()
    
def vrtxLst2VTK(VTKFile, SurfaceFile, FundiFile, LUT=[], LUTname=[]):  # new version, activated 03/15/2011
    '''Load a face list file and a surface file to map faces onto the surface and save the result into VTK format
    
    Parameters
    ============
    
    LUT    : list of LUTs
        LUT[i] is the i-th LUT, e.g., distance transformation values, curvature.convexity values.
        LUT[i] is a list of float numbers
        So more than one LUTs may be written to the final VTK file 
        
    LUTname    : list of strings
        LUTname[i] is the name of the lookup table to be inserted into VTK file
     
    '''
    Fp = open(VTKFile,'w')
    Vertex, Face = fileio.readSurf(SurfaceFile)
    FundiList = loadFundiList(FundiFile)
    writeVrtxFundi(Fp, Vertex, FundiList)
    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            wrtVrtxLUT(Fp, LUT[i], LUTName=LUTname[i])
    
    Fp.close()

def seg2VTK(VTKFile, SurfaceFile, FundiFile, LUT=[], LUTname=[]):  # new version, activated 03/15/2011
    '''Load a fundus curve segment list file and a surface file to map curve segments onto the surface and save the result into VTK format

    Parameters
    ============
    
    LUT    : list of LUTs
        LUT[i] is the i-th LUT, e.g., distance transformation values, curvature.convexity values.
        LUT[i] is a list of float numbers
        So more than one LUTs may be written to the final VTK file 
        
    LUTname    : list of strings
        LUTname[i] is the name of the lookup table to be inserted into VTK file
         
    '''
    Fp = open(VTKFile,'w')
    Vertex, Face = fileio.readSurf(SurfaceFile)
    FundiList = loadSegFundi(FundiFile)
    writeSegFundi(Fp, Vertex, FundiList)
    
    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            wrtVrtxLUT(Fp, LUT[i], LUTName=LUTname[i])
    
    Fp.close()

# End of Level 2 functions ---

# Level 3 functions, called by Level 2 functions --- 
def writeSeg(Fp, FundiList):
    '''Write a line segment with two ending points into VTK format 
    
    Parameters
    ===========
    
    FundiList: list of strings (NOT integers)
       each element of FundiList is a string containing IDs of two vertexes, like "1 3" 
    
    '''
    Fp.write("LINES " + str(len(FundiList)) + " " + str(len(FundiList)*3) + "\n")
    [Fp.write("2 " + Vrtx) for Vrtx in FundiList]
        
def loadFundiList(filename):
    '''Load the file storing face/vertex IDs, which are fundi faces/vertexes. 
    '''
    f = open(filename, 'r')
    lines = f.readlines()
    for i in xrange(0,len(lines)):
        lines[i] = int(lines[i][:-1])
    return lines
    
def wrtFcFtr(Fp, Vertex, Face, FundiList):
    '''Load IDs of faces (in my output format) and original surface file to output a face-from feature in VTK format at STDIO 
    '''
    writeHeader(Fp, 'patient 290 fundi')
    writePoint(Fp, Vertex)
    
    Fundi = []
    for i in xrange(0, len(FundiList)):
        Fundi.append(Face[FundiList[i]])
    
    writeFace(Fp, Fundi)

def writeVrtxFundi(Fp, Vertex, FundiList):
    '''Load IDs of fundi faces (in my output format) and original surface file to output fundi in VTK format
    '''
    writeHeader(Fp, 'vtk output by Forrest')
    writePoint(Fp, Vertex)
       
    writeVrtx(Fp, FundiList)

def wrtVrtxLUT(Fp, LUT, LUTName=''):
    '''write per-VERTEX values as a scalar LUT into a VTK file
    
    This function is called by fcLst2VTK
    
    Parameters
    ==========
    
    LUT    : list of floats
    
    ''' 
    Fp.write('POINT_DATA ' + str(len(LUT)) +'\n')
    Fp.write('SCALARS '+ LUTName +' float\n')
    Fp.write('LOOKUP_TABLE '+ LUTName +'\n')
    for Value in LUT:
        Fp.write(str(Value) + '\n')    
    Fp.write('\n')
        

def writeSegFundi(Fp, Vertex, FundiList):
    '''Load IDs of fundi curve segments (in my output format) and original surface file to output fundi in VTK format 
    '''
    
    writeHeader(Fp, 'Created by Mindboggle')
    writePoint(Fp, Vertex)
    writeSeg(Fp, FundiList)
    
def loadSegFundi(Filename):
    '''Load the file storing fundi as curve segments 
    '''
    
    Fp = open(Filename, 'r')
    lines = Fp.readlines()
    Fp.close()
    
    Segs = []
    for line in lines:
        Segs.append(line)   # I do NOT convert strings to integers because later we write strings into VTK files. Also no need to split. So break line is also included.
        
    return Segs