# This is a bunch of Python library to read FreeSurfer format files
# including surface file, curvature and convexity fiels.  
# The function readSurf reads in surface files, while readCurv reads in
# both curvature (.curv) and convexity (.sulc) files  

import struct, os

def readSurf(filename):
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
    
def readCurv(filename):
    '''Read FreeSurfer Curvature and Convexity files
    '''
    f = open(filename, "rb")
    
    f.seek(3) # skip the first 3 Bytes "Magic" number
    
    s = f.read(8)  # get the VertexCount and FaceCount
    VertexCount, FaceCount = struct.unpack(">ii", s)
#    print "# of Vertexes:", VertexCount, ", # of Faces:", FaceCount
    
    Curvature = [0.0]
    
    s = f.read(8)
    ValsPerVertex, Curvature[0] = struct.unpack(">if", s)
    
    VertexCount -= 1  # because the first curvature value has been loaded
    
    while VertexCount > 1:
        s = f.read(8)
        VertexVal1, VertexVal2  =  struct.unpack(">ff", s)
        Curvature += [VertexVal1, VertexVal2]
        VertexCount -= 2    
    
    if VertexCount != 0: # number of vertexes is even (NOT ODD!!!)
        f.seek(-4, os.SEEK_CUR)       # backward 4 Bytes from current position
        s = f.read(8)
        VertexVal1, VertexVal2 = struct.unpack(">ff", s)
        Curvature.append(VertexVal2)
            
#    if f.read() == '':
#        print "Loading per-vertex file ", filename, " succeeded. Da-da!"
#        return 0, Curvature
#    elif f.read() != '':
#        print "Loading per-vertex file ", filename, "does not reach the EOF. Oops!"
#        return -1, Curvature
#    else:
#        print "Loading per-vertex file ", filename, " reaches an unknown error. Debug! "
#        return -2, Curvature

    f.close()

    return Curvature

def writeList(File, List):
    '''Write a list in to a file, each line of which is a list element
    '''
    Fp = open(File,'w')
    for Element in List:
        Fp.write(str(Element) + '\n')
    Fp.close()

def loadVrtxNbrLst(Filename):
    '''Load neighbor list of vertexes from a file

    Input
    ======
        Filename: string
            the file from which neighbor list will be loaded

    '''
    NbrLst = []
    Fp = open(Filename, 'r')
    lines = Fp.readlines()
    for line in lines:
        NbrLst.append([int(i) for i in line.split()])
    Fp.close()
    return NbrLst

def loadFcNbrLst(Filename):
    '''Load neighbor list of faces from a file

    Input
    ======
        Filename: string
            the file from which neighbor list will be loaded

    '''
    NbrLst = []
    Fp = open(Filename, 'r')
    lines = Fp.readlines()
    for line in lines:
        six = [int(i) for i in line.split()]
        NbrLst.append([six[0:3], six[3:6]])
    Fp.close()
    return NbrLst



def writeFundiSeg(Filename, Paths):
    '''Write fundi as curve segments, each line contains segments consisting the path from a fundus vertex to the nearest the other fundus vertex. 
    '''
    Fp = open(Filename, 'w')
    for Path in Paths:
        if len(Path)>1:
            [Fp.write(str(Path[i]) + '\t' + str(Path[i+1]) + '\n') for i in xrange(0,len(Path)-1)]
            
    Fp.close()
    
def writeDTMap(Filename, Maps):
    '''Write distance transform map into file, each line for a connected component 
    '''
    
    Fp = open(Filename, 'w')
    for Map in Maps:
        [Fp.write(str(Dist) + '\t') for Dist in Map]
        Fp.write('\n')
    Fp.close()
    
def wrtLists(Filename, Lists):
    '''Output list of lists, each line in the file contains one element (also a list) of the top-level list
    
    Parameters
    ==========
    
    Lists : List of lists (2-D so far)
        Each element of Lists is a 2-D list of equal or unequal size  
    
    Notes 
    ======
    
    2-D lists are seperated by a delimiter which is 4 dashes now, i.e., \n----\n 
    
    '''
    
    Fp = open(Filename, 'w')
    for List in Lists:
        for Row in List:
            for Element in Row:
                Fp.write(str(Element) + '\t')
            Fp.write('\n')
        Fp.write('----\n')
    Fp.close()
    
def readLists(Filename):
    '''The reversing function of wrtLists
    
    Assume all data are supposed to be integers.  --- Change if floats are needed. 
    
    '''
    Fp = open(Filename, 'r')
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

def readFltLsts(Filename):
    '''Read in float type lists 
    
    '''
    Fp = open(Filename, 'r')
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

## Deprecated functions
#def writeCmpnt(Cmpnts, Filename):
#    '''Write connected components into a file, one line for each connected component
#    '''
#    Fp = open(Filename, 'w')
#    for Cmpnt in Cmpnts:
#        for Element in Cmpnt:
#            Fp.write(str(Element) + '\t')
#        Fp.write('\n')
#    Fp.close()
#    
#def loadCmpnt(Filename):
#    '''Read in connected components from a file, one line for each connected component
#    '''
#    Fp = open(Filename, 'r')
#    Cmpnts = []
#    lines = Fp.readlines()
#    Fp.close()
#    
#    for line in lines:
#        Cmpnt = [int(i) for i in line.split()]
#        Cmpnts.append(Cmpnt)
#        
#    return Cmpnts