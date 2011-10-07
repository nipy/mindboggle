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

# For more details, please visit the wiki page:
# http://code.google.com/p/mindboggle-utils/wiki/readFreeSurfer

# version 0.1 Forrest Sheng Bao http://fsbao.net
# Last update: 2011-09-20 


import os
import struct 

def readSurf(Filename):
    '''Read in a FreeSurfer Triangle Surface in Binary Format. 
    
    Parameters
    ===========
    Filename : string
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
    >>> Vrtx, Face = rfs.readSurf('lh.pial')
    >>> len(Vrtx)
    130412
    >>> len(Face)
    260820
    >>> Vrtx[10]
    [-7.902474880218506, -95.6839370727539, -21.856534957885742]
    >>> Face[10]
    [2, 39, 3]

    
    
    '''
    f = open(Filename, "rb")
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
    
def readCurv(Filename):
    '''Read in a FreeSurfer curvature (per-vertex) files.
    
    Parameters
    ==========
    
    Filename : string
        A binary FreeSurfer curvature (pre-vertex) file
    
    Outputs
    ========
    
    Curvature : list of floats
        Each element is the curvature value of a FreeSurfer mesh vertex. 
        Elements are ordered by orders of vertexes in FreeSurfer surface file. 
    
    Example
    ========
    
    >>> import readFreeSurfer as rfs
    >>> Curv = rfs.readCurv('lh.curv')
    >>> len(Curv)
    130412
    >>> Curv[10]
    -0.37290969491004944
    
    '''
    f = open(Filename, "rb")
    
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