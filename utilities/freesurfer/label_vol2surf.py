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

# version 0.1 Last update 2011-09-07

# two external librabies needed 
import nibabel
import numpy

# three system libraries needed
import sys
import struct
import os

# two nested functions to read FreeSurfer surface and curvature file
def readSurf(filename):
    f = open(filename, "rb")
    f.seek(3)  # skip the first 3 Bytes "Magic" number
    
    s = f.read(50)   # the second field is string of creation information of variable length
    End2 = s.find('\n\n',0)  # end of the second field is a '\n\n'
    
    f.seek(3+End2+2)  # jump to immediate Byte after the creating information  
    
    s = f.read(8)
    VertexCount, FaceCount = struct.unpack(">ii", s)
        
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
    
def readCurv(filename):
    '''Read FreeSurfer Curvature and Convexity files
    '''
    f = open(filename, "rb")
    
    f.seek(3) # skip the first 3 Bytes "Magic" number
    
    s = f.read(8)  # get the VertexCount and FaceCount
    VertexCount, FaceCount = struct.unpack(">ii", s)
    
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

    f.close()

    return Curvature
# end of two nested functions

# now begins the "real" code 

if len(sys.argv) < 4:
        print "Usage: python project.py volume_file surface_file curvature_file vtk_output"
	print "Example: python project.py LBottom.nii.gz lh.orig lh.curv LBottom.orig.vtk"
	print "For more usage information, check the wiki: http://code.google.com/p/mindboggle-utils/wiki/annot2vtk"
        exit(-1)

# get voxel coordinates
Img =  nibabel.load(sys.argv[1])

Voxel = Img.get_data()

Affine = Img.get_affine()

Idx = numpy.nonzero(Voxel>0)
X, Y, Z = Idx

if len(X) != len(Y) or len(X) != len(Z) or len(Z) != len(Y):
        print "error in nonzero voxel dimension"
        exit(-1)

else:
        Vec = numpy.concatenate((numpy.array(Idx), numpy.ones((1,len(X)))))
        Loc = numpy.dot(Affine, Vec)[0:3]

# end of get voxel coordinates

# begin vertex coordinates and curvature
Vertexes, Faces = readSurf(sys.argv[2])
Curvature = readCurv(sys.argv[3])
# end of vertex coordinates and curvature

# project voxels onto surfaces
Projected = [] # to store the IDs of projected vertexes
for i in xrange(0, len(X)): # for each voxel
	MiniDist, MiniVrtx = 10000, 0 
	for j in xrange(0, len(Vertexes)): # compare the voxel with all vertexes on the surface
		if Curvature[j] > 0:  # only compare with vertexes of positive curvatures
			Vec1 = numpy.array([Loc[0, i], Loc[1,i], Loc[2, i]])
			Vec2 = numpy.array(Vertexes[j])
			Dist = numpy.dot( (Vec1 - Vec2), (Vec1 - Vec2) )
			Dist = numpy.sqrt(Dist)
			if Dist < MiniDist:
				MiniDist = Dist 
				MiniVrtx = j

	Projected.append(MiniVrtx)
# end of project voxels onto surfaces

# output projected voxels 
Fp = open(sys.argv[4], 'w')
Fp.write('# vtk DataFile Version 2.0\n')
Fp.write('created by vol2surf_label_transfer.py of Mindboggle-utils http://code.google.com/p/mindboggle-utils/ \n')
Fp.write('ASCII\n')

Fp.write('DATASET POLYDATA\n')
Fp.write("POINTS " + str(len(Projected)) + " float \n")
for Vrtx in Projected:
	Fp.write( str(Vertexes[Vrtx][0]) + " " +  str(Vertexes[Vrtx][1]) + " " + str(Vertexes[Vrtx][2]) +" \n" )

Fp.write("VERTICES " + str(len(Projected)) + " " + str(len(Projected) + 1) + " \n")
Fp.write(str(len(Projected)) + " ")
for i in xrange(0,len(Projected)):
	Fp.write(str(i) + " ") # i am pretty sure it starts from -0
Fp.write("\n")

Fp.write('POINT_DATA ' + str(len(X)) + ' \n')
Fp.write('SCALARS label integer\n')
Fp.write('LOOKUP_TABLE label \n')
for i in xrange(0,len(X)):
	Fp.write(str(Voxel[X[i], Y[i], Z[i]]) + '\n')

# end of output projected voxels
