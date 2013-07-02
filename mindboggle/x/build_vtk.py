#!/usr/bin/env python
"""
Build simple VTK files to test realignment of label borders to polylines.

Multiple triangular meshes will be split horizontally into labeled vertices,
and the label boundary will be flanked or intersected by polylines
representing label-delimiting features (such as brain cortical fundus curves).

Author:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np
from mindboggle.utils.io_vtk import write_vtk

N = 99  # divisible by 3
build_labels = False
build_lines = True
verbose = 0

#========
# Points
#========
Points = np.zeros((N**2, 3))
X = np.array(range(N))
Points[:, 0] = X.repeat(N)
count = 0
for i1 in range(N):
    for i2 in range(N):
        Points[count, 1] = i2
        count += 1
if verbose:
    print('Points:\n{0}'.format(Points))

#==========
# Vertices
#==========
Vertices = range(len(Points))

#=======
# Faces
#=======
nfaces = 2 * (N-1)**2
Faces = np.zeros((nfaces, 3), dtype=np.int)
count = 0
for row in range(N - 1):
    for col in range(N - 1):
        index = row * N + col
        # Upper right triangles
        Faces[count, :] = [index, index + 1, index + N + 1]
        count += 1
        # Lower left triangles
        Faces[count, :] = [index, index + N, index + N + 1]
        count += 1
if verbose:
    print('Faces:\n{0}'.format(Faces))

#========
# Labels
#========
if build_labels:

    Labels = np.ones(N**2)
    Labels[np.round((N**2)/2)::] = 2
    Labels = [int(x) for x in Labels]
    if verbose:
        print('Labels:\n{0}'.format(Labels))

    # Write VTK file
    vtk_file = 'test_labels.vtk'
    write_vtk(vtk_file, Points, Vertices, [], Faces, [Labels], ['Labels'])

#=======
# Lines
#=======
#-------------------
# Test 1:
#-------------------
if build_lines:

    vert_offset = 1/8
    vert_inset = vert_offset * N +  N * np.round(N / 8)
    offset_pairs = np.array([[1/6, 1/3], [1/3, 1/3]])
    for itest, offset_pair in enumerate(offset_pairs):

        offset1, offset2 = N * np.round(offset_pair * N)

        Lines = np.zeros((2 * (N - 1), 2), dtype=np.int)
        Labels = np.ones(N**2)
        index = np.round(vert_inset)
        for i in range(offset1, N - 1):
            Lines[i, :] = [index, index + 1]
            Labels[index] = 3
            index += 1
            Labels[index] = 3
        index = np.round(N**2 - vert_inset)
        for i in range(N - 1, 2 * (N - 1) - offset2):
            Lines[i, :] = [index, index + 1]
            Labels[index] = 4
            index += 1
            Labels[index] = 4
        Labels = [int(x) for x in Labels]
        if verbose:
            print('Lines:\n{0}'.format(Lines))

        # Write VTK file
        vtk_file = 'test_linespacing{0}.vtk'.format(itest + 1)
        write_vtk(vtk_file, Points, Vertices, [], Faces, [Labels], ['Labels'])
