__author__ = 'eli'

import vtk_operations as vo
import pyvtk
import numpy as np

# Create vtk data file

fname = '/home/eli/Desktop/testing1.vtk'
num_nodes = 10000
nodes = np.zeros((num_nodes,3))
num_rows = num_cols = int(np.sqrt(nodes.shape[0]))
print num_rows, num_cols

for i in xrange(num_nodes):
	x = i/num_rows
	y = np.mod(i,num_rows)
	nodes[i] = [x,y,0]

mesh = np.zeros((2*(num_rows-1)*(num_cols-1),3))
print mesh.shape
i = 0
for x in xrange(num_rows-1):
	for y in xrange(num_cols-1):
		mesh[i] = [num_rows*x+y,num_rows*x+y+1,num_rows*x+y+num_cols]
		i += 1
		mesh[i] = [num_rows*x+y+1,num_rows*x+y+num_cols,num_rows*x+y+num_cols+1]
		i += 1

labels = np.hstack((np.zeros(nodes.shape[0]/2),np.ones(nodes.shape[0]/2)))

vo.write_all(fname, nodes, mesh, labels)

# Create fundus file

fname = '/home/eli/Desktop/testing1_fundi.vtk'

polylines = np.zeros((num_rows-1,2))
prev = .2*num_nodes
for i in xrange(num_rows-1):
	next = .2*num_nodes+i+1
	polylines[i] = [prev,next]
	prev = next

vo.write_all(fname,nodes,polylines,np.zeros(num_nodes)+2)

def __consensus__():
	""" Testing consensus labeling scheme..."""

	consensus = np.ones(num_rows)
	consensus = np.cumsum(consensus) + .6*num_nodes
	print consensus
	return consensus