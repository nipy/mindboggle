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
labels[.26*num_nodes:.5*num_nodes:num_rows] = 2

vo.write_all(fname, nodes, mesh, labels)

# Create fundus file

fname = '/home/eli/Desktop/testing1_fundi.vtk'

polylines1 = np.zeros((num_rows-1,2))
prev = .25*num_nodes
for i in xrange(num_rows-1):
	next = .25*num_nodes+i+1
	polylines1[i] = [prev,next]
	prev = next

polylines2 = np.zeros((num_rows-1,2))
prev = .25*num_nodes - num_rows
for i in xrange(num_rows-1):
	next = .25*num_nodes - num_rows +i+1
	polylines2[i] = [prev,next]
	prev = next

polylines = np.vstack((polylines1,polylines2))
print polylines.shape

vo.write_all(fname,nodes,polylines,None)

def __consensus__():
	""" Testing consensus labeling scheme..."""

	consensus = np.ones(num_rows)
	consensus = np.cumsum(consensus) + .6*num_nodes
	print consensus
	return consensus

# Blah... Delete when ready...

	#		if self.same_fundus:
	#			self.propagate_labels(realign=True, current_segment= self.label_boundary_segments[(a,b)], fundus_value = 1)
	#
	#		elif intersection[0] == intersection[1] or -1 in intersection:
	#			# It intersects at the same place, or not at all:
	#			# We must now run the label propagation algorithm.
	#			self.propagate_labels(realign=True, current_segment= self.label_boundary_segments[(a,b)], fundus_value = 0)
	#
	#		#################################################################
	#		# Now let us do the co-segment, which is defined by classes (b,a)
	#		#################################################################
	#
	#		if self.label_boundary_segments[(b,a)]:
	#			# partition_points = set.intersection(set(self.label_boundary_segments[(a,b)]),set(self.fundal_nodes))
	#			endpoint = [-1, -1]
	#			for node in self.label_boundary_segments[(b,a)]:
	#				ns = self.neighbors(node)
	#				boundary_ns = set.intersection(set(ns), set(self.label_boundary_segments[(b,a)]))
	#				if len(boundary_ns) == 1:
	#					if endpoint[0] == -1:
	#						endpoint[0] = node # This is an endpoint. Construct segments starting from here.
	#					else:
	#						endpoint[1] = node
	#						break
	#				intersection = [-1,-1] # This array holds the indices of the first two (maximally) intersections of label boundary with fundi
	#			for i in xrange(2):
	#				current_node = endpoint[i]
	#				avoid = {current_node}
	#				while current_node not in self.fundal_nodes:
	#					current_neighbors = set.intersection(set(self.neighbors(current_node)), set(self.label_boundary_segments[(b,a)]))
	#					try:
	#						next_node = set.difference(current_neighbors, avoid).pop()
	#					except:
	#						print 'you have reached the exception clause. The intersection should be -1'
	#						current_node = -1
	#						break
	#					avoid.add(next_node)
	#					current_node = next_node
	#				intersection[i] = current_node # It might be -1, if there are no intersections.
	#				print "The intersection points are: ", intersection
	#				self.same_fundus = False
	#			if -1 not in intersection and intersection[0] != intersection[1]:
	#				# It found 2 different intersections. We must therefore now assess whether these two intersections are
	#				# part of the same fundus.
	#					pointer = intersection[0] # This is the first intersection point.
	#				print "First pointer is: ", pointer
	#				avoid = set() # This will be an array of rows in self.Fundi to avoid
	#				rows = list(np.nonzero([pointer in row for row in self.Fundi])[0]) # This is a list of rows to explore
	#				print "And the list of rows to explore is: ", rows
	#					while rows:
	#					path_to_follow = rows[0]
	#					print "Following path! ", path_to_follow
	#					avoid.add(path_to_follow)
	#					tmp = set(list(self.Fundi[path_to_follow]))
	#					print "fundal nodes which are being analyzed are: ", tmp
	#					pointer = tmp.remove(pointer)
	#					print "pointer is now: ", pointer
	#					if pointer == intersection[1]:
	#						# Bingo!
	#						print 'Bingo! Both intersections are part of the same fundus!'
	#						self.same_fundus = True
	#						break
	#					rows = rows + list(set.difference(set(np.nonzero([pointer in row for row in self.Fundi])[0]),avoid))
	#					print "Rows is now: ", rows
	#				if self.same_fundus:
	#						print 'Propagating labels in the case of YES 2-fundi intersection.'
	#						self.propagate_labels(realign=True, current_segment= self.label_boundary_segments[(a,b)], fundus_value = 1)
	#				elif intersection[0] == intersection[1] or -1 in intersection:
	#					# It intersects at the same place, or not at all:
	#					# We must now run the label propagation algorithm.
	#					print 'Propagating labels in the case of NO 2-fundi intersection.'
	#					self.propagate_labels(realign=True, current_segment= self.label_boundary_segments[(a,b)], fundus_value = 0)
	#			i += 1