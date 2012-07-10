""" Python Module for the Spectral Analysis of Shapes:
###################################################

Overview of Functions:
#######################

(A) Initialize Object **

(B) Import Data
  1- add_nodes **
  2- add_mesh **
  3- add_labels **
  4- import_vtk **
  5- import_fundi **
  6- set_id **

(C) Pre-Processing of Data
  1- compute_mesh_measure **
  2- compute_angles **
  3- compute_smallest_angles **
  4- remove_isolated **
  5- refine_mesh **
  6- fix_triangles **
  7- initialize_labels (Done, check!)
  8- check_well_formed **
  9- create_vtk **
  10- pre_process **

(D) Processing of Data
  1- compute_lbo
  2- propagate_labels

(E) Post-Processing of Data

(F) Analysis of Data

(G) Visualization of Data
"""

###########################
# Imports:
import numpy as np
import pyvtk
from time import time
from subprocess import Popen, PIPE, STDOUT
from scipy.sparse import csr_matrix, lil_matrix
import pylab
import os
import atlas_functions
import create_unittest_files

import vtk_operations as vo
import compute_weights as cw
import graph_operations as go
import pickle

# np.set_printoptions(threshold='nan')

###################################
# Base Class:

class Shape:
	"""
	Shape Class.
	1) Import data into object from either a vtk file or manually.
	2) Construct vtk if necessary.
	3) Pre-process data if necessary.
	4) Run LBO code.
	5)
	"""

	# 'Initialize Object' Method

	def __init__(self, id_tag='PyShape Object'):
		"""Initialize attributes of shape object."""
		self.id = str(id_tag)

		self.Nodes = self.Mesh = self.Labels = self.vtk = 0
		self.has_nodes = self.has_mesh = self.has_labels = self.has_vtk = 0
		self.num_nodes = self.num_faces = 0

		# For computing eigenspectrum of shape
		self.eigenvalues = self.eigenvectors = 0

		# For label propagation
		self.assigned_labels = 0
		self.preserved_labels = 0
		self.Fundi = self.fundal_nodes = 0
		self.border = 0

		# For constructing the neighbors matrix
		self.neighbors_constructed = 0

	"""
#########################################
# ---------------------------------------
#     'Import Data' Methods
# ---------------------------------------
#########################################
	"""

	def add_nodes(self, nodes):
		"""Add 3D coordinates of nodes as 2d array."""
		# Check to make sure that Nodes were inputted as a 2D array.
		# Check to make sure that Nodes are of dimension <= 3

		nodes = np.asarray(nodes)
		if nodes.ndim != 2:
			print 'Please Enter Data as a 2D Array.'
		elif nodes.shape[1] > 3:
			print 'Please Provide Data of Dimension <= 3.'
		else:
			self.Nodes = nodes
			self.has_nodes = 1
			self.num_nodes = self.Nodes.shape[0]
		return 0

	def add_mesh(self, mesh):
		"""Add triangular meshing as 2d array"""

		mesh = np.asarray(mesh)
		if mesh.ndim !=2:
			print 'Please Enter Data as a 2D Array.'
		elif mesh.shape[1] < 2 or mesh.shape[1] > 4:
			print 'Please Provide Polylines, Triangles or Tetrahedra.'
		elif not all([len(set(list(i)))==mesh.shape[1] for i in mesh]):
			print 'Some of your faces reference the same node multiple times.'
		elif np.amax(mesh) >=  self.num_nodes:
			print 'Meshing refers to non-existent nodes. Be advised.'
		else:
			self.Mesh = mesh
			self.has_mesh = 1
			self.num_faces = self.Mesh.shape[0]
		return 0

	def add_labels(self, labels):
		"""Add labels to nodes as 1d array."""

		labels = np.asarray(labels)
		if labels.ndim != 1:
			print 'Please enter labels as a 1D array.'
		elif labels.shape[0] != self.num_nodes:
			print 'You have not provided the appropriate number of labels.'
		else:
			self.Labels = np.asarray(labels)
			self.has_labels = 1
			self.assigned_labels = np.array(self.Labels.size)

		self.set_manual_classes = np.sort(np.asarray(list(set(self.Labels))))
		self.num_manual_classes = len(self.set_manual_classes)

		return 0

	def import_vtk(self, fname, check=1):
		"""Import all data from vtk file."""

		if not isinstance(fname, str):
			print 'Please enter the file name as a string.'
		else:
			Data = pyvtk.VtkData(fname)
			self.Nodes = np.asarray(Data.structure.points)
			self.Mesh = np.asarray(Data.structure.polygons)

			if 'lh' in fname:
				if check:
					conf = raw_input("Shall I change the id_tag to its region and number?[y/n]")
					if conf=='y':
						self.id = 'lh' + Data.header
			elif 'rh' in fname:
				if check:
					conf = raw_input("Shall I change the id_tag to its region and number?[y/n]")
					if conf=='y':
						self.id = 'rh' + Data.header

			self.has_nodes = self.has_mesh = 1
			self.num_nodes = self.Nodes.shape[0]
			self.num_faces = self.Mesh.shape[0]

			if Data.point_data.data != []:
				self.Labels = np.asarray(Data.point_data.data[0].scalars)
				self.has_labels = 1
				self.assigned_labels = np.array(self.Labels.size)
				self.set_manual_classes = np.sort(np.asarray(list(set(self.Labels))))
				self.num_manual_classes = len(self.set_manual_classes)

			self.has_vtk = 1
			self.vtk = open(fname, 'r')

		return self.id

	def import_fundi(self, fname):
		"""Import fundus lines from a vtk file."""
		Data = pyvtk.VtkData(fname)

		new_nodes = np.asarray(Data.structure.points)
		if new_nodes.shape != self.Nodes.shape:
			print 'Nodes in the fundus file do not match nodes in the original file!'
		try:
			self.Fundi = np.asarray(Data.structure.lines)
		except:
			print 'The file does not contain polylines. Please import a different file.'
			return

		self.fundal_nodes = np.asarray(list(set(self.Fundi.flatten())))

		if np.amax(self.Fundi) >= self.num_nodes:
			print 'The fundi reference nodes which are not in the file. Please consider.'

		self.fundi_file = fname

		return self.id

	"""
#########################################
# ---------------------------------------
#     'Pre-Processing of Data' Methods
# ---------------------------------------
#########################################
	"""

	def compute_mesh_measure(self, total=False):
		"""Computes the surface area of a shape object, and finds the area of each triangle."""

		# Check that nodes and meshing have been inputted.
		# Check that if shape is composed of polylines, the area is 0.
		# If shape is composed of tetrahedra, respond that method will not currently.

		if not(self.has_nodes and self.has_mesh):
			print 'Please input both the nodes and meshing of the shape.'
			return

		self.measure = np.zeros(self.Mesh.shape[0])
		if self.Mesh.shape[1] == 2:
			print 'The meshing comprises polylines. Length will be outputted.'
			i = 0
			for line in self.Mesh:
				self.measure[i] = np.linalg.norm(self.Nodes[line[0]] - self.Nodes[line[1]])
				i += 1
		elif self.Mesh.shape[1] == 3:
			print 'The meshing comprises triangles. Area will be outputted.'
			i = 0
			for triangle in self.Mesh:
				a = np.linalg.norm(self.Nodes[triangle[0]] - self.Nodes[triangle[1]])
				b = np.linalg.norm(self.Nodes[triangle[1]] - self.Nodes[triangle[2]])
				c = np.linalg.norm(self.Nodes[triangle[2]] - self.Nodes[triangle[0]])
				s = (a+b+c)/2.0

				self.measure[i] = np.sqrt(s*(s-a)*(s-b)*(s-c))
				i += 1
		elif self.Mesh.shape[1] == 4:
			print 'The meshing comprises tetrahedra. Computation currently unavailable.'
			self.measure = 0

		self.surface_area = sum(self.measure)

		if total:
			return self.surface_area
		else:
			return self.measure

	def compute_angles(self):
		"""Computes the angles for each triangle."""
		# Currently only available for triangles.

		if self.Mesh.shape[1] != 3:
			print 'Sorry, this method only works for triangles.'
			return

		if not(self.has_nodes and self.has_mesh):
			print 'You have yet to add nodes and meshing!'

		self.angles = np.zeros((self.num_faces, 3))
		i = 0
		for triangle in self.Mesh:
			a = np.linalg.norm(self.Nodes[triangle[0]] - self.Nodes[triangle[1]])
			b = np.linalg.norm(self.Nodes[triangle[1]] - self.Nodes[triangle[2]])
			c = np.linalg.norm(self.Nodes[triangle[2]] - self.Nodes[triangle[0]])
			self.angles[i,0] = np.arccos((b**2+c**2-a**2)/(2.0*b*c))
			self.angles[i,1] = np.arccos((a**2+c**2-b**2)/(2.0*a*c))
			self.angles[i,2] = np.arccos((a**2+b**2-c**2)/(2.0*a*b))
			i += 1
		return self.angles

	def compute_smallest_angles(self, threshold=0.03):
		"""Find triangles in meshing with smallest angles."""
		# Currently only available for triangles.

		if self.Mesh.shape[1] != 3:
			print 'Sorry, this method only works for triangles.'
			return

		if not(self.has_nodes and self.has_mesh):
			print 'You have yet to add nodes and meshing!'
			return

		self.compute_angles()
		minima = np.amin(self.angles, 1)

		self.smallest_angles = np.arange(self.angles.size)[minima<threshold]
		return self.smallest_angles

	def remove_isolated(self):
		"""Remove any vertices which are not connected to others via the meshing."""


		if not(self.has_nodes and self.has_mesh):
			print 'You have yet to enter the nodes and meshing!'
			return

		verts = set(np.arange(self.num_nodes))
		meshing = set(self.Mesh.ravel())

		self.isolated = list(set.difference(verts, meshing))

		self.Nodes = np.delete(self.Nodes,self.isolated,0)
		self.num_nodes = self.Nodes.shape[0]

		# Update mesh numbering
		self.isolated = sorted(self.isolated, reverse=True)
		for i in self.isolated:
			for j in xrange(i, self.num_nodes+1):
				self.Mesh[self.Mesh==j] = j - 1

		self.num_faces = self.Mesh.shape[0]
		return self.isolated

	def refine_mesh(self, depth = 1, which_fraction=1):
		"""Refine the meshing of the shape object.
		Option to refine only the largest triangles exists.
		Select which_fraction=float.
		Option to refine multiple times, by depth=int."""

		# Check to make sure that the algorithm works for not only triangles!
		# Or if it does, document that.
		# Check that any time you change the number of nodes, you update the num_nodes attr.
		if not(self.has_nodes and self.has_mesh):
			print 'You have yet to enter the nodes and meshing!'
			return

		while(depth > 0):

			if which_fraction != 1:
				Areas = self.compute_mesh_measure()
				sortedIndex = np.argsort(Areas)
			else:
				sortedIndex = np.arange(self.num_faces)

			num_to_refine = which_fraction * self.num_faces
			old_triangles = []
			threshold = 1E-8

			for i in xrange(num_to_refine-1,-1,-1):
				# Find i'th largest triangle:
				ind = int(np.nonzero(sortedIndex==i)[0])

				# Record index of old triangle to delete later:
				old_triangles.append(ind)

				# Get vertices of this triangular mesh:
				v0, v1, v2 = self.Nodes[self.Mesh[ind,0]], self.Nodes[self.Mesh[ind,1]], self.Nodes[self.Mesh[ind,2]]

				# Find midpoints of each edge:
				mid01 = (v0+v1)/2.0
				mid02 = (v0+v2)/2.0
				mid12 = (v1+v2)/2.0

				# Add vertices to the list of nodes:
				#############################################################
				# Check to make sure vertices aren't already in list of nodes:
				dist = self.Nodes - mid01
				duplicates = [np.linalg.norm(dist[j]) for j in xrange(dist.shape[0])]
				minimum, minindex = np.amin(duplicates), np.argmin(duplicates)

				if minimum < threshold:
					# Duplicate! Assign new vertex the number of old vertex
					ind01 = minindex
				else:
					self.Nodes = np.vstack((self.Nodes,mid01))
					ind01 = self.Nodes.shape[0] - 1

				dist = self.Nodes - mid02
				duplicates = [np.linalg.norm(dist[j]) for j in xrange(dist.shape[0])]
				minimum, minindex = np.amin(duplicates), np.argmin(duplicates)

				if minimum < threshold:
					# Duplicate! Assign new vertex the number of old vertex
					ind02 = minindex
				else:
					self.Nodes = np.vstack((self.Nodes,mid02))
					ind02 = self.Nodes.shape[0] - 1

				dist = self.Nodes - mid12
				duplicates = [np.linalg.norm(dist[j]) for j in xrange(dist.shape[0])]
				minimum, minindex = np.amin(duplicates), np.argmin(duplicates)

				if minimum < threshold:
					# Duplicate! Assign new vertex the number of old vertex
					ind12 = minindex
				else:
					self.Nodes = np.vstack((self.Nodes,mid12))
					ind12 = self.Nodes.shape[0] - 1
				#############################################################

				# Add 4 new triangles:
				self.Mesh = np.vstack((self.Mesh,np.array([[self.Mesh[ind,0],ind01,ind02]])))
				self.Mesh = np.vstack((self.Mesh,np.array([[self.Mesh[ind,1],ind01,ind12]])))
				self.Mesh = np.vstack((self.Mesh,np.array([[self.Mesh[ind,2],ind02,ind12]])))
				self.Mesh = np.vstack((self.Mesh,np.array([[ind12,ind01,ind02]])))

			# Delete triangles which were refined:
			for old in sorted(old_triangles, reverse=True):
				self.Mesh = np.delete(self.Mesh,[old*3,old*3+1,old*3+2]).reshape(-1,3)

			self.num_nodes = self.Nodes.shape[0]
			self.num_faces = self.Mesh.shape[0]

			depth -= 1

		return 0

	def fix_triangles(self, method='delete', threshold=0.03):
		"""Handle the ill-shaped triangles of low quality.
		First attempt: delete the triangles."""

		if not(self.has_nodes and self.has_mesh):
			print 'You have yet to enter the nodes and meshing!'
			return

		# Find triangles with angles below the threshold.
		self.low_quality_triangles = self.compute_smallest_angles()

		if method=='delete':
			# Delete those triangles from the meshing.
			bad_triangles = sorted(self.low_quality_triangles, reverse=True)
			for t in bad_triangles:
				self.Mesh = np.delete(self.Mesh, t, 0)
			self.num_faces = self.Mesh.shape[0]

		return sorted(self.low_quality_triangles)

	# Works. Just update doc strings.
	def initialize_labels(self, keep='fundi', fraction=.05):
		"""Initialize a set of labels to serve as the seeds for label propagation.
		Options include: 'border' for nodes connected to fundi.
						 'fundi' for nodes which are part of fundi.
						 'both' for both the fundi and the borders.
						 'label_boundary' for the nodes which comprise the label boundary.
						 'random' for preserving a <fraction> of random nodes."""

		if not self.has_labels:
			print 'The object does not have any labels. Please add them.'
			return

		self.assigned_labels = np.zeros(self.num_nodes) - 1
		self.preserved_labels = np.zeros(self.num_nodes)

		# To preserve the border nodes, find all nodes which are part of a triangle with
		# fundal nodes, and record their index in the array self.preserved_labes
		if keep in ['border','both']:
			print 'Preserving border nodes...'
			for triangles in self.Mesh:
				node0, node1, node2 = triangles[0], triangles[1], triangles[2]
				num_nodes_in_fundi = (node0 in self.fundal_nodes) + (node1 in self.fundal_nodes) + (node2 in self.fundal_nodes)
				if num_nodes_in_fundi > 0:
					self.preserved_labels[triangles] += 1
			self.preserved_labels[self.fundal_nodes] = 0

		# To preserve the fundi nodes, find all nodes which are part of a fundus, and
		# record their index in the array self.preserved_labels
		if keep in ['fundi','both']:
			print 'Preserving fundal nodes...'
			self.preserved_labels[self.fundal_nodes] = 1

		# To preserve the nodes which comprise the label boundary, call find_label_boundary()
		if keep == 'label_boundary':
			print 'Preserving nodes of the label boundary'
			self.find_label_boundary()
			self.preserved_labels[self.label_boundary] = 1

		# To preserve a fraction of random nodes, keep every 1/fraction'th label.
		if keep == 'random':
			print 'Preserving random nodes...'
			if fraction > 1:
				print 'Please enter a fractional number less than or equal to 1.'
				return

			randoms = np.array([np.mod(i, int(1.0/fraction)) for i in xrange(self.num_nodes)])
			self.preserved_labels[randoms==0] = 1

		# Reassign all positive numbers to 1
		self.preserved_labels[self.preserved_labels > 0] = 1

		# Assign the preserved labels to self.assigned_labels.
		self.assigned_labels[self.preserved_labels==1] = self.Labels[self.preserved_labels==1]

		# self.assigned_labels now contains the 'true' labels of a subset of the nodes.
		# main work of the function is complete.
		# now, provide some statistics for what was done.

		self.num_labels_preserved = len(self.preserved_labels[self.preserved_labels==1])
		self.percent_labels_preserved = (self.num_labels_preserved+0.0)/self.num_nodes * 100

		print 'The percentage of preserved labels: {0}'.format(self.percent_labels_preserved)

		return self.assigned_labels

	def get_label_matrix(self):
		"""Constructs an n x C matrix of labels.
		Input: Array of n labels. -1 corresponds to no label.
		Output n x C matrix. Row corresponds to node, column corresponds to class.
		1 in column = membership in that class. -1 = absence. 0 = unlabeled data."""

		# Remove duplicates
		self.set_of_labels = np.sort(np.asarray(list(set(self.assigned_labels))))

		# If all data is labeled, insert -1 at beginning of list for consistency of later methods.
		if -1 not in self.set_of_labels:
			self.set_of_labels = np.insert(self.set_of_labels, 0, -1)

		# Number of classes and nodes
		C = len(self.set_of_labels) - 1
		n = self.Labels.shape[0]

		# Relabel the classes 0 through C, 0 now indicating no class.
		for i in self.set_of_labels[2:]:
			self.assigned_labels[np.nonzero(self.assigned_labels == i)] = np.nonzero(self.set_of_labels == i)
		self.assigned_labels[np.nonzero(self.assigned_labels == 0)] = 1
		self.assigned_labels[np.nonzero(self.assigned_labels == -1)] = 0

		# Create a dictionary mapping new class labels to old class labels:
		self.label_mapping = dict([(i, int(self.set_of_labels[i+1])) for i in xrange(-1,C)])
		print "Label Mapping: ", self.label_mapping

		# Construct n x C Matrix
		self.label_matrix = np.zeros((n, C))

		for i in xrange(n):
			if self.assigned_labels[i] > 0:
				self.label_matrix[i, :] = -1
				self.label_matrix[i, self.assigned_labels[i] - 1] = 1

		self.num_classes = C

		return self.label_matrix, self.label_mapping

	def find_label_boundary(self, realigned_labels = False):
		""" Find the set of nodes which comprise the label boundary.

			Parameters:
			===========
			realigned_labels: boolean (choice of using realigned labels, as opposed to manual labels)

			Returns:
			=======
			self.label_boundary: numpy array (of indices of nodes which comprise the label boundary)
			self.label_boundary_file: string (VTK file name containing highlighted label boundary)
			--OR--
			self.Rlabel_boundar: numpy array (of indices of nodes which comprise the realigned label boundary)
			self.Rlabel_boundary_file: string (VTK file name containing highlighted realigned label boundary)

			Explanation:
			============

			I will define a label boundary as the set of all nodes
			whose neighbors are not all from the same class.
			Thus, any node which is connected to two (or more) nodes from different classes
			will be a boundary node, and will help constitute the label boundary.

			To do so, we simply go triangle by triangle and see if they are all from the same class.
			If yes, do nothing.
			If no, then add the nodes from the majority class (it'll be 2 to 1) to the label boundary.
			First, let us define an empty array to store the label boundaries.
			It will have zeros and ones. We can then define a new array to store all and only those nodes
			which are part of the label boundary.

			Notes
			=====
			Method seems highly redundant and inelegant. Find better way to incorporate realignment protocol into this method.

		"""

		self.label_boundary = self.Rlabel_boundary = np.zeros(self.num_nodes)

		""" Now, let us go triangle by triangle. Set the node's index in label_boundary to 1
		if and only if it is satisfies the definition.
		"""

		for triangle in self.Mesh:
			v0, v1, v2 = triangle[0], triangle[1], triangle[2]
			# If they are not all the same...
			if not realigned_labels:
				same_labels = [self.Labels[v1]==self.Labels[v2], self.Labels[v0]==self.Labels[v2], self.Labels[v0]==self.Labels[v1]]
			else:
				same_labels = [self.RLabels[v1]==self.RLabels[v2], self.RLabels[v0]==self.RLabels[v2], self.RLabels[v0]==self.RLabels[v1]]

			# Did it this way for option to modify it later. Perhaps reconsider 'not all' statement.
			if not all(same_labels):
				# Then label those nodes as part of the boundary.
				if not realigned_labels:
					self.label_boundary[triangle] = 1
				else:
					self.Rlabel_boundary[triangle] = 1

		# We can now output a file to show the boundary.
		if not realigned_labels:
			filename = '/home/eli/Desktop/label_boundaries_'+self.id+'.vtk'
			vo.write_all(filename,self.Nodes,self.Mesh,self.label_boundary)
			self.label_boundary_file = filename
		else:
			filename = '/home/eli/Desktop/realigned_label_boundaries_'+self.id+'.vtk'
			vo.write_all(filename,self.Nodes,self.Mesh,self.Rlabel_boundary)
			self.Rlabel_boundary_file = filename

		# Reformat label_boundary to include only the indices of those nodes in the boundary.
		if not realigned_labels:
			self.label_boundary = np.nonzero(self.label_boundary==1)[0]
		else:
			self.Rlabel_boundary = np.nonzero(self.Rlabel_boundary==1)[0]

		if not realigned_labels:
			return self.label_boundary, self.label_boundary_file
		else:
			return self.Rlabel_boundary, self.Rlabel_boundary_file

	def find_label_boundary_by_class(self):
		""" Divides the label boundary into classes.

			Parameters:
			===========

			Returns:
			========
			self.label_boundary_by_class: dict (key: int class, value: list of nodes)

			Explanation:
			============
			Helper method for find_label_boundary_segment. No real need for visualization.

		"""

		try:
			self.label_boundary
		except AttributeError:
			self.find_label_boundary() # get the initial boundaries, without regard for class.

		self.label_boundary_by_class = {}
		setA = set(self.label_boundary)

		for Class in self.set_manual_classes:
			setB = set(np.nonzero(self.Labels==Class)[0])
			setC = setA.intersection(setB)

			self.label_boundary_by_class[Class] = list(setC)

		return self.label_boundary_by_class

	def find_label_boundary_segments(self, completed=''):
		""" Breaks up the label boundary into segments.

			Parameters:
			===========
			completed: string (name of pickled file containing dict of segment boundaries, or '' for not yet completed)

			Returns:
			========
			self.label_boundary_segments: dict (key: tuple of classes, value: set of nodes)
			self.segment_file: string (pickled file containing the dictionary, for future, and faster, use)
			self.highlighted_segment_file: string (VTK file with boundary segments highlighted according to class)

			Explanation:
			============
			This method will output a dictionary which will store label boundary segments (and subsegments).
			The key will be a tuple consisting of the class it is in, along with the class it is adjacent to.
			The value will be the set of nodes which comprise the segment.

			As this method is currently very slow, I will dump the results into a file for future unpickling.

		"""
		if not completed:
			self.find_label_boundary_by_class()

			# Initialize dictionary for later ease of use (though there's probably a simpler way to do the concatenation).
			self.label_boundary_segments = {}
			for a in self.set_manual_classes:
				for b in self.set_manual_classes:
					self.label_boundary_segments[(a,b)]=[]


			# Populate the dictionary with nodes
			for Class in self.set_manual_classes:
				print Class
				for node in self.label_boundary_by_class[Class]:
					neighbors = self.neighbors(node)
					A = set(self.Labels[neighbors])
					B = set(list([self.Labels[node]]))
					neighbor_classes = set.difference(A,B)
					for c in neighbor_classes:
						self.label_boundary_segments[(Class,c)] += [node]

			# Trim results - delete any dict entry which has no nodes
			for key, value in self.label_boundary_segments.items():
				if value == []:
					self.label_boundary_segments.pop(key)

			# Print Results
			for key in self.label_boundary_segments.keys():
				print "For classes: ", key, "  ", self.label_boundary_segments[key]

			# Dump pickled results into a file
			self.segment_file = open('/home/eli/Desktop/label_boundary_segments_' + self.id + '.pkl', 'wb')
			pickle.dump(self.label_boundary_segments, self.segment_file)
			self.segment_file.close()

		else:
			self.segment_file = open(completed, 'rb')
			self.label_boundary_segments = pickle.load(self.segment_file)

			# Trim results - delete any dict entry which has no nodes
			for key, value in self.label_boundary_segments.items():
				if value == []:
					self.label_boundary_segments.pop(key)

		# Output the results to a VTK file
		self.highlighted_segment_file = '/home/eli/Desktop/highlighted_segments' + self.id + '.vtk'
		color = 1
		colored_segments = np.zeros(self.Labels.shape)
		for value in self.label_boundary_segments.values():
			colored_segments[value] = color
			color += 1
		vo.write_all(self.highlighted_segment_file, self.Nodes,self.Mesh,colored_segments)

		return self.label_boundary_segments, self.segment_file.name, self.highlighted_segment_file

	def find_endpoints(self,segment):
		""" Find the endpoints of a given label boundary segment.

		Parameters
		==========
		segment: list (containing nodes which comprise the label boundary segment)

		Returns
		=======
		endpoint: list (of size 2, containing both endpoints)

		Problems
		========
		1) There needn't be only two endpoints - the boundary might not be contiguous.
		2) There needn't be even two enpoints - the "end" might be triangular: the degree of the nodes will be two.
		Ah.
		"""

		endpoint = [-1, -1]
		print 'Finding endpoints...'

		for node in segment:
			ns = self.neighbors(node)
			boundary_ns = np.intersect1d(ns, segment)
			if len(boundary_ns) == 1:
				if endpoint[0] == -1:
					endpoint[0] = node
				else:
					endpoint[1] = node
					break

		print "Endpoints are: ", endpoint

		if endpoint == [-1,-1]:
			print segment
			raw_input("Check this out.")

		return endpoint

	def find_intersections(self, segment, endpoint):
		""" Find the terminal nodes which are at the intersection of the fundi and manual label boundary.

		Parameters
		==========
		segment: list (nodes comprising the label boundary segment)
		endpoint: list (two endpoints of the label boundary segment)

		Returns
		=======
		intersection: list (two outermost nodes of the label boundary segment which intersect the fundi)
						If there are not two intersection, return -1s in place

		Explanation
		===========

		"""

		print 'Finding intersection of segment with fundi...'
		intersection = [0,0]

		for i in xrange(2):
			pointer = endpoint[i]
			used_nodes = [pointer]
			neighbors = []
			while pointer not in self.fundal_nodes:
				tmp0 = np.intersect1d(self.neighbors(pointer),segment)
				tmp1 = np.setdiff1d(tmp0,used_nodes)
				neighbors = neighbors + list(tmp1)

				if not neighbors:
					pointer = -1
					print 'No neighbors left to explore.'
					break

				pointer = neighbors.pop()
				used_nodes.append(pointer)

			intersection[i] = pointer

		print "The intersection points are: ", intersection

		if np.product(intersection) < 0:
			print segment
			labels = np.zeros(self.Labels.shape)
			labels[segment] = 100
			labels[endpoint] = 200
			vo.write_all('/home/eli/Desktop/debug_intersections.vtk',self.Nodes, self.Mesh, labels)
			raw_input("Check this out...")

		return intersection

	def get_realignment_matrix(self):
		""" Constructs a label matrix for realignment protocol.

		Parameters
		==========
		[self.label_boundary_segments: dict (contains segments with keys of class membership tuples)]

		Returns
		=======
		self.realignment_matrix: np array (n x num_segments matrix of initial labels, used in propagation)
		self.realignment_mapping: dict (key-value pairs of new-old labels)

		Explanation
		===========
		Method takes the dictionary self.label_boundary_segments and converts it into the label matrix self.realignment_matrix
		for label propagation.
		We have separate segments which should each be assigned separate labels.
		The label propagation algorithms all work on some n x C label matrix.
		It starts with a "sparsely populated" matrix, and fills it out with probabilities.

		The mapping of "keys" to new labels will be stored in a dictionary

		Now, we can simply run through the dictionary, and assign to each node which is part of a segment a label.

		"""

		# Step 0. Find num segments
		self.num_segments = len(self.label_boundary_segments)

		# Step 1. Construct n x self.num_segments matrix of labels - concurrently produce label mapping dictionary
		self.realignment_matrix = np.zeros((self.num_nodes,self.num_segments))

		self.realignment_mapping = {}
		class_label = 0
		for key, value in self.label_boundary_segments.items():
			self.realignment_mapping[class_label] = key
			self.realignment_matrix[value,:] = -1
			self.realignment_matrix[value,class_label] = 1
			class_label += 1

		print "Mapping is: ", self.realignment_mapping


		return self.realignment_matrix, self.realignment_mapping

	def analyze_label_fundi_intersections(self, completed=''):
		""" Find fundal nodes which intersect label boundary, and determine whether they belong to the same fundus curve.

		Runs functions: self.find_label_boundary_segments() - if necessary

		Parameters
		==========
		completed: string (name of pickled file containing dict of segment boundaries, or '' for not yet completed)

		Returns
		=======
		self.label_boundary_segments: dict (updated dict, with code at end signifying what type of label propagation to use)
			0: two fundi intersections from same fundus curve - do fill operation, in essence.
			1: all other cases - do label propagation, see other methods for details.
			This list may change as the algorithm for propagation gets more sophisticated.

		"""

		try:
			self.label_boundary_segments
		except AttributeError:
			self.find_label_boundary_segments(completed=completed)

		for key, segment in self.label_boundary_segments.items():
			if segment:
				endpoint = self.find_endpoints(segment)
				intersection = self.find_intersections(segment, endpoint)

				self.same_fundus = False
				if -1 not in intersection and intersection[0] != intersection[1]:
					# It found 2 different intersections. We must therefore now assess whether these two intersections are
					# part of the same fundus.

					pointer = intersection[0] # This is the first intersection point.
					print "First pointer is: ", pointer
					row_avoid = [] # This will be an array of rows in self.Fundi to avoid
					node_avoid = [pointer] # This will be an array of nodes to avoid

					rows = list(np.nonzero([pointer in row for row in self.Fundi])[0]) # This is a list of rows to explore
					print "And the list of rows to explore is: ", rows

					while rows:
						path_to_follow = rows.pop(0)

						print "Following path: ", path_to_follow

						row_avoid.append(path_to_follow)

						tmp = list(self.Fundi[path_to_follow])
						pointer = set.difference(set(tmp), set(node_avoid)).pop()
						node_avoid.append(pointer)
						print "pointer is now: ", pointer

						if pointer == intersection[1]:
							# Bingo!
							print 'Bingo! Both intersections are part of the same fundus!'
							self.same_fundus = True
							break

						rows = rows + list(set.difference(set(np.nonzero([pointer in row for row in self.Fundi])[0]),set(row_avoid)))
						print "Rows is now: ", rows

				self.label_boundary_segments[key] = np.append(segment,int(self.same_fundus))

		return self.label_boundary_segments

	def realign_boundary(self, skip_file = '/home/eli/label_boundary_segments_PyShape Object.pkl'):
		""" Realign the label boundary to coincide with fundi.

			Runs functions: find_label_boundary_segments()

	        Parameters
		    ==========
		    skip_file: string (name of pickled file containing dict of label_segments to save time)

	        Returns
		    =======
		    output_files: vtk file depicting initial (manual) label boundaries
		            vtk file of fundi
		            vtk file of NEW realigned boundaries
		            vtk file highlighting changed regions

        """

		# Step 1. Find Manual Label Boundary Segments, Construct VTK File for Output
		# The VTK file is stored in self.label_boundary_file
		self.find_label_boundary()

		# Step 2. Construct VTK File of Fundi for Output
		# This has already been obtained, as self.fundi_file

		# Step 3. Find Label Boundary Segments
		self.find_label_boundary_segments(completed=skip_file)

		# Step 4. Propagate Labels from Label Boundary Segments, With Desired Conditions In Place
		self.propagate_labels(realign=True)

		# Step 5.

		return self.label_boundary_file, self.fundi_file

	def neighbors(self, node):
		""" This method will accomplish the simple task of taking a node as input and returning
		an np array of the node's neighbors, as defined by self.Mesh.
		"""
		# First check to see if the neighbors matrix was constructed.
		if not self.neighbors_constructed:

			print 'Constructing neighborhood function.'

			""" Then we must construct it now."""
			self.Neighbors = lil_matrix((self.num_nodes, self.num_nodes))

			for row in self.Mesh:
				self.Neighbors[row[0], row[1]] = self.Neighbors[row[1], row[0]] = 1
				self.Neighbors[row[0], row[2]] = self.Neighbors[row[2], row[0]] = 1
				self.Neighbors[row[1], row[2]] = self.Neighbors[row[2], row[1]] = 1

			self.Neighbors = self.Neighbors.tocsr()
			self.neighbors_constructed = 1

		return np.nonzero(self.Neighbors[node])[1]

	def get_consensus_labels(self,surface_file, annot_files):
		""" Call Arno's function to obtain consensus labels for the brain of interest.

		Parameters
		==========
		brain_vtk_file: string (name of vtk file of current brain)
		similar_brains: list (of strings, each string corresponding to file location of comparable brains)

		Returns
		=======
		self.consensus_labels: np array (of indices of nodes which were labeled homogeneously)

		Explanation
		===========
		Consensus labels will be used to prevent realigning the label boundary across consensus labels.

		"""
		# a,b,c,self.consensus_labels = atlas_functions.majority_vote_label(surface_file, annot_files)

		# For now, and to test it out...
		self.consensus_labels = create_unittest_files.__consensus__()

		return self.consensus_labels

	def check_well_formed(self):
		"""Check whether the inputted data is well formed."""
		# Check that number of labels corresponds to number of nodes.
		# Check that numbers in meshing don't exceed number of nodes.

		if not self.has_nodes:
			print 'There are no nodes!'
		if not self.has_mesh:
			print 'There are no faces!'

		if self.has_labels and self.has_nodes:
			if self.Labels.size != self.num_nodes:
				print 'There is a mismatch betweeen the number of labels provided \
						and the number of nodes in the shape.'
				print 'There are {0} nodes and {1} labels. Please fix'.format(self.Nodes.shape[0],self.Labels.size)

		if self.has_nodes and self.has_labels:
			max_mesh_num = np.amax(self.Mesh)
			if max_mesh_num >= self.num_nodes:
				print 'The meshing contains reference to a non-existent node. Please fix.'

		return 0

	def create_vtk(self, fname, label = 'Labels', header='Shape Analysis by PyShape'):
		"""Create vtk file from imported data."""
		if not(self.has_nodes and self.has_mesh):
			print 'You have yet to enter the nodes and meshing!'
			return

		if not self.has_labels:
			self.Labels = None

		self.vtk = vo.write_all(fname, self.Nodes, self.Mesh, self.Labels, label_type=label, msg=header)
		print 'vtk file was successfully created at: ', self.vtk.name
		self.has_vtk = 1

		return self.vtk.name

	def pre_process(self, fname):
		"""Full pre-processing of the shape object for Reuter's Code."""
		self.remove_isolated()
		self.fix_triangles()
		self.check_well_formed()
		self.create_vtk(fname)

	"""
#########################################
# ---------------------------------------
#     'Processing of Data' Methods
# ---------------------------------------
#########################################
	"""

	def compute_lbo(self, num=500, check=0, fname='/home/eli/Neuroscience-Research/Analysis_Hemispheres/Testing.vtk'):
		"""Computation of the LBO using ShapeDNA_Tria software."""
		# Check that everything has been done properly
		# Create vtk file

		if not(self.has_nodes and self.has_mesh):
			print 'You have yet to enter the nodes and meshing!'
			return

		proceed = 'y'
		if check:
			proceed = raw_input('Has the data been pre-processed?[y/n] ')

		if proceed == 'n':
			print 'Then consider running check_well_formed(), remove_isolated(), fix_triangles()'
			return

		if not self.has_vtk:
			print 'Creating a vtk file for visualization and data processing'
			self.vtk = self.create_vtk(fname)

		# Run Reuter's code:
		outfile = fname[:-4]+'_outfile'
		execute = str('./shapeDNA-tria/shapeDNA-tria --mesh ' + self.vtk.name + ' --num ' + str(num) + ' --outfile ' + outfile + ' --ignorelq')
		params = ' --mesh ' + self.vtk.name + ' --num ' + str(num) + ' --outfile /home/eli/Desktop/outfile_' + self.id

		process = Popen(execute, shell=True, stdout = PIPE, stderr = STDOUT)
		# print repr(process.communicate()[0])
		if self.num_nodes < 5000:
			time.sleep(7)
		else:
			time.sleep(16)
		f = open(outfile)

		self.eigenvalues = np.zeros(num)
		add = False
		i = 0

		for line in f:
			if add:
				line = line.strip('{')
				line = line.replace('}','')
				line = line.replace('\n','')
				vals = line.split(';')
				if vals[-1] is '':
					vals.pop()
				try:
					vals = map(float,vals)
				except:
					vals = [-1]
					print 'Could not properly convert line'

				self.eigenvalues[i:i+len(vals)] = vals
				i += len(vals)
			elif 'Eigenvalues' in line:
				add = True

			if i == num:
				break

		return self.eigenvalues

	def propagate_labels(self,method='weighted_average', realign=False,
	                     kernel=cw.rbf_kernel, sigma=10, vis=True, alpha=1, diagonal=0, repeat=1, max_iters=200, tol=.001, eps=1e-7):
		""" Main function to propagate labels.

		Parameters
		==========
		method: string (choice of algorithm)
		realign: boolean (whether label propagation will be used for realigning boundaries)
		kernel: function (used in constructing affinity matrix)
		sigma: float (parameter in gaussian kernel)
		vis: boolean (whether to show progress of algorithm)
		max_iters: int (number of times to repeat the algorithm before breaking)
		tol: float (threshold to assess convergence of the algorithm)
		diagonal: int (value for diagonal entries in weight matrix, optionally change them)
		eps: float (for numerical stability in some algorithms)
		alpha: float (clamping factor)

		Returns
		=======
		self.probabilistic_assignment: np array (n x C matrix of probabilities that node belongs to a given class)

		"""

		# Step 1. Construct Affinity Matrix - compute edge weights:
		self.aff_mat = cw.compute_weights(self.Nodes,self.Mesh,kernel=kernel,sigma=sigma, add_to_graph=False)

		# Step 2. Transform column of labels into L x C Matrix, one column per class
		if not realign:
			self.get_label_matrix()
		else:
			self.get_realignment_matrix()

		# Step 3. Propagate Labels!
		if method == "weighted_average":
			print 'Performing Weighted Average Algorithm! Parameters: max_iters={0}'.format(str(max_iters))
			self.weighted_average(realign, max_iters, tol, vis=vis)
			# Produces self.probabilistic assignment matrix.

#		elif method == "jacobi_iteration":
#			print 'Performing Jacobi Iteration Algorithm! Parameters: max_iters={0}'.format(str(max_iters))
#			prob_matrix = self.jacobi_iteration(alpha, max_iters, tol, eps)

#		elif method == "Label_Spreading":
#			print 'Performing Label Spreading Algorithm! Parameters: alpha={0}, max_iters={1}'.format(str(alpha), str(max_iters))
#			Graph = label_spreading(G, Label_Matrix, aff_mat, label_mapping, data['num_changed'],
#									repeat, alpha, max_iters, tol)
#		elif method == "Label_Propagation":
#			Graph = label_propagation(G, Label_Matrix, aff_mat, label_mapping, data['num_changed'],
#									  repeat, alpha, eps, max_iters, tol)
#		elif method == "Nearest_Neighbor":
#			print 'Performing 1_nearest_neighbor algorithm!'
#			Graph = nearest_neighbor(G, Label_Matrix, label_mapping)

		else:
			print 'That algorithm is not available.'
			return

		return self.probabilistic_assignment

	"""
#########################################
# ---------------------------------------
#     'Post-Processing of Data' Methods
# ---------------------------------------
#########################################
	"""

	def assign_max_prob_label(self):
		""" Assings hard labels to nodes by finding the highest probability label of self.probabilistic_assignment.
			Outputs results to a vtk file for visualization.

		Parameters
		==========

		Returns
		=======
		self.max_prob_label: np array (of size num_nodes which contains the most likely label for each node)
		self.max_prob_file: string (vtk file containing highest probability labels)


		Explanation
		===========
		This method takes self.probabilistic_assignment and determines the most likely labeling of a node.
		It outputs a separate array, max_prob_label, which contains one label for each node.
		The labels are those of the initial numbering.

		"""

		""" First, we must check that the label propagation algorithm has been called. """

		try:
			a = self.probabilistic_assignment[0]
		except:
			print 'First call propagate_labels().'
			return

		""" Now, we go row by row (i.e. one node at a time), and find the column with the maximum value."""

		max_col = np.argmax(self.probabilistic_assignment,axis=1)

		"""	Let us define an array called max_prob_label which will contain the final values."""

		self.max_prob_label = np.zeros(self.num_nodes)

		""" Now, we use the dict label_mapping to convert this matrix back to the original labeling.
		max_col[i] is the artificial label number."""

		for i in self.num_nodes:
			self.max_prob_label[i] = self.label_mapping[max_col[i]]

		""" self.max_prob_label is now complete. Let us now visualize..."""

		self.max_prob_file = '/home/eli/Desktop/max_prob_visualization.vtk'

		vo.write_all(self.max_prob_file, self.Nodes, self.Mesh, self.max_prob_label)

		return self.max_prob_label, self.max_prob_file

	def assign_realigned_labels(self, filename = '/home/eli/Desktop/max_prob_visualization_realignment.vtk'):
		""" Determines which nodes should be reassigned a new label from realignment propagation; visualizes results.

		Parameters
		==========
		filename: string (output vtk file for visualizing realigned boundary)

		Returns
		=======
		self.RLabels: np array (of size num_nodes which contains the most likely realigned label for each node)
		self.RLabels_file: string (vtk file containing highest probability labels)

		Explanation
		===========
		This method takes self.probabilistic_assignment and determines which nodes should be relabeled.
		It outputs a separate array, self.RLabels, which contains one label for each node.
		The labels are those of the initial numbering.
		Care must be taken to only reassign labels when appropriate.

		First try: reassign labels to anything which is in (0,1] to the second value in the tuple in the key of
		the dictionary self.label_boundary_segments.

		"""

		""" First, we must check that the label propagation algorithm has been called. """

		try:
			a = self.probabilistic_assignment[0]
		except:
			print 'First call propagate_labels().'
			return

		self.RLabels = self.Labels.copy()

		# Collect consensus labels...
		self.get_consensus_labels(0,0)

		# go column by column, find those entries which meet criterion, map to tuple, select second entry
		counter = 0
		nodes_to_change = {}

		# SOMEHOW YOU NEED ACCOUNTABILITY HERE - IF TWO DIFFERENT CLASSES ARE TRYING TO RELABEL THE SAME NODES...
		# YOU NEED A DECISION PROCESS WHICH WILL RESOLVE THE DISPUTE!

		for column in self.probabilistic_assignment.T:
			nodes_to_change[counter] = list(np.nonzero(column > .501)[0])
			if set.intersection(set(nodes_to_change[counter]), set(self.consensus_labels)):
				print 'You are trying to cross consensus labels!'
				nodes_to_change = []
			print nodes_to_change
			self.RLabels[nodes_to_change] = self.realignment_mapping[counter][1]
			counter += 1

		self.RLabels_file = filename

		vo.write_all(self.RLabels_file, self.Nodes, self.Mesh, self.RLabels)

		return self.RLabels, self.RLabels_file



	"""
#########################################
# ---------------------------------------
#     'Analysis of Data' Methods
# ---------------------------------------
#########################################
	"""

	def assess_percent_correct(self):
		""" This method compares the results of label propagation to the "ground truth" labels found
		in the original vtk file."""

		self.percent_labeled_correctly = (np.sum(self.max_prob_label == self.Labels) + 0.0) / self.num_nodes
		return self.percent_labeled_correctly

	# This method is devoid of substance.
	def assess_secondary_probabilities(self):
		""" This method finds out which and how many nodes would have been correctly labeled using
		not the maximum probability, but the 'next' best ones.
		"""

		""" The way we will accomplish is to sort the probabilities in self.probabilistic_labels.
		For each node, we will find out which column reports the actual node"""

	def count_assigned_members(self, label):
		""" Method which returns the number of members in a class. Label used is artificial ordering."""
		self.num_assigned_members = sum(map(int,self.label_matrix[:,label]==1))
		return self.num_assigned_members

	def count_real_members(self, label):
		""" Method which returns the number of members in a class, as per the manual labeling.
		Label used is the label name in the vtk file!"""
		self.num_real_members = sum(np.asarray(map(int,self.Labels==label)))
		return self.num_real_members

	def count_current_members(self, label):
		""" This method counts the number of nodes with a given label, after the
		label propagation algorithm has been run.
		"""
		a = self.probabilistic_assignment[:,label] > 0
		b = map(int, a)
		self.num_current_members = sum(b)
		return self.num_current_members

	"""
#########################################
# ---------------------------------------
#     'Visualization of Data' Methods
# ---------------------------------------
#########################################
	"""

	def highlight(self, class_label):
		"""
		This method highlights a set of nodes which belong to the specified class.
		To accomplish this, we find all the nodes which have that label -
		this can be found in self.Labels -
		and then we create a new array where those indices are labeled 1 and all others are labeled -1.
		"""

		indices = np.asarray(map(int,self.Labels==class_label)) * 2 - 1

		""" That's it. Now we just write the file."""

		filename = '/home/eli/Desktop/highlighted_'+self.id+'_'+str(class_label)+'.vtk'

		vo.write_all(filename, self.Nodes, self.Mesh, indices)

	"""
#########################################
# ---------------------------------------
# 			  Helper Methods
# ---------------------------------------
#########################################
	"""

	def weighted_average(self, realign, max_iters, tol, vis=True):
		"""Performs iterative weighted average algorithm to propagate labels to unlabeled nodes.

		Parameters
		==========
		realign: boolean (whether the propagation is for realignment or not)
		max_iters: int (number of iterations of algorithm to perform before breaking)
		tol: float (threshold at which to consider algorithm complete)
		vis: boolean (whether or not to visualize the progress of the algorithm)

		Returns
		=======
		self.probabilistic_assignment: np array (n x C matrix of probabilities that node belongs to a given class)



		Features: Hard label clamps, probabilistic solution. See: Zhu and Ghahramani, 2002.

		"""

		""" The first approach to be considered in the semi-supervised learning case
		is to propagate labels on the graph.
		A simple algorithm of this sort has been propoosed by Zhu and Ghahramani (2002),
		and starts (like all such algorithms) with a set of n nodes,
		l of which are labeled, and u unlabeled.
		The algorithm takes as its input the affinity matrix W (self.aff_mat).
		From the affinity matrix, one may construct the diagonal degree matrix,
		which is a measure of the total weight (or number of edges) which are attached to a node."""

		self.DDM = go.compute_diagonal_degree_matrix(self.aff_mat, inverse=True)

		""" Next, we must initialize a vector to represent the results of the label propagation algorithm.
		It will contain l labels and u 0's.
		This has already been done by the function initialize_labels,
		and is called self.assigned_labels.
		We will just check to make sure this has been accomplished."""

		if isinstance(self.assigned_labels,int):
			print 'Please initialize the labels by calling self.initialize_labels()'
			return

		""" Now, we can actually proceed to perform the iterative algorithm.
		At each timestep, the labels will be updated to reflect the weighted average
		of adjacent nodes. An important caveat of this algorithm
		is that the labeled nodes remain fixed, or clamped.
		They should not be changed, and will need to be reset.
		We accomplish the reset by recalling that self.preserved_labels
		stores the indexes of those nodes whose labels were preserved,
		and self.Labels contains the actual labels.
		The algorithm repeates itself until either convergence or max_iters
		(which will prevent excessive computation time).
		We must also take care to solve the multi-label problem.
		To do so, we employ a one-vs-all framework, where each label is considered independently,
		and set against the rest of the labels.
		More specifically, self.label_matrix is an n x C matrix, where each row represents a node
		and each column represents class membership. We can go column by column and process the algorithm
		iteratively. So, we'd start at the first column and see which nodes get labeled.
		Then we'd move to the next column and label more nodes.
		Because it is possible (likely) that some nodes will not receive any label,
		and also to account for probabilistic labeling, we will assign a probability
		of a node receiving a label. Then we can report these probabilities.
		So, to begin, let us first construct this probabilistic label assignment:
		This matrix will store a 1 for 100% probability, 0 for 0%, and fractional values for the rest.
		We will rename self.label_matrix for this purpose."""

		if not realign:
			self.probabilistic_assignment = self.label_matrix
		else:
			# For realignment, we need to use a different labeling matrix,
			# one which has each segment assigned a different labels
			self.probabilistic_assignment = self.realignment_matrix

		""" We will later change the -1s to 0s.
		As nodes get labeled, we assign a confidence measure to the labeling and store the value
		in this matrix.
		Now, let us go column by column, and run the weighted averaging algorithm.
		For each column, you're studying one class. Therefore, when updating self.probabilistic_assignment,
		you'll be working with one column at a time too.
		If a label gets node, keep the fractional value, do not simply round to 1 to assign membership."""

		i = 0 # record of class number
		for column in self.probabilistic_assignment.T:

			t0 = time()
			print 'Working on class: ', i
			print 'Number of members initially in this class: ', np.nonzero(column==1)[0].size

			# Set up indices and values to be clamped during propagation
			if not realign:
				restore_indices = self.preserved_labels==1
				restore_values = column[restore_indices]
			else:
				restore_indices = np.hstack((self.label_boundary,self.fundal_nodes))
				restore_values = column[restore_indices]
				# Try clamping the other label boundaries at 0...
				restore_values[restore_values==-1] = 0
				#pylab.plot(restore_values)
				#pylab.show()

			Y_hat_now = csr_matrix(column).transpose()
			converged = False
			counter = 0
			while not converged and counter < max_iters:
				""" The option will exist to visualize the proceedings of the algorithm.
				The results of a number of the iterations will be sent to vtk files which can then be visualized.
				For the visualization, we will construct two types of vtk files.
				The first will be the actual (manual) labels, as found in self.Labels,
				with the class of interest highlighted (=1), and the others blanked out (=-1)
				The other vtk files will be the result of the algorithm.
				"""
				if vis and not realign:
					"""First, we'll find out which class/label we're working with, by calling label_mapping.
					We'll then send that class to the method highlight() which will do the actual work of
					creating the vtk."""
					label = self.label_mapping[i]

					if not counter: # No need to do this more than once :-)
						self.highlight(label)

					""" Next, we'll construct vtk files, assuming that the iteration step is one we care about.
					For our purposes, let's see the early iterations in high density, once low density in the middle,
					and the last iteration before convergence or max_iters.
					So, the numbers of interest will be when counter is between 0 and 10, and at 100.
					We'll also see max_iters/2, and max_iters (or convergence).
					"""
					# Actually, just see if counter is a multiple of a number of interest.
					# set_of_interest = np.array([0,101,202,max_iters-2,max_iters-1])

					""" Let's define a file for output.
					We have the nodes and meshing, and we have the labels which are found in column.todense().flatten()."""

					filename = '/home/eli/Desktop/'+self.id+'_'+str(label)+'_'+str(counter)+'.vtk'

					if not np.mod(counter,1000):
						LABELS = np.zeros(self.num_nodes)
						LABELS[:] = Y_hat_now.todense().T.flatten()
						vo.write_all(filename, self.Nodes, self.Mesh, LABELS)

				Y_hat_next = (self.DDM * self.aff_mat * Y_hat_now).todense() # column matrix
				Y_hat_next[restore_indices, 0] = restore_values # reset
				converged = (np.sum(np.abs(Y_hat_now.todense() - Y_hat_next)) < tol) # check convergence
				# print 'Iteration number {0}, convergence = {1}'.format(str(counter),str(np.sum(np.abs(column.todense() - tmp))))
				Y_hat_now = csr_matrix(Y_hat_next)
				counter += 1

			# Print out the number of iterations, so that we get a sense for future runs.
			# It is also an indication of whether the algorithm converged.

			if counter == max_iters:
				print 'The algorithm did not converge.'
			else:
				print 'The algorithm converged in {0} iterations.'.format(str(counter))

			print 'Done in {0} seconds'.format(str(time()-t0))

			self.probabilistic_assignment[:,i] = Y_hat_now.todense().flatten()

			#print 'There were {0} nodes initially preserved in this class'.format(str(self.count_assigned_members(i)))
			#print 'The file actually had {0} nodes in this class'.format(str(self.count_real_members(self.label_mapping[i])))
			#print 'Using only those nodes which crossed the threshold, there are now: '.format(str(self.count_real_members(i)))

			#pylab.plot(self.probabilistic_assignment[:,i])
			#pylab.show()

			i += 1

		""" Before reporting the probabilistic assignment, we change all -1's, which were used
		to indicate non-membership in a class, into 0's, which signify 0 probability that the
		node belongs to that class.
		Note: because the labels were initially numbered -1 and 1, there will be 'probabilities' below 0.
		So, to obtain a sort of probability distribution which preserves order, we will add 1 to each number,
		and then divide by 2. Thus -1 --> 0, 1 --> 1 and everything else keeps its order."""

		""" Uncomment this section when fixed!!"""
		self.probabilistic_assignment += 1
		self.probabilistic_assignment /= 2

		""" self.probabilistic_assignment is now complete."""
		return self.probabilistic_assignment

	### WORK ON THIS ALGORITHM NEXT. CLARIFY SLIGHT AMBIGUITY IN WORDING.
	def jacobi_iteration(self, alpha, max_iters, tol, eps):
		"""Performs label propagation inspired from Jacobi iterative algorithm
		to propagate labels to unlabeled nodes.
		Features: Soft label clamps (alpha), probabilistic solution.
		See: Chapelle, ch. 11 (algorithm 11.2)."""

		"""The next approach to be considered in the semi-supervised learning case
		is to propagate labels on the graph, using a modified version of the above algorithm.
		The main differences are soft clamping, forcing the diagonals to be equal to 0,
		and the introduction of an error term (eps) for numerical stability.

		We start with a set of n nodes, l of which are labeled, and u unlabeled.
		The algorithm takes as its input the affinity matrix W (self.aff_mat).
		From the affinity matrix, we construct the diagonal degree matrix,
		which is a measure of the total weight (or number of edges) which are attached to a node."""

		self.DDM = go.compute_diagonal_degree_matrix(self.aff_mat, inverse=True)

		"""Next, we must initialize a vector to represent the results
		of the label propagation algorithm. It will contain l labels and u 0's.
		This has already been done by the function initialize_labels(),
		and is called self.assigned_labels.
		We will just check to make sure this has been accomplished."""

		if isinstance(self.assigned_labels,int):
			print 'Please initialize the labels by calling self.initialize_labels()'
			return

		""" Now, we can actually proceed to perform the iterative algorithm.
		At each timestep, the labels will be updated to reflect the weighted average
		of adjacent nodes. An important caveat of this algorithm,
		is that the labeled nodes do not remain fixed, or clamped.
		The algorithm repeates itself until either convergence or max_iters
		(which will prevent excessive computation time).
		We must again take care to solve the multi-label problem.
		So, to begin, let us first construct this probabilistic label assignment:
		This matrix will store a 1 for 100% probability, 0 for 0%, and fractional values for the rest.
		We will rename self.label_matrix for this purpose.
		We will later change the -1s to 0s."""

		self.probabilistic_assignment = self.label_matrix

		""" Before proceeding, let us check that the parameters are valid"""

		if not (alpha < 1 and alpha > 0 and eps > 0 and isinstance(max_iters, int) and max_iters > 0 and tol > 0):
			print 'You have failed to properly input parameters. Alpha must be strictly between 0 and 1, eps \ \
					and tol must be greater than 0 and max_iters must be an integer.'
			return

		""" As nodes get labeled, we assign a confidence measure to the labeling and store the value in this matrix.
		Now, let us go column by column, and run the weighted averaging algorithm.
		For each column, you're studying one class. Therefore, when updating self.probabilistic_assignment,
		you'll be working with one column at a time too.
		If a label gets node, keep the fractional value, do not simply round to 1 to assign membership."""

		i = 0 # record of class number

		for column in self.probabilistic_assignment.T:
			self.labeled_indices = column[self.preserved_labels==1]
			column = csr_matrix(column).transpose()
			converged = False
			counter = 0
			while not converged and counter < max_iters:
				tmp = self.DDM * self.aff_mat * column # column matrix
				tmp = tmp.tolil() # store results of iteration
				tmp[self.labeled_indices,0] = self.labeled_indices # reset
				converged = (np.abs(column - tmp).sum() < tol) # check convergence
				print 'convergence=', np.abs(column-tmp).sum()
				column = csr_matrix(tmp)
				counter += 1

			# Print out the number of iterations, so that we get a sense for future runs.
			# It is also an indication of whether the algorithm converged.

			if counter == max_iters:
				print 'The algorithm did not converge.'
			else:
				print 'The algorithm converged in {0} iterations.'.format(str(counter))
			i += 1

		""" Before reporting the probabilistic assignment, we change all -1's, which were used
		to indicate non-membership in a class, into 0's, which signify 0 probability that the
		node belongs to that class."""

		self.probabilistic_assignment[self.probabilistic_assignment==-1] = 0

		return self.probabilistic_assignment

	"""
#########################################
# ---------------------------------------
#           TESTS / DEBUGGING
# ---------------------------------------
#########################################
	"""

shape = Shape()
f1 = '/home/eli/mindboggle/mindboggle/propagate/realignment_test/testdatalabels.vtk'
f2 = '/home/eli/mindboggle/mindboggle/propagate/realignment_test/testdatafundi.vtk'
f3 = '/home/eli/Desktop/label_boundary_segments_PyShape Object.pkl'

g1 = '/home/eli/Desktop/testing1.vtk'
g2 = '/home/eli/Desktop/testing1_fundi.vtk'
g3 = ''

def test():
	""" This test is for the realignment task."""
	t0 = time()
	shape.import_vtk(g1)
	shape.import_fundi(g2)
	print 'Imported Data in: ', time() - t0

	shape.initialize_labels(keep='label_boundary')
	shape.find_label_boundary_segments(completed=g3)
	shape.propagate_labels(realign=True, max_iters=6000)
	shape.assign_realigned_labels()
	shape.find_label_boundary(realigned_labels=True)

test()