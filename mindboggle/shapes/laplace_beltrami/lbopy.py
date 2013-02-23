""" Python Module for the Spectral Analysis of Shapes:

Authors:
    - Eliezer Stavsky, 2012  (eli.stavsky@gmail.com)

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License
"""

###########################
# Imports:
import numpy as np
import pyvtk
from time import time, sleep
from subprocess import Popen, PIPE, STDOUT

import vtk_operations as vo

###################################
# Base Class:

class lbopy:
	"""Class for the Laplace-Beltrami spectral analysis of shapes.

	Usage:
	=====
	Load vtk file of brain surface: self.load_surface()
	Preprocess data: self.preprocess()
	Compute LBO spectrum: self.compute_lbo()

	Returns:
	=======
	set of eigenvalues: self.eigenvalues.
	"""

	def __init__(self, id_tag='PyShape Object'):
		"""Initialize attributes of shape object."""
		self.id = str(id_tag)

		self.Vertices = self.Faces = self.Labels = self.vtk = 0
		self.has_vertices = self.has_faces = self.has_labels = self.has_vtk = 0
		self.num_vertices = self.num_faces = 0

		# For computing eigenspectrum of shape
		self.eigenvalues = self.eigenvectors = 0

# ---------------------------------------
#     VTK I/O Methods
# ---------------------------------------

	def assign_vertices(self, vertices):
		"""Add 3D coordinates of vertices as 2D array.

		Parameters:
		==========
		vertices: 2D array (of 3 dimensional vertices)
		"""

		vertices = np.asarray(vertices)
		if vertices.ndim != 2:
			print 'Please Enter Data as a 2D Array.'
		elif vertices.shape[1] > 3:
			print 'Please Provide Data of Dimension <= 3.'
		else:
			self.Vertices = vertices
			self.has_vertices = 1
			self.num_vertices = self.Vertices.shape[0]
		return 0

	def assign_faces(self, faces):
		"""Add triangular faces as 2D array.

		Parameters:
		==========
		faces: 2D array (of indices of vertices which comprise the faces)
		"""

		faces = np.asarray(faces)
		if faces.ndim !=2:
			print 'Please Enter Data as a 2D Array.'
		elif faces.shape[1] < 2 or faces.shape[1] > 4:
			print 'Please Provide Polylines, Triangles or Tetrahedra.'
		elif not all([len(set(list(i)))==faces.shape[1] for i in faces]):
			print 'Some of your faces reference the same vertex multiple times.'
		elif np.amax(faces) >=  self.num_vertices:
			print 'Faces reference non-existent vertices. Be advised.'
		else:
			self.Faces = faces
			self.has_faces = 1
			self.num_faces = self.Faces.shape[0]
		return 0

	def assign_labels(self, labels):
		"""Assign labels to the vertices as 1D array.

		Parameters:
		==========
		labels: list or np array (of labels)
		"""

		labels = np.asarray(labels)
		if labels.ndim != 1:
			print 'Please enter labels as a 1D array.'
		elif labels.shape[0] != self.num_vertices:
			print 'You have not provided the appropriate number of labels.'
		else:
			self.Labels = np.asarray(labels)
			self.has_labels = 1

		return 0

	def load_surface(self, fname):
		"""Import vertices, faces and optionally labels from a vtk file.

		Parameters:
		==========
		fname: string (vtk file)
		"""

		if not isinstance(fname, str):
			print 'Please enter the file name as a string.'
		else:
			Data = pyvtk.VtkData(fname)
			self.Vertices = np.asarray(Data.structure.points)
			self.Faces = np.asarray(Data.structure.polygons)

			self.has_vertices = self.has_faces = 1
			self.num_vertices = self.Vertices.shape[0]
			self.num_faces = self.Faces.shape[0]

			if Data.point_data.data != []:
				self.Labels = np.asarray(Data.point_data.data[0].scalars)
				self.has_labels = 1
				self.assigned_labels = np.array(self.Labels.size)
				self.set_manual_classes = np.sort(np.asarray(list(set(self.Labels))))
				self.num_manual_classes = len(self.set_manual_classes)

			self.has_vtk = 1
			self.vtk = open(fname, 'r')

		return self.id

	def check_well_formed(self):
		"""Check whether the inputted data is well formed."""

		if not self.has_vertices:
			print 'There are no vertices!'
		if not self.has_faces:
			print 'There are no faces!'

		if self.has_labels and self.has_vertices:
			if self.Labels.size != self.num_vertices:
				print 'There is a mismatch betweeen the number of labels provided \
						and the number of vertices in the shape.'
				print 'There are {0} vertices and {1} labels. Please fix'.format(self.Vertices.shape[0],self.Labels.size)

		if self.has_vertices and self.has_labels:
			max_faces_num = np.amax(self.Faces)
			if max_faces_num >= self.num_vertices:
				print 'The facesing contains reference to a non-existent vertex. Please fix.'

		return 0

	def write_vtk(self, fname, label = 'Labels', header='Shape Analysis by PyShape'):
		"""Create vtk file from imported data.

		Parameters:
		==========
		fname: string (name of vtk file to be created)
		label: string (label type)
		header: string (vtk file header)

		Returns:
		=======
		self.vtk.name: string (name of created vtk file)
		"""
		if not(self.has_vertices and self.has_faces):
			print 'You have yet to enter the vertices and faces!'
			return

		if not self.has_labels:
			self.Labels = None

		self.vtk = vo.write_all(fname, self.Vertices, self.Faces, self.Labels, label_type=label, msg=header)
		print 'vtk file was successfully created at: ', self.vtk.name
		self.has_vtk = 1

		return self.vtk.name

# ---------------------------------------
#     'Pre-Processing of Data' Methods
# ---------------------------------------

	def compute_faces_measure(self, total=False):
		"""Computes the surface area of a shape object, and finds the area of each triangle.

		Parameters:
		==========
		total: boolean (whether you'd like to return the total surface area)

		Returns:
		=======
		self.surface_area: float (of total surface area)
		--OR--
		self.measure: np array (of area of each triangle in mesh)
		"""

		if not(self.has_vertices and self.has_faces):
			print 'Please input both the vertices and facesing of the shape.'
			return

		self.measure = np.zeros(self.Faces.shape[0])
		if self.Faces.shape[1] == 2:
			print 'The facesing comprises polylines. Length will be outputted.'
			i = 0
			for line in self.Faces:
				self.measure[i] = np.linalg.norm(self.Vertices[line[0]] - self.Vertices[line[1]])
				i += 1
		elif self.Faces.shape[1] == 3:
			print 'The facesing comprises triangles. Area will be outputted.'
			i = 0
			for triangle in self.Faces:
				a = np.linalg.norm(self.Vertices[triangle[0]] - self.Vertices[triangle[1]])
				b = np.linalg.norm(self.Vertices[triangle[1]] - self.Vertices[triangle[2]])
				c = np.linalg.norm(self.Vertices[triangle[2]] - self.Vertices[triangle[0]])
				s = (a+b+c)/2.0

				self.measure[i] = np.sqrt(s*(s-a)*(s-b)*(s-c))
				i += 1
		elif self.Faces.shape[1] == 4:
			print 'The facesing comprises tetrahedra. Computation currently unavailable.'
			self.measure = 0

		self.surface_area = sum(self.measure)

		if total:
			return self.surface_area
		else:
			return self.measure

	def compute_angles(self):
		"""Computes the angles for each triangle.

		Returns:
		=======
		self.angles: np array (num_faces x 3 array of angles size for each triangle)
		"""

		if self.Faces.shape[1] != 3:
			print 'Sorry, this method only works for triangles.'
			return

		if not(self.has_vertices and self.has_faces):
			print 'You have yet to add vertices and facesing!'

		self.angles = np.zeros((self.num_faces, 3))
		i = 0
		for triangle in self.Faces:
			a = np.linalg.norm(self.Vertices[triangle[0]] - self.Vertices[triangle[1]])
			b = np.linalg.norm(self.Vertices[triangle[1]] - self.Vertices[triangle[2]])
			c = np.linalg.norm(self.Vertices[triangle[2]] - self.Vertices[triangle[0]])
			self.angles[i,0] = np.arccos((b**2+c**2-a**2)/(2.0*b*c))
			self.angles[i,1] = np.arccos((a**2+c**2-b**2)/(2.0*a*c))
			self.angles[i,2] = np.arccos((a**2+b**2-c**2)/(2.0*a*b))
			i += 1
		return self.angles

	def compute_smallest_angles(self, threshold):
		"""Find triangles in facesing with smallest angles.

		Parameters:
		==========
		threshold: float (below which a triangle is consider of poor quality)

		Returns:
		=======
		self.smallest_angles: np array (of triangles which fall below threshold)
		"""

		if self.Faces.shape[1] != 3:
			print 'Sorry, this method only works for triangles.'
			return

		if not(self.has_vertices and self.has_faces):
			print 'You have yet to add vertices and faces!'
			return

		self.compute_angles()
		minima = np.amin(self.angles, 1)

		self.smallest_angles = np.arange(self.angles.size)[minima<threshold]
		return self.smallest_angles

	def remove_isolated(self):
		"""Remove any vertices which are not connected to others via the faces.

		Returns:
		=======
		self.isolated: list (of isolated vertices)
		"""

		if not(self.has_vertices and self.has_faces):
			print 'You have yet to enter the vertices and facesing!'
			return

		verts = set(np.arange(self.num_vertices))
		facesing = set(self.Faces.ravel())

		self.isolated = list(set.difference(verts, facesing))

		self.Vertices = np.delete(self.Vertices,self.isolated,0)
		self.num_vertices = self.Vertices.shape[0]

		# Update faces numbering
		self.isolated = sorted(self.isolated, reverse=True)
		for i in self.isolated:
			for j in xrange(i, self.num_vertices+1):
				self.Faces[self.Faces==j] = j - 1

		self.num_faces = self.Faces.shape[0]
		return self.isolated

	def refine_faces(self, depth = 1, which_fraction=1):
		"""Refine the faces of the shape object, with option to refine only the largest triangles.

		Parameters:
		==========
		depth: int (number of times to refine the triangles)
		which_fraction: float (fraction of largest triangles to refine. 1 refines all.)
		"""

		if not(self.has_vertices and self.has_faces):
			print 'You have yet to enter the vertices and faces!'
			return

		while depth:

			if which_fraction != 1:
				Areas = self.compute_faces_measure()
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

				# Get vertices of this triangular faces:
				v0, v1, v2 = self.Vertices[self.Faces[ind,0]], self.Vertices[self.Faces[ind,1]], self.Vertices[self.Faces[ind,2]]

				# Find midpoints of each edge:
				mid01 = (v0+v1)/2.0
				mid02 = (v0+v2)/2.0
				mid12 = (v1+v2)/2.0

				# Add vertices to the list of vertices:
				#############################################################
				# Check to make sure vertices aren't already in list of vertices:
				dist = self.Vertices - mid01
				duplicates = [np.linalg.norm(dist[j]) for j in xrange(dist.shape[0])]
				minimum, minindex = np.amin(duplicates), np.argmin(duplicates)

				if minimum < threshold:
					# Duplicate! Assign new vertex the number of old vertex
					ind01 = minindex
				else:
					self.Vertices = np.vstack((self.Vertices,mid01))
					ind01 = self.Vertices.shape[0] - 1

				dist = self.Vertices - mid02
				duplicates = [np.linalg.norm(dist[j]) for j in xrange(dist.shape[0])]
				minimum, minindex = np.amin(duplicates), np.argmin(duplicates)

				if minimum < threshold:
					# Duplicate! Assign new vertex the number of old vertex
					ind02 = minindex
				else:
					self.Vertices = np.vstack((self.Vertices,mid02))
					ind02 = self.Vertices.shape[0] - 1

				dist = self.Vertices - mid12
				duplicates = [np.linalg.norm(dist[j]) for j in xrange(dist.shape[0])]
				minimum, minindex = np.amin(duplicates), np.argmin(duplicates)

				if minimum < threshold:
					# Duplicate! Assign new vertex the number of old vertex
					ind12 = minindex
				else:
					self.Vertices = np.vstack((self.Vertices,mid12))
					ind12 = self.Vertices.shape[0] - 1
				#############################################################

				# Add 4 new triangles:
				self.Faces = np.vstack((self.Faces,np.array([[self.Faces[ind,0],ind01,ind02]])))
				self.Faces = np.vstack((self.Faces,np.array([[self.Faces[ind,1],ind01,ind12]])))
				self.Faces = np.vstack((self.Faces,np.array([[self.Faces[ind,2],ind02,ind12]])))
				self.Faces = np.vstack((self.Faces,np.array([[ind12,ind01,ind02]])))

			# Delete triangles which were refined:
			for old in sorted(old_triangles, reverse=True):
				self.Faces = np.delete(self.Faces,[old*3,old*3+1,old*3+2]).reshape(-1,3)

			self.num_vertices = self.Vertices.shape[0]
			self.num_faces = self.Faces.shape[0]

			depth -= 1

		return 0

	def fix_triangles(self, method='delete', threshold=0.03):
		"""Handle the ill-shaped triangles of low quality.

		Parameters:
		==========
		method: string (how should the ill-shaped triangles be dealth with?)
		threshold: float (below which a triangle is consider of poor quality)

		Returns:
		=======
		sorted(self.low_quality_triangles): np array (sorted list of bad triangles)
		"""

		if not(self.has_vertices and self.has_faces):
			print 'You have yet to enter the vertices and faces!'
			return

		# Find triangles with angles below the threshold.
		self.low_quality_triangles = self.compute_smallest_angles(threshold=threshold)

		if method=='delete':
			# Delete those triangles from the facesing.
			bad_triangles = sorted(self.low_quality_triangles, reverse=True)
			for t in bad_triangles:
				self.Faces = np.delete(self.Faces, t, 0)
			self.num_faces = self.Faces.shape[0]

		return sorted(self.low_quality_triangles)

	def pre_process(self, fname):
		"""Full pre-processing of the shape object.

		Parameters:
		==========
		fname: string (name of vtk file to be created)
		"""

		self.remove_isolated()
		self.fix_triangles()
		self.check_well_formed()
		self.write_vtk(fname)

# ---------------------------------------
#     'Processing of Data' Methods
# ---------------------------------------

	def compute_lbo(self, fname, num=500):
		"""Computation of the LBO using ShapeDNA_Tria software.

		Parameters:
		==========
		num: int (number of eigenvalues to compute)
		fname: string (vtk file...

		Returns:
		=======
		self.eigenvalues: np array (of first num eigenvalues)
		"""

		if not(self.has_vertices and self.has_faces):
			print 'You have yet to enter the vertices and faces!'
			return

		if not self.has_vtk:
			print 'Creating a vtk file for visualization and data processing'
			self.vtk = self.write_vtk(fname)

		# Run Reuter's code:
		outfile = fname[:-4]+'_outfile'
		execute = str('./shapeDNA-tria/shapeDNA-tria --faces ' + self.vtk.name + ' --num ' + str(num) + ' --outfile ' + outfile + ' --ignorelq')
		params = ' --faces ' + self.vtk.name + ' --num ' + str(num) + ' --outfile /home/eli/Desktop/outfile_' + self.id

		process = Popen(execute, shell=True, stdout = PIPE, stderr = STDOUT)
		# print repr(process.communicate()[0])
		if self.num_vertices < 5000:
			sleep(7)
		else:
			sleep(16)
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

def __main__():

	file_name = '/home/eliezer/mindboggle/mindboggle/propagate/realignment_test/testdatalabels.vtk'
	output_vtk = '/home/eliezer/Desktop/processed_file.vtk'
	test_shape = lbopy()
	test_shape.load_surface(file_name)
	test_shape.pre_process(output_vtk)
	test_shape.compute_lbo(fname=output_vtk)

"""
Log of tasks awaiting completion...

- In fix_triangles(), devise alternate method of handling poor-quality triangles.
- Bypass Reuter's code.
- Incorporate unit testing

"""

"""
Development zone...

- For handling poor-quality triangles:

	The following tests should be done:
		1) Visualize the reformatted vtk file to ensure that the shape is not massively disturbed.
		2) Check to see that the angles get larger.
		3) Check for convergence - i.e. that the problem of poor-quality triangles gets resolved. 

	Problem: determine which vertex should be moved, and the minimum amount it needs to be moved to 
		 increase the size of the angle above threshold. 
		 This should have an analytic solution.

	Set up the problem: We have the indices of three vertices, along with their 3D coordinates. 
		 We can change the coordinates of the vertices we deem necessary. 

	Note: it is possible that the minimum amount of movement will occur by slightly disturbing two different vertices.
	      No renumbering will be necessary.
		
	Solution: Will need to prove this, but it seems that the solution is to move the vertex at the base of the angle
		  along the line which connects it to the midpoint of the other two vertices.

	The algorithm could be made more "subtle" by distributing the labor of vertex movement to all of the vertices.
	You'd need to find the path which each vertex should take, and then have each one take a smaller trip.
	The path that the other two vertices should take is along the vector that is orthogonal to the line connecting
	them to the main vertex. 

	Algorithm: 
		For a given triangle of poor-quality:
			Find which angle fails the test.
			Find which of the three vertices sits at that juncture.
			Find the coordinates of the midpoint of the line which joins the other two vertices.
			Find the vector which connects the vertex of interest to that midpoint. 
			Progressively add small fractions of that vector to the vertex until the angle is less than threshold.
		Check to see that the problem was solved correctly:
			Visualize the results by creating a new vtk file.
			Check that all angles are now below threshold.
			Check that algorithm will converge. 

"""
