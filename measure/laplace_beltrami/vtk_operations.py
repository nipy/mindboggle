""" VTK File Operations

--------------------------------------------------------
Parameters:
----------
	VTK = file location
		main file with nodes, meshes and labels
	fundi_file = file location
		file containing polylines of fundi
	
	return_files: boolean
		Would you like to have access to the lines of the file later?
	nodes: boolean
		Would you like a numpy array of the Nodes?
	meshes: boolean
		Would you like a numpy array of the Meshes?
	labels: boolean
		Would you like a numpy array of the Manual Labels?
	fundi: boolean
		Would you like a numpy array of the Fundi?
	delete: boolean
		Would you like to "delete" some of the manual labels (to test algorithm)?
	keep_fundi: boolean
		Would you like to keep the labels of the fundi?
	keep_attached: boolean
		Would you like to keep the labels of the nodes attached to fundi?
	keep_random: boolean
		Would you like to keep the labels of randomly selected nodes?

Returns:
-------
	main_file, fundi_file: files' lines
	Nodes, Meshes, Manual_Labels, Fundi, Fundi_Labels: numpy arrays
	num_changed: int
	preserved_labels: numpy array corresponding to which labels were preserved	
------------------------------------------------
"""

import numpy as np

# Default files to test performance:
default_VTK ='/Users/eli/Dropbox/PythonFiles/Data/testdatalabels.txt'
default_fundi_file ='/Users/eli/Dropbox/PythonFiles/Data/testdatafundi.txt' 

def convert_vtk_to_matrix(VTK=default_VTK, fundi_file=default_fundi_file, return_files=True,
						  nodes=True, meshes=True, labels=True, fundi=True, delete=True, 
						  keep_fundi=True, keep_attached=True, keep_random=False):
	"""Converts vtk files into (optionally) four matrices: 
	   Nodes, Meshes, Manual_Labels, Fundi_Labels."""
	results = {}
	
	f = open(VTK)
	A = [lines for lines in f]
	if return_files:
		results['main_file'] = A
	
	num_points = int(A[4].split()[1])
	num_meshes = int(A[5+num_points].split()[1])
	print 'Number of nodes in the file:', num_points
	
	if nodes:
		Nodes = np.asarray([coordinates for i in xrange(5, 5+num_points) 
			                for coordinates in map(float, A[i].split())]).reshape(-1,3)
		results['Nodes'] = Nodes
	if meshes:
		Meshes = np.asarray([indices for i in xrange(6+num_points, 6+num_points+num_meshes) 
						     for indices in map(int, A[i].split()[1:4])]).reshape(-1,3)
		results['Meshes'] = Meshes
	if labels:	
		Labels = np.asarray(map(int, A[-num_points:]))	
		Manual_Labels = Labels.copy()
		results['Manual_Labels'] = Manual_Labels
	if fundi:
		g = open(fundi_file)
		B = [lines for lines in g]
		if return_files:
			results['fundi_file'] = B
		
		num_lines = int(B[5+num_points].split()[1])
	
		Fundi = np.asarray([indices for i in xrange(6+num_points, 6+num_points+num_lines) 
						    for indices in map(int, B[i].split()[1:3])]).reshape(-1,2)
		results['Fundi'] = Fundi
	# Delete a portion of the manual labels to test the performance of the algorithm:
	if delete:
		# Array of indices of nodes whose labels are to be preserved
		print 'Deleting the labels from some of the nodes...'
		
		preserved_labels = np.zeros(num_points)
		fundal_nodes = np.zeros(num_points)
	
		# If node is part of a fundus:
		if keep_fundi or keep_attached:
			for i in xrange(num_points):
				if i in Fundi:
					fundal_nodes[i] = 1
			
		# If node is part of triangle with a fundal nodes		
		if keep_attached:
			for triangles in Meshes:
				node0, node1, node2 = triangles[0], triangles[1], triangles[2]
				num_nodes_in_fundi = (node0 in Fundi) + (node1 in Fundi) + (node2 in Fundi)
				if num_nodes_in_fundi > 0:
					preserved_labels[triangles] = 1
	
		preserved_labels[fundal_nodes == 1] = int(keep_fundi)
	
		# Change Labels:
		Labels[preserved_labels == 0] = -1
		num_changed = len(np.nonzero(preserved_labels == 0)[0])
	
		# OPTION 2. Just delete random nodes (keep 10%)
		if keep_random and not keep_fundi and not keep_attached:
			num_changed = 0
			for i in xrange(num_points):
				if np.mod(i,10) != 0:
					num_changed += 1
					Labels[i] = -1
			preserved_labels = np.ones(num_points)
			preserved_labels[Labels == -1] = 0
		
		results['num_changed'] = num_changed
		results['preserved_labels'] = preserved_labels
		results['Fundi_Labels'] = Labels
		
		print "num_changed:", num_changed
		print "percent_changed:", (num_changed+0.0)/num_points*100	
		
	return results