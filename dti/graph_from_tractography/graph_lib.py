"""
GRAPH_lib

list of python functions to build graph from connectivity maps
  
  
 ---------------------------------------------------------------------------
 FREESURFER_ROI_NAMES(file_name,seed_templ)
 From the file name retrive the ROI ID number and look for the ROI name.
 INPUT
   * file_name (name of the file containing the ROI)
   * seed_templ (template of the seed file name)
 OUTPUT
   * ROI_properties (ROI ID number and ROI name)

 ---------------------------------------------------------------------------
 GET_UNDIRECTED_GRAPH(graph, method)
 Compute the undirected weighted graph from a directed weighted one. The 
 input parameter <method> contains the method to use to convert.
 INPUT
   * graph (directed, weighted)
   * method (max, average, min, sum)
 OUTPUT
   * graph (undirected, weighted)

 --------------------------------------------------------------------------
 FIND_MODULES(graph)
 Look for modules in a given graph.
 INPUT
   * graph
 OUTPUT
   * Module list. Each element in the list is a list of the nodes in a module

 --------------------------------------------------------------------------
 GET_FILE_NAME(path)
 From a full path return the file name without extension.
 (Used to obtain the node name)
 INPUT
   * path
 OUTPUT
   * file_name

 --------------------------------------------------------------------------
 GET_NODE_NAME(path)
 From a given string return the name for the node
 (Used to obtain the node name)
 INPUT
   * path
 OUTPUT
   * file_name

 -------------------------------------------------------------------------------
 GET_ROI_MASS_CENTER(roi)
 Compute the ROI mass center as the mean voxel position inside the ROI
 INPUT
   * roi
 OUTPUT
   * mass_center
   
 -------------------------------------------------------------------------------
 GET_ROI_CENTER(roi)
 Compute the ROI center as the voxel with the shortest distance from the other 
 voxels. Differently from get_roi_mass_center, in this case the resulting point
 is inside the ROI.
 INPUT
   * roi
 OUTPUT
   * center

 --------------------------------------------------------------------------
 GET_MAP_FROM_FILE(map_file)
 Read the connectivity/probability/mask map from a file.
 INPUT
   * map_file :map file path
 OUTPUT
   * map

 --------------------------------------------------------------------------
 GET_LIST_OF_FILE_WITH_PATTERN(path,pattern)
 Read all files in a given folder and return a list of all files matching a
 given pattern.
 INPUT
   * path
   * pattern
 OUTPUT
   * file_list

 --------------------------------------------------------------------------
 FIND_FILES(path,pattern)
 Check for all files from a folder matching a given pattern
 INPUT
   * string
   * pattern
 OUTPUT
   * boolean value

 --------------------------------------------------------------------------
 ROI_PROPERTY(roi,property)
 Compute a property of a given roi. List of available properties is:
 - maximum value
 - mean value
 INPUT
   * map
   * roi
   * measure_method
 OUTPU
   * roi_property

"""

# ---------------------------------------------------------------------------
#                                                                      IMPORT
# ---------------------------------------------------------------------------

import os, sys
import string
import nibabel
import numpy as np
import nibabel as nb
import networkx as nx
from freesurfer_labels import initialize_ROI_dictionary

# ---------------------------------------------------------------------------
#                                                            GLOBAL VARIABLES
# ---------------------------------------------------------------------------
text_file_separator = ' ' # character used to separate the file names (get_list_from_text_file)
text_file_comment_char = '#' # Character used to comment a row in the text file (get_list_from_text_file)
file_sep = os.sep # Folder separator character (get_normalization_factor)

fsl_prob_log_file = 'probtrackx.log' # Name usually used by fsl for the command log file
fsl_gui_log_file  = 'fdt.log' # Name usually used by fsl for the GUI setup file
fsl_conn_map_file = 'fdt_paths' # Name usually used by fsl for the connectivity map file
fsl_waytotal_file = 'waytotal' # Name usually used by fsl for the connectivity map file

edge_weight_threshold = 0.01 # Minimum weight allowed. If lower the edge is not added to the graph

mk_list_tract_folder = 'tract_%d' + file_sep + fsl_conn_map_file + '.nii.gz'
mk_list_seed_folder = 'seeds' + file_sep + 'seed_%d_mask.nii.gz'

ROI_dictionary = initialize_ROI_dictionary()


# ---------------------------------------------------------------------------
#                                                                     METHODS
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# FREESURFER_ROI_NAMES(file_name,seed_templ)
# From the file name retrive the ROI ID number and look for the ROI name.
# INPUT
#   * file_name (name of the file containing the ROI)
#   * seed_templ (template of the seed file name)
# OUTPUT
#   * ROI_properties (ROI ID number and ROI name)

def freesurfer_ROI_names(file_name,seed_templ):
	seed_templ_file_name = seed_templ.split(".")[0]
	seed_templ_parts = seed_templ_file_name.split("%s")
	ROI_number = file_name
	for cont in range(len(seed_templ_parts)):
		ROI_number = ROI_number.replace(seed_templ_parts[cont],"")
	if ROI_number in ROI_dictionary.keys():
		ROI_properties=[ROI_number, ROI_dictionary[ROI_number]]
	else:
		ROI_properties=[ROI_number, 'Unknown']
	return ROI_properties



# ---------------------------------------------------------------------------
# GET_UNDIRECTED_GRAPH(graph, method)
# Compute the undirected weighted graph from a directed weighted one. The 
# input parameter <method> contains the method to use to convert.
# INPUT
#   * graph (directed, weighted)
#   * method (max, average, min, sum)
# OUTPUT
#   * graph (undirected, weighted)

def get_undirected_graph(graph_obj_IN,method):
	# Create a deep copy of the original graph, but undirected. Then it remove 
	# all edges from the copied graph to add them in a more controlled way.
	graph_obj = graph_obj_IN.to_undirected() 
	graph_obj.remove_edges_from(graph_obj.edges()) 
	for edge in graph_obj_IN.edges(data=True):
		nodo1 = edge[0]
		nodo2 = edge[1]
		if ((nodo1 in graph_obj.nodes())and(nodo2 in graph_obj.nodes())):
			try:
				# Try to look for the inverse edge
				# print 'Edges: '+nodo1+'-'+nodo2
				attributes1 = graph_obj_IN.edge[nodo1][nodo2]
				attributes2 = graph_obj_IN.edge[nodo2][nodo1]
				fields = attributes1.keys()
				# print '  Attributo letto'
				if (method.upper()=='MAX'):
					# print '  metodo max'
					if (attributes1['weight']>=attributes2['weight']):
						final_attributes = attributes1
					else:
						final_attributes = attributes2
					# print '  Used'
				elif (method.upper()=='MEAN'):
					#print '  metodo mean'
					final_attributes = attributes1.copy()
					for field in fields:
						print '  Calcolo la media di ' + field
						final_attributes[field] = 0.5*(attributes1[field]+attributes2[field])
						print '  '+str(attributes1[field])+'-'+str(attributes2[field])+'='+str(final_attributes[field])
					#print '  Used'
				elif (method.upper()=='MIN'):
					#print '  metodo min'
					if (attributes1['weight']>=attributes2['weight']):
						final_attributes = attributes2
					else:
						final_attributes = attributes1
					#print '  Used'
			except:
				print ' Arco inverso non trovato'
				final_attributes = graph_obj_IN.edge[nodo1][nodo2]
				#print ' Uso attributo del primo arco'
		#print ' Aggiungo arco'
		graph_obj.add_weighted_edges_from([(nodo1,nodo2,final_attributes['weight'])])
		for field in final_attributes.keys():
			if field != 'weight':
				graph_obj.edge[nodo1][nodo2][field]=final_attributes[field]
	return graph_obj




# --------------------------------------------------------------------------
# FIND_MODULES(graph)
# Look for modules in a given graph.
# INPUT
#   * graph
# OUTPUT
#   * Module list. Each element in the list is a list of the nodes in a module
def find_modules(graph_obj):
	import community as comm
	
	# Selecting the isolated nodes
	node_degree = nx.degree(graph_obj)
	isolated_nodes = []
	non_isolated_nodes = []
	for n in graph_obj.nodes():
		if node_degree[n] == 0:
			# Isolated node
			isolated_nodes = isolated_nodes + [n]
		else:
			# non-isolated node
			non_isolated_nodes = non_isolated_nodes + [n]
	
	connected_subgraph = graph_obj.subgraph(non_isolated_nodes)
	(Q_score,sub_graph_list) = comm.detect_communities(connected_subgraph)
	
	modules_list = [isolated_nodes] + sub_graph_list 
	return modules_list




# --------------------------------------------------------------------------
# GET_FILE_NAME(path)
# From a full path return the file name without extension.
# (Used to obtain the node name)
# INPUT
#   * path
# OUTPUT
#   * file_name
def get_file_name(full_path):
	full_file_name=os.path.basename(full_path)
	file_name = string.split(full_file_name,'.')[0]
	return file_name




# --------------------------------------------------------------------------
# GET_NODE_NAME(path)
# From a given string return the name for the node
# (Used to obtain the node name)
# INPUT
#   * path
# OUTPUT
#   * file_name
def get_node_name(string):
	file_name = get_file_name(string)
	return file_name




# -------------------------------------------------------------------------------
# GET_ROI_MASS_CENTER(roi)
# Compute the ROI mass center as the mean voxel position inside the ROI
# INPUT
#   * roi
# OUTPUT
#   * mass_center
def get_roi_mass_center(roi):
	indices = roi.nonzero()
	mass_center = []
	for cont in range(len(indices)):
		mass_center = mass_center + [np.round(np.mean(indices[cont]))]
	return mass_center




# -------------------------------------------------------------------------------
# GET_ROI_CENTER(roi)
# Compute the ROI center as the voxel with the shortest distance from the other 
# voxels. Differently from get_roi_mass_center, in this case the resulting point
# is inside the ROI.
# INPUT
#   * roi
# OUTPUT
#   * center

def get_roi_center(roi):
	indices = roi.nonzero()
	nVoxel = len(indices[0])
	distance_vector = np.zeros([nVoxel, 1])
	for cont1 in range(nVoxel):
		temp_vett = np.zeros([nVoxel, 1])
		pos1 = [indices[0][cont1], indices[1][cont1], indices[2][cont1]]
		for cont2 in range(nVoxel):
			pos2 = [indices[0][cont2], indices[1][cont2], indices[2][cont2]]
			temp_vett[cont2] = np.sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]))
		distance_vector = temp_vett.sum()
	pos_min = (temp_vett==temp_vett.min()).nonzero()[0]
	center = [indices[0][pos_min][0],indices[1][pos_min][0],indices[2][pos_min][0]]
	return center



# --------------------------------------------------------------------------
# GET_MAP_FROM_FILE(map_file)
# Read the connectivity/probability/mask map from a file.
# INPUT
#   * map_file :map file path
# OUTPUT
#   * map
def get_map_from_file(map_file):
	map = 0
	if os.path.isfile(map_file):
		# Reading the file

		try:
			map_file_obj = nb.load(map_file)
			map = map_file_obj.get_data()
		except:
			print 'Unexpected error while reading from %s' %(map_file)
	else:
		print 'ERROR: file %s not found' %(map_file)
	return map




# --------------------------------------------------------------------------
# GET_LIST_OF_FILE_WITH_PATTERN(path,pattern)
# Read all files in a given folder and return a list of all files matching a
# given pattern.
# INPUT
#   * path
#   * pattern
# OUTPUT
#   * file_list
def get_list_of_file_with_pattern(path,pattern):
	all_file_list=os.listdir(path)
	matching_file_list=[]
	for file in all_file_list:
		if pattern_matching(file,pattern):
			matching_file_list.append(file)
	return matching_file_list




# --------------------------------------------------------------------------
# FIND_FILES(path,pattern)
# Check for all files from a folder matching a given pattern
# INPUT
#   * string
#   * pattern
# OUTPUT
#   * boolean value

def find_files(path,pattern):
	cmd = 'find %s -name %s -type f' %(path,pattern)
	find_results = os.popen(cmd)
	file_list = []
	for line in find_results.readlines():
		file_list.append(line.rstrip())
	return file_list




# --------------------------------------------------------------------------
# ROI_PROPERTY(roi,property)
# Compute a property of a given roi. List of available properties is:
# - maximum value
# - mean value
# INPUT
#   * map
#   * roi
#   * measure_method
# OUTPU
#   * roi_property

def single_roi_analysis(map,roi,measure_method):
	if map.shape == roi.shape:
		# Dimensions are ok
		roi_indices = roi.nonzero()
		if len(roi_indices[0])>0:
			vett_values = map[roi_indices]
		
			voxel_number = len(vett_values)
			mean_value = vett_values.mean()
			max_value = vett_values.max()
			min_value = vett_values.min()
			std_value = vett_values.std()
		else:
			voxel_number = 0
			mean_value = 0
			max_value = 0
			min_value = 0
			std_value = 0
	else:
		print 'ERROR: map and roi do not have the same dimension'
		raise
		voxel_number = 0
		mean_value = 0
		max_value = 0
		min_value = 0
		std_value = 0
	if measure_method.upper() == 'MEAN':
		return [mean_value]
	elif measure_method.upper() == 'MAX':
		return [max_value]
	elif measure_method.upper() == 'MIN':
		return [min_value]
	elif (measure_method.upper() == 'STD')or(measure_method.upper() == 'SD'):
		return [std_value]
	elif measure_method.upper() == 'NVOXEL':
		return [voxel_number]
	elif measure_method.upper() == 'ALL':
		return [voxel_number, mean_value, std_value, max_value, min_value]
	else:
		print 'ERROR: measurement option unknown (%s)' %(measure_method)
		raise