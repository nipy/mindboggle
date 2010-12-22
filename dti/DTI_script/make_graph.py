"""
MAKE_GRAPH

INPUT
  * tractography folder
  * save folder 
  * graph file name prefix [optional]

"""
import os, sys
import numpy as np
import nibabel as nb
import scipy as sc
import graph_lib as gl
import networkx as nx
import matplotlib.pyplot as plt

# Variables
file_sep = os.sep
graph_type = 'dir_weighted'
soft_input = 'fsl'
connection_type = 'max'
list_name  = 'seed_list.txt'
graph_name = '%s_graph.gpi' %(graph_type)


# Check for presence of input data
if len(sys.argv)<3:
	print 'ERROR: method %s: 2 input arguments expected, %d found' %(sys.argv[0],len(sys.argv))
	raise
else:
	# Get Input data
	tract_folder = sys.argv[1]
	if tract_folder[-1] != file_sep:
		tract_folder = tract_folder + file_sep
	
	save_folder  = sys.argv[2]
	if save_folder[-1] != file_sep:
		save_folder = save_folder + file_sep
if len(sys.argv)>3:
	temp_input = os.path.splitext(sys.argv[3])
	temp_name = temp_input[0] + '%s.gpi'
	graph_name = temp_name %(graph_type)

print 'METHOD SUMMARY'
print ' - tractography data: %s' %(tract_folder)
print ' - saving folder:     %s' %(save_folder)

# Check Input data
print ' '
print ' Checking data...'
if not(os.path.exists(tract_folder)):
	print 'ERROR: tractography data non found: %s' %(tract_folder)

if not(os.path.exists(save_folder)):
	cmd = 'mkdir -m 777 ' + save_folder
	os.system(cmd)

# Creating of the txt file containing the seed list
print ' Seed file list generation...'
list_fileName = save_folder + list_name
gl.make_list(tract_folder,list_fileName)

# Building the graph
file_list = gl.get_list_from_text_file(list_fileName)
DTI_graph = gl.build_graph(file_list,graph_type,connection_type,soft_input)

# Filter fake edges
for e in DTI_graph.edges():
	if (not(e[0] in DTI_graph.nodes()))or(not(e[1] in DTI_graph.nodes())):
		try:
			DTI_graph.remove_edge(e[0],e[1])
		except:
			print(' fake edge removed')

# Saving the 
nx.write_gpickle(DTI_graph,save_folder + graph_name)