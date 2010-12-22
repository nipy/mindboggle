""" Prepare the Graph from the DTI study result.
DTI study performed by using fsl 4.1.4 and the target options.
"""

import os, sys
import numpy as np
import nibabel as nb
import scipy as sc
import graph_lib as gl
import networkx as nx
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------------------
#                                                                              IMPORT PHASE
# -----------------------------------------------------------------------------------------
# Subject folder, containing all necessary folders and files
main_folder = sys.argv[1]
if (main_folder[-1] != os.sep):
	main_folder=main_folder + os.sep

saving_file = sys.argv[2]
saving_folder = os.path.dirname(saving_file)
if (saving_folder[-1] != os.sep):
	saving_folder=saving_folder + os.sep

# -----------------------------------------------------------------------------------------
#                                                                                  SETTINGS
# -----------------------------------------------------------------------------------------
# INPUT SUBFOLDERS
seed_folder = '%sROIs/' %(main_folder)
seed_templ  = 'ROI_%s.nii.gz'

tract_folder_templ = main_folder + 'tract_from_%s/'
tract_result_templ = 'seeds_to_ROI_%s.nii.gz'

voxel_tract_threshold = 10 
connection_p_threshold = 0.0001
graph_name = 'connection_graph'
# -----------------------------------------------------------------------------------------
#                                                                              MAIN PROGRAM
# -----------------------------------------------------------------------------------------
# Input data summary
print 'TARGET GRAPH TOOL'
print ' input folder: %s' %(main_folder)
print ' output file:  %s\n' %(saving_folder)

# Read seed files and prepare a list
print ' 1) Looking for seed files...',
pattern = seed_templ %('*')
seed_file_list = gl.find_files(seed_folder,pattern)
print ' %d found' %(len(seed_file_list))

# Graph initialization
print ' 2) Graph initialization...'
DTI_graph = nx.DiGraph() # Directed weighted graph
DTI_graph.voxel_tract_threshold = voxel_tract_threshold
DTI_graph.connection_p_threshold = connection_p_threshold
DTI_graph.name = saving_file
print '    - minimum tract per voxel: %s' %(str(voxel_tract_threshold))
print '    - connection probability threshold: %s' %(str(connection_p_threshold))

# Insert node into the graph (with properties)
print ' 3) Node insertion...'
for file_path in seed_file_list:
	file_name = gl.get_file_name(file_path)
	node_name = file_name[4:]
	seed_mask = gl.get_map_from_file(file_path)
	node_pos = gl.get_roi_mass_center(seed_mask)
	DTI_graph.add_node(node_name)
	DTI_graph.node[node_name]['path'] = file_path
	DTI_graph.node[node_name]['position'] = node_pos
	print '    - node: %s' %(node_name)
	
# Read connection between seeds
print ' 4) Edge insertion...'
for seed_node in DTI_graph.nodes():
	for target_node in DTI_graph.nodes():
		normalization_map_file = tract_folder_templ %(seed_node) + tract_result_templ %(seed_node)
		normalization_map = gl.get_map_from_file(normalization_map_file)
		if seed_node != target_node:
			# Looking for a connection from the seed_node to the target_node
			print '    - Looking for a connection from '+seed_node+' to '+target_node
			connection_map_file = tract_folder_templ %(seed_node) + tract_result_templ %(target_node)
			print '      - loading connection map from '+connection_map_file
			connection_map = gl.get_map_from_file(connection_map_file)
			
			if connection_map.max()>0:
				print '      - looking for streamlines'
				# There is a connection.
				# Computing connection properties	
				seed_voxel_indices = np.nonzero(connection_map)
				numerator = 0
				denominator = 0
				n_seed_voxel = 0
				for cont in range(len(seed_voxel_indices[0])):
					x = seed_voxel_indices[0][cont]
					y = seed_voxel_indices[1][cont]
					z = seed_voxel_indices[2][cont]
					if (connection_map[x][y][z]>=voxel_tract_threshold):
						print '        '+str(connection_map[x][y][z])+' stramlines found'
						numerator = numerator + connection_map[x][y][z]
						denominator = denominator + normalization_map[x][y][z]
						n_seed_voxel = n_seed_voxel+1
					else:
						print '        only '+str(connection_map[x][y][z])+' stramlines found'
				
				if ((numerator!=0)and(denominator!=0)):
					connection_probability = float(numerator)/float(denominator)
					n_tract = numerator;
					print '        Total streamline number: '+str(n_tract)
					print '        Conncetion probability:  '+str(connection_probability)
				else:
					connection_probability = 0
					n_tract = 0
				
				if (connection_probability>=connection_p_threshold):
					print '      - Adding edge'
					# Add the edges to the graph
					DTI_graph.add_weighted_edges_from([(seed_node,target_node,connection_probability)])
					DTI_graph.edge[seed_node][target_node]['n_seed_voxel']=n_seed_voxel
					DTI_graph.edge[seed_node][target_node]['n_streamlines']=n_tract
				else:
					print '      - Edge candidate discarded'
			else:
				print '      - no streamlines found'
					
# Save the graph
print ' 5) Saving the graph...'
nx.write_gpickle(DTI_graph,saving_file)
os.chmod(saving_file,0666)