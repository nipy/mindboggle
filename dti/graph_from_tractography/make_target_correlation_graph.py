#main_folder = '/data/export/home/denis/DTI/50192/tract_freesurfer_parcellation_2'
#saving_file = '/data/export/home/denis/DTI/50192/tract_freesurfer_parcellation_2/prova.gpi'
# python make_target_correlation_graph.py /data/export/home/denis/DTI/50192/tract_freesurfer_parcellation_2 /data/export/home/denis/DTI/50192/tract_freesurfer_parcellation_2/correlation_graph.gpi

""" Prepare the Graph from the DTI study result.
DTI study performed by using fsl 4.1.4 and the target options.
"""

import os,sys
import graph_lib as gl
import networkx as nx
import numpy as np
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
 # Sulcal pits settings
 #seed_folder = '%sROIs/' %(main_folder)
 #seed_templ  = 'ROI_%s.nii.gz'

 #tract_folder_templ = main_folder + 'tract_from_%s/'
 #tract_result_templ = 'seeds_to_ROI_%s.nii.gz'
 #tract_map_templ = 'fdt_paths.nii.gz'

 #voxel_tract_threshold = 10 
 #correlation_threshold = 0.0001
 #graph_type = 'correlation_graph'

# Freesurfer parcellation settings
seed_folder = '%sROIs/' %(main_folder)
seed_templ  = 'seed_%s_mask.nii.gz'

tract_folder_templ = main_folder + 'tract_from_seed_%s/'
tract_result_templ = 'seeds_to_seed_%s_mask.nii.gz'
tract_map_templ = 'fdt_paths.nii.gz'

voxel_tract_threshold = 10 
correlation_threshold = 0.0001
graph_type = 'correlation_graph'

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
print '    - correlation threshold: %s' %(str(correlation_threshold))
DTI_graph = nx.Graph() # Undirected weighted graph
DTI_graph.correlation_threshold = correlation_threshold
DTI_graph.name = saving_file
DTI_graph.type = graph_type

# Insert node into the graph (with properties)
print ' 3) Node insertion...'
for file_path in seed_file_list:
	file_name = gl.get_file_name(file_path)
	node_name_properties = gl.freesurfer_ROI_names(file_name,seed_templ)
	node_name=node_name_properties[0]
	seed_mask = gl.get_map_from_file(file_path)
	node_pos = gl.get_roi_center(seed_mask)
	DTI_graph.add_node(node_name)
	DTI_graph.node[node_name]['ROI_name'] = node_name_properties[1]
	DTI_graph.node[node_name]['path'] = file_path
	DTI_graph.node[node_name]['position'] = node_pos
	DTI_graph.node[node_name]['way_total_map'] = tract_folder_templ %(node_name) + tract_result_templ %(node_name)
	DTI_graph.node[node_name]['fdt_map'] = tract_folder_templ %(node_name) + tract_map_templ
	print '    - node: %s' %(node_name)

# Read connection between seeds
print ' 4) Edge insertion...'
node_list = DTI_graph.nodes()
for cont1 in range(len(node_list)):
	seed_node = node_list[cont1]
	# Loading first node fdt map
	print '\n    - Loading node '+seed_node+' connection map...'
	normalization_map_file = tract_folder_templ %(seed_node) + tract_result_templ %(seed_node)
	streamlines_map_file   = tract_folder_templ %(seed_node) + tract_map_templ
	normalization_map      = np.float32(gl.get_map_from_file(normalization_map_file))
	streamlines_map        = np.float32(gl.get_map_from_file(streamlines_map_file))
	seed_vett       = streamlines_map.ravel()/normalization_map.sum()
	del(normalization_map)
	del(streamlines_map)
	del(normalization_map_file)
	del(streamlines_map_file)
	for cont2 in range(len(node_list)-cont1-1):
		target_node = node_list[cont1+cont2+1]
		# Loading second node fdt map
		print '    - Loading node '+target_node+' connection map...'
		normalization_map_file = tract_folder_templ %(target_node) + tract_result_templ %(target_node)
		streamlines_map_file   = tract_folder_templ %(target_node) + tract_map_templ
		normalization_map      = np.float32(gl.get_map_from_file(normalization_map_file))
		streamlines_map        = np.float32(gl.get_map_from_file(streamlines_map_file))
		target_vett       = streamlines_map.ravel()/normalization_map.sum()
		del(normalization_map)
		del(streamlines_map)
		del(normalization_map_file)
		del(streamlines_map_file)
		print '      - computing node correlation...'
		correlation = np.dot(seed_vett,target_vett)
		if (correlation>=correlation_threshold):
			print '        adding edge. Correlation: %s\n' %(str(correlation))
			DTI_graph.add_weighted_edges_from([(seed_node,target_node,correlation)])
		else:
			print '        discarding edge. Correlation: %s (threshold %s)\n' %(str(correlation),str(correlation_threshold))
					
# Save the graph
print ' 5) Saving the graph...'
nx.write_gpickle(DTI_graph,saving_file)
os.chmod(saving_file,0655)