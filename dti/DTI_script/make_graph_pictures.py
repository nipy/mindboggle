import show_graph as show
import networkx as nx
import graph_lib as gl
import matplotlib.pyplot as plt
import numpy as np
import os,sys

input_file = '/Users/denis/Desktop/figures/ArnoAnalysis_subjects2.txt'

# GLOBAL VARIABLES
comment_character = '#'
separator_character = ' '

# MAIN CODE
# Checking input file
if not(os.path.isfile(input_file)):
	print 'ERROR: %s file not found' %(input_file)
	raise

# Reading the imput file
graph_list={}
sep_pos_list={}
try:
	file_obj = open(input_file,'r')
except:
	print 'ERROR: can not open file %s' %(input_file)
for line in file_obj:
	if line[0] != comment_character:
		# The row contains a subject
		line = line.strip()
		elements = line.split(separator_character)
		if len(elements) != 3:
			# Error in line, the row doesn't contain the correct number of elements
			print 'ERROR: reading file %s' %(input_file)
			print '       can not interpretate row: %s' %(line)
		else:
			graph_list[elements[0]]=elements[1]
			sep_pos_list[elements[0]]=elements[2]
			
file_obj.close()	

# Preparing the figures
for subject in graph_list.keys():
	print 'Subject %s' %(subject)
	sep_pos = int(sep_pos_list[subject])
	graph = nx.read_gpickle(graph_list[subject])
	print ' sep pos: %s' %(str(sep_pos))
	
	graph_right = gl.get_undirected_graph(gl.separate_hemisphere(graph,2,sep_pos))
	graph_left = gl.separate_hemisphere(graph,1,sep_pos)
	print ' left hemisphere: %s nodes - %s edges' %(len(graph_left.nodes()),len(graph_left.edges()))
	print ' right hemisphere: %s nodes - %s edges' %(len(graph_right.nodes()),len(graph_right.edges()))
	
	# Computing the normalization factor
	print 'Computing the normalization factor...'
	pesi = np.array([])
	for e in graph_right.edges(data=True):
		peso = e[2]['weight']
		pesi = np.hstack((pesi,[100*np.array(peso)]))	
	normalization_factor = pesi.max()
	 
	module_list = gl.find_modules(graph_right)
	max_value1 = 0
	pos1 = 0
	max_value2 = 0
	pos2 = 0
	for cont in range(len(module_list)-1):
		module_subgraph = graph_right.subgraph(module_list[cont+1])
		edgeN = len(module_subgraph.edges())
		nodeN = len(module_subgraph.nodes())
		cl_coef = nx.average_clustering(module_subgraph)
		short_path = nx.average_shortest_path_length(module_subgraph)
		
		# Edge number
		#if max_value < edgeN:
		#	max_value = edgeN
		#	pos = cont+1
			
		# Cluster Coefficient
		#if (max_value < cl_coef)and(nodeN>4):
		#	max_value = nx.average_clustering(module_subgraph)
		#	pos = cont+1
		
		# Nodes and Edges
		#if (max_value < (nodeN*edgeN))and(nodeN>4):
		#	max_value = nx.average_clustering(module_subgraph)
		#	pos = cont+1
		
		# Mixed
		mark = (nodeN*cl_coef)/short_path
		
		if (max_value1 < mark)and(nodeN>4):
			temp_mark = max_value1
			temp_pos = pos1
			max_value1 = mark
			pos1 = cont+1
			max_value2 = temp_mark
			pos2 = temp_pos
		elif (max_value2 < mark)and(nodeN>4):
			max_value2 = mark
			pos2 = cont+1
			

	selected_module = pos1
	selected_module2 = pos2
	print ' - selected module: %s ' %(str(selected_module))
	
	
	show.show_graph_v2(graph_right,0,'black','dashed','#A0CBE2',normalization_factor)
	module_graph = graph_right.subgraph(module_list[selected_module])
	show.show_graph_v2(module_graph,0,'red','solid','red',normalization_factor)
	
	module_graph = graph_right.subgraph(module_list[selected_module2])
	show.show_graph_v2(module_graph,0,'green','solid','green',normalization_factor)
	
	
	#show.draw_modules(graph_right,0,True)
	#plt.title(subject+' - right hemisphere - directed')
	#plt.savefig('/Users/denis/Desktop/figures/'+subject+'_1module_right_dir.png')
	#plt.close()
	#print ' - right undir'
	a=raw_input(' -Press a key -')
	plt.close()