# Show the graph hemisphere or the whole graph
# show_graph(graph,side)
# side: 0 whole, 1 left, 2 right)
def show_graph(graphIN,side):
	import sys,os
	import networkx as nx
	import matplotlib.pyplot as plt
	import numpy as np

	DTI_graph = nx.Graph()
	
	# Input options
	graph_edge_color='black'
	graph_node_color='#A0CBE2'

	# Cropping the graph
	if side==0:
		for n in graphIN.nodes():
			DTI_graph.add_node(n)
			DTI_graph.node[n]['mass center']=graphIN.node[n]['mass center']
	elif side ==1:
		for n in graphIN.nodes():
			pos = graphIN.node[n]['mass center']
			if pos[0]<125:
				DTI_graph.add_node(n)
				DTI_graph.node[n]['mass center']=graphIN.node[n]['mass center']
	elif side ==2:
		for n in graphIN.nodes():
			pos = graphIN.node[n]['mass center']
			if pos[0]>=125:
				DTI_graph.add_node(n)
				DTI_graph.node[n]['mass center']=graphIN.node[n]['mass center']	
	
	# Adding the edges to the cropped graph
	# (Make the graph undirect by selecting the maximum weight between the two edges)
	for e in graphIN.edges(data=True):
		nodo1 = e[0]
		nodo2 = e[1]
		if (nodo1 in DTI_graph.nodes())and(nodo2 in DTI_graph.nodes()):
			try:
				peso2=DTI_graph.edge[nodo2][nodo1]['weight']
				peso1=e[2]['weight']
				if peso1>peso2:
					DTI_graph.edge[nodo2][nodo1]['weight']=peso1
			except:
				peso1=e[2]['weight']
				DTI_graph.add_weighted_edges_from([(nodo1,nodo2,peso1)])

	# Preparing the layout
	final_layout=nx.random_layout(DTI_graph)

	for n in DTI_graph.nodes():
		pos = DTI_graph.node[n]['mass center']
		XX = pos[0]
		YY = pos[1]
		ZZ = pos[2]
		final_layout[n] = [255-YY,ZZ]
	
	pesi = np.array([])
	for e in DTI_graph.edges(data=True):
		peso = e[2]['weight']
		pesi = np.hstack((pesi,[100*np.array(peso)]))
	
	normalization_factor=pesi.max()
		 
	pesi_normalizzati = pesi*(4/normalization_factor)
	
	# Drawing
	nx.draw(DTI_graph,final_layout,node_color=graph_node_color,node_size=80,edge_color=graph_edge_color,width=pesi_normalizzati,edge_cmap=plt.cm.Greys,with_labels=True)
	plt.xlim(0,255)
	plt.ylim(0,50)
	plt.axis('off')

# -------------------------------------------------------------------------------------------
# Show the graph hemisphere or the whole graph
# show_graph(graph,side)
# side: 0 whole, 1 left, 2 right)
def show_graph_v2(graphIN,side,graph_edge_colorIN,edge_styleIN,graph_node_colorIN,normalization_factorIN):
	import sys,os
	import networkx as nx
	import matplotlib.pyplot as plt
	import numpy as np

	DTI_graph = nx.Graph()
	
	# Input options
	if sys.argv>=4:
		graph_edge_color=graph_edge_colorIN
	else:
		graph_edge_color='black'
	
	if sys.argv>=5:
		graph_node_color=graph_node_colorIN
	else:
		graph_node_color='#A0CBE2'

	# Cropping the graph
	if side==0:
		for n in graphIN.nodes():
			DTI_graph.add_node(n)
			DTI_graph.node[n]['mass center']=graphIN.node[n]['mass center']
	elif side ==1:
		for n in graphIN.nodes():
			pos = graphIN.node[n]['mass center']
			if pos[0]<125:
				DTI_graph.add_node(n)
				DTI_graph.node[n]['mass center']=graphIN.node[n]['mass center']
	elif side ==2:
		for n in graphIN.nodes():
			pos = graphIN.node[n]['mass center']
			if pos[0]>=125:
				DTI_graph.add_node(n)
				DTI_graph.node[n]['mass center']=graphIN.node[n]['mass center']	
	
	# Adding the edges to the cropped graph
	# (Make the graph undirect by selecting the maximum weight between the two edges)
	for e in graphIN.edges(data=True):
		nodo1 = e[0]
		nodo2 = e[1]
		if (nodo1 in DTI_graph.nodes())and(nodo2 in DTI_graph.nodes()):
			try:
				peso2=DTI_graph.edge[nodo2][nodo1]['weight']
				peso1=e[2]['weight']
				if peso1>peso2:
					DTI_graph.edge[nodo2][nodo1]['weight']=peso1
			except:
				peso1=e[2]['weight']
				DTI_graph.add_weighted_edges_from([(nodo1,nodo2,peso1)])

	# Preparing the layout
	final_layout=nx.random_layout(DTI_graph)

	for n in DTI_graph.nodes():
		pos = DTI_graph.node[n]['mass center']
		XX = pos[0]
		YY = pos[1]
		ZZ = pos[2]
		final_layout[n] = [255-YY,ZZ]
	
	pesi = np.array([])
	for e in DTI_graph.edges(data=True):
		peso = e[2]['weight']
		pesi = np.hstack((pesi,[100*np.array(peso)]))
	
	if sys.argv>=6:
		normalization_factor=normalization_factorIN
	else:
		normalization_factor=pesi.max()
		 
	pesi_normalizzati = pesi*(4/normalization_factor)
	
	# Drawing
	nx.draw(DTI_graph,final_layout,node_color=graph_node_color,node_size=80,edge_color=graph_edge_color,style=edge_styleIN,width=pesi_normalizzati,edge_cmap=plt.cm.Greys,with_labels=False)
	plt.xlim(0,255)
	plt.ylim(0,50)
	plt.axis('off')

# -------------------------------------------------------------------------------------------

def show_graph_cliques(graphIN,side,length_cliques,cliques_colorIN):
	import sys,os
	import networkx as nx
	import matplotlib.pyplot as plt
	import numpy as np

	# Checking for input options
	if sys.argv>=3:
		cliques_color=cliques_colorIN
	else:
		cliques_colorIN='red'
	
	# Computing the normalizaiton factor
	pesi = np.array([])
	for e in graphIN.edges(data=True):
		peso = e[2]['weight']
		pesi = np.hstack((pesi,[100*np.array(peso)]))
	normalization_factor = pesi.max()
		
	# Showing the whole graph
	show_graph_v2(graphIN,side,'black','#A0CBE2',normalization_factor)
	
	# Highlighting the cliques
	for clique_edge_list in list(nx.find_cliques(graphIN)):
		if len(clique_edge_list)==length_cliques:
			clique_subGraph = graphIN.subgraph(clique_edge_list)
			show_graph_v2(clique_subGraph,side,cliques_color,'#A0CBE2',normalization_factor)

	
# -------------------------------------------------------------------------------------------
def draw_modules(graph_objIN,side,convert):
	import community as comm
	import random
	import matplotlib as mpl
	import numpy as np
	import networkx as nx
	import graph_lib as gl

	if convert:
		# Building the undirected graph
		print 'Converting the directed graph to an undirected graph...'
		graph_obj = nx.Graph()
		graph_obj.add_nodes_from(graph_objIN.nodes(data=True))
		for edge in graph_objIN.edges(data=True):
			nodo1 = edge[0]
			nodo2 = edge[1]
			if (nodo1 in graph_obj.nodes())and(nodo2 in graph_obj.nodes()):
				try:
					peso2=graph_objIN.edge[nodo2][nodo1]['weight']
					peso1=edge[2]['weight']
					if peso1>peso2:
						graph_obj.edge[nodo2][nodo1]['weight']=peso1
				except:
					peso1=edge[2]['weight']
					graph_obj.add_weighted_edges_from([(nodo1,nodo2,peso1)])
	else:
		print 'Using the directed graph...'
		graph_obj = graph_objIN

	# Computing the normalization factor
	print 'Computing the normalization factor...'
	pesi = np.array([])
	for e in graph_obj.edges(data=True):
		peso = e[2]['weight']
		pesi = np.hstack((pesi,[100*np.array(peso)]))	
	normalization_factor = pesi.max()
	
	# Computing the Modules
	# (Q_score,sub_graph_list) = comm.detect_communities(graph_obj)
	print 'Computing the modules...'
	sub_graph_list=gl.find_modules(graph_obj)
	
	# Computing the graph-edges not in a module
	remaining_graph = graph_obj.copy()
	for cont in range(len(sub_graph_list)):
		sub_graph = graph_obj.subgraph(sub_graph_list[cont])
		for edge_element in sub_graph.edges():
			if edge_element in remaining_graph.edges():
				remaining_graph.remove_edge(edge_element[0],edge_element[1])
	
	# Showing the remaining part of the graph
	print 'Drawing the part remaining part of the graph...'
	show_graph_v2(remaining_graph,side,'black','dashed','white',normalization_factor)

	# Showing modules
	print 'Drawing the modules...'
	for cont in range(len(sub_graph_list)):
		sub_graph = graph_obj.subgraph(sub_graph_list[cont])
		if cont == 0:
			# Module containing the isolated nodes
			color_RGB_edges = [1,1,1]  #[0.2, 0.2, 0.2]
			color_RGB_nodes = 'grey'
		else:
			color_RGB_edges = [random.random(), random.random(), random.random()]
			constant = 0.9
			color_RGB_nodes = mpl.colors.rgb2hex([constant*color_RGB_edges[0], constant*color_RGB_edges[1], constant*color_RGB_edges[2]])
			#color_RGB_nodes = '#A0CBE2'
		show_graph_v2(sub_graph,side,mpl.colors.rgb2hex(color_RGB_edges),'solid',color_RGB_nodes,normalization_factor)