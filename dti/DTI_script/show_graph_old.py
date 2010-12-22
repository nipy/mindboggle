# Show the graph hemisphere or the whole graph
# show_graph(graph,side)
# side: 0 whole, 1 left, 2 right)
def show_graph(graphIN,side):
	import sys,os
	import networkx as nx
	import matplotlib.pyplot as plt
	import numpy as np

	DTI_graph = nx.Graph()

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
	
	pesi_normalizzati = pesi*(4/pesi.max())
	
	# Drawing
	nx.draw(DTI_graph,final_layout,node_color='#A0CBE2',node_size=80,edge_color='black',width=pesi_normalizzati,edge_cmap=plt.cm.Greys,with_labels=False)
	plt.xlim(0,255)
	plt.ylim(0,50)
	plt.axis('off')
