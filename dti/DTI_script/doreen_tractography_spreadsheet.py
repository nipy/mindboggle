nodeList = ['1014','2014','26','58','16']

connection_graph_list = ['/mind_cart/DTI/50192/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50195/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50201/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50202/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50204/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50205/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50207/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50212/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50213/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50228/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50231/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50232/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50233/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50235/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50237/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50241/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50244/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50250/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50277/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50301/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50303/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50329/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50341/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50345/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50348/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50359/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50362/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50390/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50400/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50410/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50423/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50437/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50438/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50468/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50475/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50531/tract_doreen/connection_graph.gpi','/mind_cart/DTI/50552/tract_doreen/connection_graph.gpi']

correlation_graph_list = ['/mind_cart/DTI/50192/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50195/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50201/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50202/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50204/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50205/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50207/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50212/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50213/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50228/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50231/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50232/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50233/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50235/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50237/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50241/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50244/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50250/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50277/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50301/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50303/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50329/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50341/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50345/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50348/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50359/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50362/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50390/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50400/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50410/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50423/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50437/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50438/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50468/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50475/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50531/tract_doreen/correlation_graph.gpi','/mind_cart/DTI/50552/tract_doreen/correlation_graph.gpi']

nSubject = len(connection_graph_list)
weight_table = np.zeros((nSubject,25),'float32')
streamlines_table = np.zeros((nSubject,25),'float32')
voxel_table = np.zeros((nSubject,25),'float32')
subList = []
edge_list_text = []
for contSub in range(nSubject):
	subFile = correlation_graph_list[contSub]
	subID = subFile[15:20]
	subList.append(subID)
	graph = nx.read_gpickle(correlation_graph_list[contSub])
	edge_list = graph.edges()
	pos =0
	for contNode1 in range(len(nodeList)):
		node1 = nodeList[contNode1]
		for contNode2 in range(len(nodeList)):
			node2 = nodeList[contNode2]
			if ((node1 != node2)and((node1,node2)in edge_list)):
				weight_table[contSub][pos] = graph.edge[node1][node2]['weight']
				#streamlines_table[contSub][pos] = graph.edge[node1][node2]['n_streamlines']
				#voxel_table[contSub][pos] = graph.edge[node1][node2]['n_seed_voxel']
			edge_list_text.append(node1+' -> '+node2)
			pos = pos+1
			

nR = weight_table.shape[0]
nC = weight_table.shape[1]
testo_W = []
#testo_S = []
#testo_V = []
for r in range(nR):
	riga_W = ''
	#riga_S = ''
	#riga_V = ''
	for c in range(nC):
		riga_W = riga_W + str(weight_table[r][c]) + '\t'
		#riga_S = riga_S + str(streamlines_table[r][c]) + '\t'
		#riga_V = riga_V + str(voxel_table[r][c]) + '\t'
	testo_W.append(riga_W)
	#testo_S.append(riga_S)
	#testo_V.append(riga_V)
	
text_file_name = '/mind_cart/DTI/FA_spreadsheets/Doreen_correlation_tables.txt'
fobj = open(text_file_name,'w')
fobj.writelines(['CONNECTION GRAPH TRABLES\n'])
fobj.writelines(['\nSubject list\n'])
fobj.writelines(['%s\n' %(line) for line in subList])
fobj.writelines(['\nEdge list\n'])
fobj.writelines(['%s\n' %(line) for line in edge_list_text])
fobj.writelines(['\nEdge Weight\n'])
fobj.writelines(['%s\n' %(line) for line in testo_W])
#fobj.writelines(['\nEdge Number of Streamlines\n'])
#fobj.writelines(['%s\n' %(line) for line in testo_S])
#fobj.writelines(['\nEdge Number of Voxels\n'])
#fobj.writelines(['%s\n' %(line) for line in testo_V])
fobj.close()
