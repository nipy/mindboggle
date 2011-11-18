# Temp Queries
g.addEdge(g.v(2),g.v(39),'CONTAINS')

# Add subject specific nodes
# Subject1 from CUMC12
node = g.v(7)
g.addEdge(node,g.addVertex(["name":"Surfaces","type":"surfaces"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Pits","type":"pits"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Fundi","type":"fundi"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Basins","type":"basins"]),'CONTAINS')

# Add surface sub types to surfaces node
node = g.v(25)
g.addEdge(node,g.addVertex(["name":"White Matter Surface","type":"wm"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Pial Surface","type":"pial"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Inflated Surface","type":"inflated"]),'CONTAINS')

# Subject2 from CUMC12
node = g.v(8)
g.addEdge(node,g.addVertex(["name":"Surfaces","type":"surfaces"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Pits","type":"pits"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Fundi","type":"fundi"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Basins","type":"basins"]),'CONTAINS')

# Add surface sub types to surfaces node
node = g.v(32)
g.addEdge(node,g.addVertex(["name":"White Matter Surface","type":"wm"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Pial Surface","type":"pial"]),'CONTAINS')
g.addEdge(node,g.addVertex(["name":"Inflated Surface","type":"inflated"]),'CONTAINS')

# Add labeling protocols
node = g.v(2)
g.addEdge(node,g.addVertex(["name":"Label Protocol","type":"label_protocol"]),'CONTAINS')

node = g.v(39)
g.addEdge(node,g.addVertex(["name":"FreeSurfer","type":"freesurfer"]),'CONTAINS')

node = g.v(41)
g.addEdge(node,g.addVertex(["name":"unknown","type":"label","label_id":1]),'HAS LABEL')
g.addEdge(node,g.addVertex(["name":"bankssts","type":"label","label_id":2]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"caudalanteriorcingulate","type":"label","label_id":3]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"caudalmiddlefrontal","type":"label","label_id":4]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"corpuscallosum","type":"label","label_id":5]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"cuneus","type":"label","label_id":6]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"entorhinal","type":"label","label_id":7]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"fusiform","type":"label","label_id":8]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"inferiorparietal","type":"label","label_id":9]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"inferiortemporal","type":"label","label_id":10]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"isthmuscingulate","type":"label","label_id":11]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"lateraloccipital","type":"label","label_id":12]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"lateralorbitofrontal","type":"label","label_id":13]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"lingual","type":"label","label_id":14]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"medialorbitofrontal","type":"label","label_id":15]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"middletemporal","type":"label","label_id":16]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"parahippocampal","type":"label","label_id":17]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"paracentral","type":"label","label_id":18]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"parsopercularis","type":"label","label_id":19]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"parsorbitalis","type":"label","label_id":20]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"parstriangularis","type":"label","label_id":21]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"pericalcarine","type":"label","label_id":22]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"postcentral","type":"label","label_id":23]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"posteriorcingulate","type":"label","label_id":24]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"precentral","type":"label","label_id":25]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"precuneus","type":"label","label_id":26]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"rostralanteriorcingulate","type":"label","label_id":27]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"rostralmiddlefrontal","type":"label","label_id":28]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"superiorfrontal","type":"label","label_id":29]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"superiorparietal","type":"label","label_id":30]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"superiortemporal","type":"label","label_id":31]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"supramarginal","type":"label","label_id":32]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"frontalpole","type":"label","label_id":33]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"temporalpole","type":"label","label_id":34]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"transversetemporal","type":"label","label_id":35]),'HAS LABEL');
g.addEdge(node,g.addVertex(["name":"insula","type":"label","label_id":36]),'HAS LABEL');

# Add basins for subject 1
node = g.v(31);
g.addEdge(node,g.addVertex(["name":"basin","type":"basin","basin_id":1]),'CONTAINS');
g.addEdge(node,g.addVertex(["name":"basin","type":"basin","basin_id":2]),'CONTAINS');
g.addEdge(node,g.addVertex(["name":"basin","type":"basin","basin_id":3]),'CONTAINS');
# Add basins for subject 2
node = g.v(35);
g.addEdge(node,g.addVertex(["name":"basin","type":"basin","basin_id":1]),'CONTAINS');
g.addEdge(node,g.addVertex(["name":"basin","type":"basin","basin_id":2]),'CONTAINS');
g.addEdge(node,g.addVertex(["name":"basin","type":"basin","basin_id":3]),'CONTAINS');

# Add relationships
node = g.v(88);
g.addEdge(node,g.v(89),'IS SIMILAR');
g.addEdge(node,g.v(88),'IS SIMILAR');
g.addEdge(node,g.v(60),'HAS LABEL');
g.addEdge(node,g.v(64),'HAS LABEL');
g.addEdge(node,g.v(73),'HAS LABEL');
g.addEdge(node,g.v(52),'HAS LABEL');
g.addEdge(node,g.v(85),'HAS LABEL');
g.addEdge(node,g.v(76),'HAS LABEL');
g.addEdge(node,g.v(78),'HAS LABEL');
g.addEdge(node,g.v(82),'HAS LABEL');
g.addEdge(node,g.v(67),'HAS LABEL');
g.addEdge(node,g.v(75),'HAS LABEL');


node = g.v(89);
g.addEdge(node,g.v(75),'HAS LABEL');
g.addEdge(node,g.v(67),'HAS LABEL');
g.addEdge(node,g.v(73),'HAS LABEL');
g.addEdge(node,g.v(78),'HAS LABEL');
g.addEdge(node,g.v(52),'HAS LABEL');
g.addEdge(node,g.v(76),'HAS LABEL');
g.addEdge(node,g.v(64),'HAS LABEL');
g.addEdge(node,g.v(60),'HAS LABEL');
g.addEdge(node,g.v(85),'HAS LABEL');
g.addEdge(node,g.v(62),'HAS LABEL');
g.addEdge(node,g.v(82),'HAS LABEL');
