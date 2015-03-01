# The test code for Le Troter's sulcus pruning algorithm based on path probability

import cPickle
import pyvtk

def load_fundi(Filename):
    '''Load fundi from a VTK file
    
    
    Returns
    ==========
    
        Fundi : list of 2-tuples of integers
            Each element is a line segment in fundi, could belong to any component 

        Vertexes: list of 3-tuples of floats
            Each element is a vertex's coordinates 

    '''
    
    VTKReader = pyvtk.VtkData(Filename)
    Fundi =  VTKReader.structure.lines
    Vertexes =  VTKReader.structure.points
    return Fundi, Vertexes
        
def load_component(Filename):
    '''Load vertex components of a hemisphere
    
    Returns
    ==========
    
        VrtxCmpnts : list of lists of integers
            each element of VrtxCmpnt is a list of vertexes that are in the same sulcal component
     
    '''
    
    
    Fp = open(Filename, 'r')
    VrtxCmpnt = cPickle.load(Fp)
    Fp.close()
    return VrtxCmpnt

def fundus_per_component(VrtxCmpnts, Fundi):
    '''Return fundus vertexes within each component
    
    Parameters
    ===========
    
        VrtxCmpnts : list of lists of integers
            each element of VrtxCmpnt is a list of vertexes that are in the same sulcal component

        Fundi : list of 2-tuples of integers
            Each element is a line segment in fundi, could belong to any component 
    
    Returns
    =========
    
        Fundus_per_cmp : list of lists of 2-tuples of integers
            Each element is a list of fundus line segments within the same component. 
            Le Troter's algorithm runs on each of this. 
    
    '''
    
    Fundus_per_cmp = [[] for i in xrange(0, len(VrtxCmpnts)) ]
    for Fundus in Fundi:
        [V1, V2] = Fundus  # V1 and V2 must belong to the same component 
        for CmpntID, Cmpnt in enumerate(VrtxCmpnts):
            if V1 in Cmpnt:
                Fundus_per_cmp[CmpntID].append(Fundus)
                
    return Fundus_per_cmp    

def other(Tuple, Element):
    '''Given a 2-tuple and one element, return the other element. 
    '''
    
    if Element == Tuple[0]:
        return Tuple[1]
    else:
        return Tuple[0]
    
def init_graph(Fundus_per_cmp):
    '''Initialize the graph for running le troter's algorithm
    
    Parameters
    ===========
    
        Fundus_per_cmp : list of lists of 2-tuples of integers
            Each element is a list of fundus line segments within the same component.
        
        Degree : dictionary of integer keys and integer values
            Keys are nodes in initial graph and values are their degrees
            
        Merged : Boolean
            True is a short segment is merged with another short segment

        Nbr : a dictionary 
            key is node ID and value is IDs of neighbors of the node 

        Graph : a 4-tuple of 2 lists and 2 dictionaries 
            The 1st list includes all nodes in the initial graph, 
            while the 2nd list includes all edges in the initial graph.
            The 1st list is a list of integers and the 2nd is a list of 2-tuples of integers.
            The 3rd element is a *Degree* dictionary
            The 4th element is an *Nbr* dictionary 

        Edge_Comp : list of 2-tuples of tuples of integers
            Each element is a 2-tuple. The first is a 2-tuple containing two vertexes on the original mesh.
            During the running of this function, those two vertexes are terminal nodes of an edge under expansion.  
            After running this function, those two vertexes are nodes on initial graph. 
            The second is a list of integers representing a path on original mesh forming the expanding edge specified by the first tuple element.  
            
    Returns
    =========
    
        Graphs : list of Graph 
            Each element corresponds to one Graph of a sulcus patch/component 
                    
    Notes
    =========
    
        In the initial graph, there should be no node of the degree 2.
        
        How do we build edges for initial graph? We start with original segments on the mesh.
        For any two connected segments, if their joint has a degree of 2, replace them by a new edge spanning over their non-joint vertexes. 
        E.g., given A-B and B-C, replace them by A-C. 
     
    '''
    def find_Comp(Edge_Comp, Edge):
        '''The nested function that returns the path forming edge (N1, N2)
          
        '''
        for Comp in Edge_Comp:
            if Comp[0] == Edge:
                return Comp[1]
            
        print "cannot find such Edge in Edge_Comp", Edge, Edge_Comp
        return -1
        
    Graphs = []
    for Fundus in Fundus_per_cmp:
        if len(Graphs) ==5: 
            print "Check the creation of pgraph 6"
        Nodes, Edges = [], []
        Degree, Nbr = {}, {}
        Edge_Comp = []

        # Get all nodes in initial graph
        for Segment in Fundus: 
            Nodes += Segment
        Nodes = list(set(Nodes))
        # End of Get all nodes in initial graph
        
        # Get all edges in initial graph
        # 1. Find all nodes of degree 1 and greater than 2. 
        for Node in Nodes:
            Degree[Node] = 0
        for Segment in Fundus:
            [V1, V2] = Segment
            Degree[V1] += 1
            Degree[V2] += 1
        
        # 2. extend small segments into edges 
        for Segment in Fundus:
            [V1, V2] = Segment
            Merged = False
            if Degree[V1] == 2: # find another edge who has V1. 
                for EdgeID, Edge in enumerate(Edges):
                    if V1 in Edge:
                        Edges[EdgeID] = [other(Edge, V1), V2]
                        Merged = True 
                        
                        Old_Path = find_Comp(Edge_Comp, Edge)
#                        print Old_Path
                        if V1 == Old_Path[0]:
                            New_Path = [V2] + Old_Path
                        elif V1 == Old_Path[-1]:  
                            New_Path = Old_Path + [V2] 
                        else:
                            print "Error"
                            exit()
                        New_Edge =  [other(Edge, V1), V2]
                        Edge_Comp.remove([Edge, Old_Path])
                        Edge_Comp.append([New_Edge, New_Path])
                        
                        break 
            elif Degree[V2] == 2: # find another edge who has V1. 
                for EdgeID, Edge in enumerate(Edges):
                    if V2 in Edge:
                        Edges[EdgeID] = [V1, other(Edge, V2)]
                        Merged = True
                        
                        Old_Path = find_Comp(Edge_Comp, Edge)
                        if V2 == Old_Path[0]:
                            New_Path = [V1] + Old_Path
                        elif V2 == Old_Path[-1]:  
                            New_Path = Old_Path + [V1]
                        else:
                            print "Error"
                            exit()
                        New_Edge =  [V1, other(Edge, V2)]
                        Edge_Comp.remove([Edge, Old_Path])
                        Edge_Comp.append([New_Edge, New_Path])                      
                         
                        break
            if not Merged:  # there is not edge currently in Edges that shares terminal vertexes with this Segment
                Edges.append(Segment)
                Edge_Comp.append([Segment, Segment])
        # End of Get all edges in initial graph
        
        # 3. Get nodes and Degrees for initial graph
        NewNodes = [] 
        for Edge in Edges: 
            NewNodes += Edge
        NewNodes = list(set(NewNodes))
        
        NewDegree = {}
        for Node in NewNodes:
            NewDegree[Node] = 0
        for Edge in Edges:
            [V1, V2] = Edge
            NewDegree[V1] += 1
            NewDegree[V2] += 1
        
        # 4. generate neighbor list for initial graph
        for Node in NewNodes:
            Nbr[Node] = []
             
        for Edge in Edges: 
#            print Edge
            [V1, V2] = Edge 
            Nbr[V1].append(V2)
            Nbr[V2].append(V1)
        # End of generate neighbor list for initial graph
        
        Graph = [NewNodes, Edges, NewDegree, Nbr, Edge_Comp]
        Graphs.append(Graph)
        
    return Graphs
    
    
def path_prob2(Graphs):
    '''Compute path probability using Dijkstra algorithm  
    '''
    pass
    
def path_prob(Graphs):
    '''Compute path probability using my own algorithm
    
    Parameters
    ===========
        Graphs : list of Graph 
            Each element corresponds to one *Graph* of a sulcus patch/component 
            See definition of the same variable in function init_graph 
    
        UnAssigned : queue of integers 
            A queue of nodes who has not been assigned path probability (=0)
            
        Degree_1 : integer
            number of nodes of degree 1 in a graph. 
    
    Returns
    ========
    
        Path_Prob : dictionary
            keys are node ID and values are their path probability  
    
    Notes
    =======
        For temporary visualization purpose, we assign edge probability to all fundus vertexes on the edge. 
    
    '''
    Path_prob = {}
    
    for Graph in Graphs:        
        import collections
        UnAssigned = collections.deque([])
        Path_prob_local = {} # Path_prob for nodes in this graph 
        [Nodes, Edges, Degree, Nbr, Edge_Comp] = Graph
        
        Counter = 0 # number of iteration has been run on this graph  
        
#        Degree_1= 0
#        for Node in Nodes:
#            if Degree[Node] == 1:
#                Degree_1 += 1
                
        # 1. initialize Path_prob 
        for Node in Nodes:
            if Degree[Node] == 1:
                Path_prob_local[Node] = 1 #it is debating whether to let it be len(Nodes) or Degree_1 or just 1. 
            else:
                Path_prob_local[Node] = 0
                UnAssigned.append(Node)
                
        # 2. assigned Path_prob to nodes of degree greater than 3 
        while len(UnAssigned) >0:
            Candidate = UnAssigned.popleft() # pick a node from the list of unassigned nodes
            
            Counter += 1
            
            Some_UnAssigned = False # a flag denoting whether some lower-degree neighbors of Candidate are not assigned  
            for Node in Nbr[Candidate]:
                if Node in list(UnAssigned) and Degree[Node] < Degree[Candidate]:
                    Some_UnAssigned = True 
                    break
                
#            UnAssigned.remove(Candidate)
            if not Some_UnAssigned: 
                for Node in Nbr[Candidate]:
                    Path_prob_local[Candidate] += Path_prob_local[Node]
            else: # put Candidate to the back  of the UnAssigned list
                UnAssigned.append(Candidate) 
            
#            print UnAssigned
            
            if Counter > len(Nodes)*3:
                print "Too many iterations "
                break
            
        # 3. Assign path probability to vertexes between any pair of nodes (i.e., edges)
        # This step is only for visualization. 
        for Comp in Edge_Comp:
            [N1, N2] = Comp[0] 
            Path = Comp[1]
#            print Path
            Avg_Prob = (Path_prob_local[N1] + Path_prob_local[N2])/2.0 
            for i in xrange(1,len(Path)-1):
                Path_prob_local[Path[i]] = Avg_Prob 
                
        Path_prob.update(Path_prob_local)
        
    return Path_prob     

    
def Prob2File(Path_prob, Fundi, Vertexes):
    '''Visualize Path Probability and other debugging information to fundi
    '''
    
    Structure = pyvtk.PolyData(points=Vertexes, lines=Fundi)
    # convert dictionary-type Path_prob into a list for pyvtk to use
    Path_prob_list = [0 for i in xrange(len(Vertexes))]
    for Key,Value in Path_prob.iteritems():
        Path_prob_list[Key] = Value
    
    Pointdata = pyvtk.PointData(pyvtk.Scalars(Path_prob_list,name='Path_Prob'))
#                                pyvtk.Scalars(Parts,name='hierarchy'))
    pyvtk.VtkData(Structure, Pointdata).tofile('path_prob_test.vtk','ascii')
    
def prune(Graphs, Path_prob):
    '''Prune fundi according to path probability
    
    We prune branches (one end has degree of 1) in ascending order of path probability, 
    until there is no branching nodes on the fundi. 
    
    Parameters
    ============
    
        Graphs : list of Graph 
            Each element corresponds to one Graph of a sulcus patch/component
            Check the output from function init_graph 
    
        Path_prob : dictionary
            keys are node ID and values are their path probability
    
    Return 
    ==========
        NewGraphs : list of Graph
            Same structure as Graphs, but some branches are removed. 
     
    Notes
    =======
    
        2012-05-12 Currently we may not prune enough. 
     
    '''
    NewGraphs = list(Graphs) # deep copy, branch remove on this. 
    for GraphID, Graph in enumerate(Graphs):
        if GraphID == 6:
            print "debug point"
        
        [Nodes, Edges, Degree, Nbr, Edge_Comp] = Graph
        
        if Nodes == []:  # this graph/component has no fundus
            continue 
        
        # 1. form a dictionary of Path_prob for nodes on initial graph only. 
        
        Path_prob_local = {}
        for Node in Nodes:
            Path_prob_local[Node] = Path_prob[Node]
            
        # 2. sort the local Path_prob dict and form a stack of Nodes,
        #    where the node of high path prob are in the bottom. 
        from operator import itemgetter
        Sorted_Path_prob_list =  sorted(Path_prob_local.items(),key=itemgetter(1))
        from collections import deque   
        Exam = deque([])
        for Pair in Sorted_Path_prob_list:
            if Degree[Pair[0]] ==1:
                Exam.append(Pair[0])
        
        # 3. Form a list of nodes whose degree is greater than 2.  
        Forks = []
        for Node in Nodes:
            if Degree[Node] > 2:
                Forks.append(Node)
            
        if len(Forks) == 1:
            continue
            
        # 4. remove nodes until this graph has no node of degree greater than 2
        Counter  = 0 
        
        while len(Forks) > 0:
            Candidate = Exam.popleft()
            for Edge in Edges:
                if Candidate in Edge:
                    NewGraphs[GraphID][1].remove(Edge)  # 1 is Edges
                    Fork = other(Edge, Candidate)
                    NewGraphs[GraphID][2][Fork] -= 1  # 2 is Degree
                    if  NewGraphs[GraphID][2][Fork] <3:  # should be < 3 
                        Forks.remove(Fork)
                    if NewGraphs[GraphID][2][Fork] ==1: 
                        Exam.append(Fork)
                    break 
                Counter += 1
            if Counter > len(Nodes)+200:
                print "Something may go wrong in pruning"
                break
            
    return NewGraphs

def getComp(Edge_Comp, Edge):
    '''Return components of an edge 
    
    Parameters
    ============
    
    Edge_Comp : list of 2-tuples of tuples of integers
            Each element is a 2-tuple. The first is a 2-tuple containing two vertexes on the original mesh.
            During the running of this function, those two vertexes are terminal nodes of an edge under expansion.  
            After running this function, those two vertexes are nodes on initial graph. 
            The second is a list of integers representing a path on original mesh forming the expanding edge specified by the first tuple element.
            
    Edge : list of 2 integers
        
        
    
    Returns
    ===========
    
    Comp : list of integers
    
    '''
    
    for [Edge_s, Comp] in Edge_Comp:
        if Edge_s == Edge:
            return Comp
        
    return -1
    
def graph2File(Graphs,Vertexes):
    '''Write Graphs to a file for visualizing pre- or post-prunning fundi. 
    '''
    Lines = []
    for Graph in Graphs:
        [Nodes, Edges, Degree, Nbr, Edge_Comp] = Graph
        for Edge in Edges:
            Comp = getComp(Edge_Comp, Edge)
            for i in xrange(len(Comp)-1):
                Lines.append([Comp[i], Comp[i+1]])
    
    Structure = pyvtk.PolyData(points=Vertexes, lines=Lines)
    # convert dictionary-type Path_prob into a list for pyvtk to use
    pyvtk.VtkData(Structure).tofile('new_fundi.vtk','ascii')

def test_run(CmpntFile, FundiFile):
    '''Run letroter algorithm on a test hemihsphere 
    '''
    Fundi, Vertexes = load_fundi(FundiFile)
    VrtxCmpnts = load_component(CmpntFile)
    Fundus_per_cmp = fundus_per_component(VrtxCmpnts, Fundi)
    Graphs = init_graph(Fundus_per_cmp)
    Path_prob = path_prob(Graphs)
#    Prob2File(Path_prob, Fundi, Vertexes)
    NewGraphs = prune(Graphs, Path_prob)
    graph2File(Graphs, Vertexes)
    
    return Path_prob, Graphs, NewGraphs

import sys
Path_prob, Graphs, NewGraphs = test_run(sys.argv[1], sys.argv[2])