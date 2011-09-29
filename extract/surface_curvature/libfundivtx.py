# All functions for extracting fundi as vertex clouds 

import fileio, libvtk, libbasin
from numpy import mean, std, abs, matrix, zeros, flatnonzero, sign, array, argmin
import os, sys

sys.setrecursionlimit(10000)

def judgeNode1(V0, CurvatureDB, Threshold=0):
    """Check whether a vertex (node) satisfies the zero-order criterion that its curvature value is negative
    
    Input
    ======
    
        V0: integer
                the ID of a vertex, indexing from 0
                
        CurvatureDB: list
                len(CurvatureDB) == number of vertexes in the hemisphere
                CurvatureDB[i]: integer, the curvature of the i-th vertex 
    
    """
    
#    if abs(CurvatureDB[V0] - Threshold) < 0.01:  # for output curvature contour
    if CurvatureDB[V0] > Threshold:
        return True
    else:
        return False

def getNeighbor(V0, NbrDB=[], NbrFile = '', FaceDB=[]):
    """ Given a vertex, return all its neighbors as a list. 
    
    Input
    ======
    
        V0: integer
                the ID of a vertex, indexing from 0
                
        NbrDB/FaceDB: list
                len(DB) == number of vertexes in the hemisphere
                DB[i]: integer, the curvature of the i-th vertex 

        NbrFile: string
                the file name of a file storing neighbors of vertexes
    
    """ 
    
    Nbr = []
    
#    for FaceID, Face in enumerate(DB):
#        Tmp = list(Face)   # I wanna remove elements from Face, but I do not wanna change it                
#        if V0 in Tmp:
#            Tmp.remove(V0)
#            [V1, V2] = Tmp
#            if not (V1 in Neighbor):
#                Neighbor.append(V1)
#            if not (V2 in Neighbor):
#                Neighbor.append(V2)
#            if len(Neighbor) == 3:
#                return Neighbor

    if NbrDB != []:
        return NbrDB[V0]
                
    return Nbr # if this face has no neighbor co-edging the two vertexes
    
    
def judgeNode2(V0, CurvDB, NbrDB, Lb = -0.03, Ub = 0.03, Level = 1):
    """Check whether a vertex satisfies the first-order criterion
    
    Input
    ======
        V0: integer
                the ID of a vertex, indexing from 0
        
        Neighbor: list
                all neighbors of vertex V0
                
        CurvDB: list
                len(CurvatureDB) == number of vertexes in the hemisphere
                CurvatureDB[i]: integer, the curvature of the i-th vertex 

        NbrDB: list
                len(NbrDB) == number of vertexes in the hemisphere
                NbrDB[i]: list, the neighbors of the i-th vertex                

        Lb/Ub: float
                the lower or upper boundary of the derivative test
                
        Level: integer
            There are many ways to judge whether a vertex can satisfy the 1st-order derivative test, e.g., comparing each neighbor
            with the vertex or comparing the average value of all neighbors with the vertex. Level allows us to choose among different
            ways. 
            
            Level = 1: A face has to satisfy a condition with each of its neighbors/directions
            Level = 2: A face has to satisfy a condition with the average of all its neighbors/directions 
    
    """
    Neighbor = NbrDB[V0]
    
    
    if len(Neighbor)==0 :
        return False
    else:
        Diff = [CurvDB[i] - CurvDB[V0] for i in Neighbor]

# option 1: requires all neighbors have higher value than the node    
    if Level == 1:    
        Judge = [Ub > D > Lb for D in Diff ] 


        for J in Judge:
            if not J:
                return False
        return True

# option 2: only require 80% of neighbors have higher value than this node.
    elif Level == 2:
        Judge = [Ub > D > Lb for D in Diff ] 
        FalseCount = 0
        for J in Judge:
            if not J:
                FalseCount += 1
            
        if FalseCount <= 0.5 * len(Neighbor):
            return True
        else:
            return False
    
#option 3: requires the average value of neighbors be highre than the node
    elif Level == 3:
        Diff = mean(Diff)
    
        if Ub > Diff > Lb:
            return True
        else:
            return False
    

def judgeNode3(V0, CurvDB, NbrDB, Lb = 0, Level = 1):
    """Check whether a vertex satisfies the second-order criterion
        
    Input
    ======
        V0: integer
                the ID of a vertex, indexing from 0
        
        Neighbor: list
                all neighbors of vertex V0
                
        CurvDB: list
                len(CurvatureDB) == number of vertexes in the hemisphere
                CurvatureDB[i]: integer, the curvature of the i-th vertex 

        NbrDB: list
                len(NbrDB) == number of vertexes in the hemisphere
                NbrDB[i]: list, the neighbors of the i-th vertex                

        Lb: float
                the lower boundary of 2nd-order derivative test
        
    """
    
    FirstNbr = NbrDB[V0] # first-hop neighbors of V0
    SecondNbr = [NbrDB[i] for i in FirstNbr] # second-hop neighbors of V0
    AvgSecond = [mean([CurvDB[i] for i in Nbr]) for Nbr in SecondNbr]# average curvatures of second-hop neighbors
    
    Diff = [AvgSecond[i] - 2*CurvDB[FirstNbr[i]] + CurvDB[V0] for i in xrange(0,len(FirstNbr))]

# option 1: requires all second-order derivatives to be greater than Lb 
    if Level ==1 :
        Judge = [ D > Lb for D in Diff ] 

        for J in Judge:
            if not J:
                return False

        return True

# option 2: requires the average of all second-order derivative to be greater than Lb
    elif Level == 2: 
        Judge = [D > Lb for D in Diff ] 
        FalseCount = 0
        for J in Judge:
            if not J:
                FalseCount += 1
            
        if FalseCount <= 0.5 * len(FirstNbr):
            return True
        else:
            return False
    elif Level == 3:
        Diff = mean(Diff)
    
        if Diff > Lb:
            return True
        else:
            return False
    
def rawFundi(VertexNum, CurvDB, NbrLst, Skip3 = True, Skip2=False):
    """Given a list of Vertexes, Faces and curvature values on vertexes, find out vertexes (as two lists) that are either candidate or strict fundi
    
    Notes
    =====
    
    Desent fundus vertex clouds can be obtained by the following settings:
    
    1. t1Lb, t2Lb, t2Ub, t3Lb = mean(CurvDB) + 1*std(CurvDB), -1*std(CurvDB)*0.02, 1*std(CurvDB)*0.02, -1*std(CurvDB)*0.01
       and
       t2@L3, t3@L3   
    
    """
    
    print "Extracting fundus vertex clouds"
    
    Strict, Candidate = [], [] # to store IDs for strict and candidate faces 
    
#    t1Lb, t2Lb, t2Ub, t3Lb = mean(CurvDB) + 1*std(CurvDB), -1*std(CurvDB)*0.02, 1*std(CurvDB)*0.02, -1*std(CurvDB)*0.01
#    t1Lb, t2Lb, t2Ub, t3Lb = 0.15, -1*std(CurvDB)*0.01, 1*std(CurvDB)*0.01, -1*std(CurvDB)*0.005
    t1Lb, t2Lb, t2Ub, t3Lb = mean(CurvDB)+1*std(CurvDB), -1*std(CurvDB)*0.01, 1*std(CurvDB)*0.01, -1*std(CurvDB)*0.005

#    t1Lb, t2Lb, t2Ub, t3Lb = 0, -1*std(CurvDB)*0.01, 1*std(CurvDB)*0.01, -1*std(CurvDB)*0.005

    
#empirical result: 
# if Level = 3, the mean approach 
# 1. t2Lb and t2Ub has to be 0.01*std for option 3 in t2, o/w, fundi areas are not included in strict set. Used to be 0.02. But that's visual illusion. 
# 2. t1Lb has to be 1*std, cann't be higher. O/w, 
    
#    t1Lb, t2Lb, t2Ub, t3Lb = mean(CurvDB) + 1*std(CurvDB), 0, 1*std(CurvDB)*0.1, 0
#    t1Lb, t2Lb, t2Ub, t3Lb = std(CurvDB) + min(CurvDB), -1*(max(CurvDB) - min(CurvDB))*0.01, 1*(max(CurvDB) - min(CurvDB))*0.01, -1*(max(CurvDB) - min(CurvDB))*0.01
#    print t1Lb, t2Lb, t2Ub, t3Lb#, max(CurvDB), min(CurvDB), std(CurvDB)
    
    for i in xrange(0, VertexNum):
#        if judgeNode1(i, CurvDB) and judgeNode2(i, CurvDB, NbrLst, Lb = -0.03, Ub = 0.03):
#            if judgeNode3(i, CurvDB, NbrLst, Lb = -0.017):
        if judgeNode1(i, CurvDB, t1Lb): 
            if Skip2:
                Candidate.append(i)
            else:
                if judgeNode2(i, CurvDB, NbrLst, Lb = t2Lb, Ub = t2Ub, Level = 3):
                    if Skip3:
                        Candidate.append(i)
                    else:
                        if judgeNode3(i, CurvDB, NbrLst, Lb = t3Lb, Level = 3):
                            Strict.append(i)
                        else:
                            Candidate.append(i)

                
    return Strict, Candidate

def lineUp(Special, NbrLst, VrtxCmpnts, VtxCoords, CurvFile, CurvatureDB):
    '''Form an MST of a set of vertexes and prune unwanted subtrees
    
    Parameters
    ===========
    
        Special : list of integers
            A list of special vertexes to be connected by MST 
            It can be all vertexes on original mesh
            
        NbrLst : list of list of integers
            neighbor list of vertexes
            
        VrtxCmpnts : list of list of integers
            each element of VrtxCmpnt is a list of vertexes that are in the same sulcal component
      
        Curve : list of integers
            vertexes on the curve segment form fundi  NOT LINED UP really now 04/12/2011
            
        VtxCoords : list of 3-tuples of doubles
            each element is the coordinate of a vertex element 
            
        CurvatureDB : list of doubles
            each element is the curvature value of a vertex 
            
        Links : list of list of 2-tuples of integers
            each element is a list of 2-tuples containing vertexes (in global ID) at two ends of an edge on MSTs
            Since the original mesh has been downsampled, each link in *Links* is a real edge on the mesh. 

        SpecialGroup: list of list of integers 
            each element is a list of Special vertexes in each connected component

        Ring : list of integers
            IDs of vertexes that have 0-curvature
            
        Curv: float
            curvature value of a vertex
            
        Idx: integer
            ID of a vertex
      
    Notes
    =======
    
        last updated: 2011-07-17 Forrest
        The function name lineUp does NOT reveal its purpose. Originally, this function is used to connect
        special vertexes via MST algorithm. But now MST is spanned on both special vertexes and common
        vertexes.  

    '''

    # Step 1: downsample the graph 
    Ring = []
    for Idx, Curv in enumerate(CurvatureDB):
        if abs(Curv) <0.01:
            Ring.append(Idx)
    
    NewVrtxCmpnts, NewNbrLst, SpecialGroup = downsample(Special+Ring, VrtxCmpnts, NbrLst, 0.5) 
        
    # step 2: prepare the distance matrix for FUNDUS vertex, for each fundus vertex, only 2 shortest edges are left
    #         Columns and rows for non-fundus vertexes are all zeros.
    DistFile = CurvFile + '.dist'
   
    Dists = dist(NewVrtxCmpnts, NewNbrLst, CurvatureDB, DistFile=DistFile)

    # Dists = zeroout(Dists, VrtxCmpnt, FundiList)  # only needed if MST only spans over special vertexes
    # Now Dists are weighted adjacency matrix of fundus vertexes in each component
    
#    fileio.wrtLists(DistFile+'.reduced', Dists)  # optional line, for debugging
#    Dists = fileio.readLists(DistFile+'.reduced')  # for fast debugging   only 
    # End of step 2
    
    # Step 3: use MST to connect all vertexes in NewVrtxCmpnts and dump into VTK format
    #FundusCmpnts = filterCmpnt(Nodes, VrtxCmpnts) # no need if Nodes are all vertexes on the mesh
    #FundusCmpnts = VrtxCmpnts # if MST will span via all nodes on the mesh
    #Links = mst(Dists, VrtxCmpnts, FundusCmpnts, NbrLst)  # deactivated Forrest 2011-07-21 because filters is done in downsample() 

    Links = mst(Dists, NewVrtxCmpnts, SpecialGroup, NbrLst) # Running version 
    #Links = mst(Dists, NewVrtxCmpnts, NewVrtxCmpnts, NbrLst) # debugging version
        
    return Links, NodeColor
    
def prune(Path, Degree, TreeNbr, Terminal, Branching, Special, VrtxCmpnt):
    '''Prune an MST by deleting edges
    
    Parameters
    ===========
    Path : list of lists (2-tuple) of integers
        Each element of *Path* is a list of the two terminals of each pair of connected fundus vertexes
        Vertexes indexed LOCALLY, i.e., they are referred by their ID in currenct component   
    
    Degree : list of integers 
        Degrees of nodes in a component
        
    TreeNbr : list of list of integers
        Each element is a list of neighbors of a node. All LOCAL IDs.
        
    Terminal : list of integers
        Local IDs of nodes that are terminals. 
        
    Branching : list of integers
        Local IDs of nodes that are branching nodes. 
    
    Special : list of integers
        Local IDs of nodes that are special vertexes connected by MST 
        
    Trace : list of 2-tuples of integers
        Edges that are visited along the Path from one node to another
        
    Visited : list of integers
        Nodes that are visited along the Path from one node to another
        
    At : integer
        a node
        
    Previous : integer
        a node
        
    VrtxCmpnt: list of integers
        A list of vertexes in this component in global ID
        
    Returns
    ========
    Path : list of lists (2-tuple) of integers
        Each element of *Path* is a list of the two terminals of each pair of connected fundus vertexes
        Vertexes indexed LOCALLY, i.e., they are referred by their ID in currenct component
        This one is reduced 
        
    ''' 
#    print "\t # of terminals:", len(Terminal)
    for T in Terminal:
#        print "\t\t\t Tracing begins at ", T, 
        if len(Branching) < 1:
            break
        Trace= [ ]# store the trace visited from a terminal node to a branching node
        Visited = [] # store nodes that have been visited  
        At = T # the node of current position in tracing
        while(not At in Branching and not At in Special):
            Visited.append(At)
            for Nbr in TreeNbr[At]:
                if not Nbr in Visited: # search toward the mainstream 
                    Trace.append((Nbr, At))
                    Previous = At  # the node visited before At in tracing
                    At = Nbr
                    break
        if not At in Special: # after the while loop, At stops at a branching node.
            TreeNbr[At].remove(Previous)     
            for Pair in Trace:
                (Src, Dst) = Pair
                if Pair in Path:
                    Path.remove(Pair)
                else:  # it is possible the order of nodes is reversed in Path
                    Path.remove((Dst, Src)) 
    return Path
    
def bfs(Seed, FundiList, NbrLst, Cmpnt):
    '''Use BFS to find the nearest vertex that is also in the cloud and in the same basin component. 
    Record the path which will be returned to lineUp(). The path is the fundus segment between the 
    seed to its nearest fundus neighbor on the surface. A special data structure is maintained to track
    the grow of bfs. 
    
    Parameters
    ===========
    
    Seed : integer
        the vertex to start/continue bfs
        
    FundiList : list of integers
        Vertexes that are in skeletonied vertex cloud of fundi
            
    NbrLst : list of list of integers
        neighbor list of vertexes
            
    Cmpnt : list of integers
        each element of VrtxCmpnt is a list of vertexes that are in the same basin component
    
    Support : list of integers
        Support[i] is the node from which BFS first visit node Cmpnt[i] (Cmpnt[i] is an id for a vertex in a hemisphere.)
        Once BFS reaches another fundus vertex from the seed vertex, we can extract the path by tracing back using Support.
        
    Path : list of integers
        Each integer in *Path* represents a vertex on the shortest path from *Seed* to any vertex in *FundiList*.
     
    '''
    Support = [None for i in xrange(0, len(Cmpnt))]
    NotVisited = list(Cmpnt)
    Queue = [Seed]
    Path = []
    PreviousNode = Seed
    while Queue != []:
#        print Seed, Queue
        Visit = Queue.pop()     # the node where BFS is at now. 
#        NotVisited.remove(Visit)
        if Visit in Cmpnt:
            if Visit in NotVisited:
                NotVisited.remove(Visit)
                Nbrs = NbrLst[Visit]
                for Nbr in Nbrs:
                    if Nbr in NotVisited and not Nbr in Queue: 
                        Support[Cmpnt.index(Nbr)] = Visit
                        Queue.insert(0, Nbr)
                        if Nbr in FundiList:
                            Queue = [] # leave the BFS to build the path
                            Path = [Nbr] 
                            PreviousNode = Visit
                            break
    
    while PreviousNode != Seed:
        Path.append(PreviousNode)
        CheckPoint = PreviousNode
        PreviousNode = Support[Cmpnt.index(CheckPoint)]
        
    Path.append(PreviousNode)  # A Path is a list of vertexes. fileio.writeFundiSeg() will break vertexes strings into edges on cortical surface  
    
    return Path


# the function below is deactivated by Forrest on 2011-07-10
#def bfsReach(Edges, Cmpnts, NbrLst):  # activated Forrest 2011-05-25 23:23
#    Using BFS to find the shortest (in Hamiltonian distance) path for all fundus vertexes pairs in Edges
#    
#    Edges : list of list (3-tuple) of integers
#        Each element is a 3-tuple. the first two are vertex IDs. 
#        Followed is the length of the distance between them. therefore, BFS should stop when the length is reached.
#        The last is the ID of the component.  
#
#    NbrLst : list of list of integers
#        neighbor list of vertexes 
#    
#    Cmpnts : list of lists of integers
#        Cmpnts[i] is a list of vertexes in the i-th connected component
    
    
#    Curve = []
#    
#    for Edge in Edges:
#        (Src, Dst, Cost, CID) = Edge
#        if Cost ==1:
#            Curve.append([Src, Dst])
#        else:
#            Curve.append( bfs(Src, [Dst], NbrLst, Cmpnts[CID]) ) 
#    return Curve


def bfsReach(Edge, Cmpnt, NbrLst): # New version of bfsReach to connect only one MST-determined vertex pair  
    '''This is a new version of bfsReach, activated on 2011-07-10. An older version was activated on 2011-05-25.
    
    The 2011-07-10 version replaces the loop of 2010-05-25 version by the loop body in the older version.  
    
    Parameters
    ===========
    
    Edges : list of list (3-tuple) of integers
        Each element is a 3-tuple. the first two are vertex IDs. 
        Followed is the length of the distance between them. therefore, BFS should stop when the length is reached.
        The last is the ID of the component.  

    NbrLst : list of list of integers
        neighbor list of vertexes 
        
    Cmpnt : list of integers
        a list of vertex IDs in this coonnected component
        
    Curve : list of list of integers
        Each element of *Curve* is a list of vertex IDs forming one branch on the MST
         
        
    ''' 
    
    (Src, Dst, Cost, CID) = Edge
    if Cost ==1:
        return [Src, Dst]
    else:
        return bfs(Src, [Dst], NbrLst, Cmpnt)
#    return [Src, Dst] 
    
def filterCmpnt(Vrtxs, Cmpnts):
    '''Leave vertexes in Cmpnts that are also in Vrtxs only. Thus, use Vrtxs (vertexes of certian property, i.e., fundi) to filter Cmpnts. 
    
    One useful application of this function is to remove vertexes in Cmpnts that are not fundus vertexes, which are in Vrtxs. 
    
    Notes
    ======
    
    Called turnCmpnt before 2011-05-21 17:50
     
    '''
    NewCmpnts = []
    for CID, Cmpnt in enumerate(Cmpnts):
        NewCmpnts.append([])
        for Vrtx in Cmpnt:
            if Vrtx in Vrtxs:
                NewCmpnts[-1].append(Vrtx)
    
    return NewCmpnts

def dist(VrtxCmpnts, VrtxNbrLst, CurvatureDB, DistFile=''): 
    '''Compute/load weighted adjacency matrixes of all connected sulcal components. 
        
    Parameters
    ===========
    
    Adj : list of list of integers/doubles
        adjacency matrix of a connected component
        
    Dist : list of list of integers 
        Dist[i][j] is the shortest distance between vertex i and vertex j
        Dist[i] is the Adj of the i-th connected component
        
    DistFile : string
        A path to filed extracted distance matrix
        
    CID : integer
        a component ID
    
    Cmpnt : list of integers
        global vertex IDs of all vertexes in a connected sulcal component 
        
    VrtxIdx : integer
        local (inner-componnent) ID of a vertex
         
    NbrIdx : integer
        local (inner-componnent) ID of a vertex's Nbr 
        
    CurvatureDB : list of doubles
            each element is the curvature value of a vertex 
    
    Notes
    =======
    
    Before HBM, this function uses graph power to weigh links between nodes. 
    Now (2011-07-17) it uses curvature to weigh. 
      
    '''    
   
    print "computing/loading weighted adjacency matrixes"
 
#    if DistFile != '' and os.path.exists(DistFile):
#        return fileio.readFltLsts(DistFile)
    # the block above is disabled because now computing adjacency matrix isn't time-consuming.
    # Forrest 2011-07-17
    
    Dists = [] 
    for CID, Cmpnt in enumerate(VrtxCmpnts):
        Num = len(Cmpnt)
        print "\t component", CID+1, ": size", Num
        if Num > 1:
            
            Adj = matrix(zeros((Num, Num)))
            for VrtxIdx, Vrtx in enumerate(Cmpnt):
                for Nbr in VrtxNbrLst[Vrtx]:
                    if Nbr in Cmpnt:
                        NbrIdx = Cmpnt.index(Nbr)                        
                        LinkWeight = -1. * (CurvatureDB[Vrtx]  + CurvatureDB[Nbr]) 
                         
                        # add a double check here to ensure the matrix is diagonally symmetric
                        if   Adj[VrtxIdx, NbrIdx] == 0:
                            Adj[VrtxIdx, NbrIdx] = LinkWeight # 
                            # Adj[NbrIdx, VrtxIdx] = LinkWeight # write only once for checking later
                        elif Adj[VrtxIdx, NbrIdx] != 0 and Adj[VrtxIdx, NbrIdx] != LinkWeight:
                            print "error, Adj is not symmetric."
                        elif Adj[NbrIdx, VrtxIdx] != 0 and Adj[NbrIdx, VrtxIdx] != LinkWeight:
                            print "error, Adj is not symmetric."
             
            #Dist = [[int(i) for i in Row] for Row in list(array(Adj))]
            Dist = [[i for i in Row] for Row in list(array(Adj))]
#            for Row in list(array(Dist)):
#                for Element in Row:
#                    print Element
        
        else:
            Dist = [[1]]            
            
        Dists.append(list(Dist)) # I forgot why this line is needed Forrest 2011-07-17
        
#    if DistFile != '':
#        fileio.wrtLists(DistFile, Dists)
    return Dists

def zeroout(Dists, VrtxCmpnts, Fundi):
    '''Zero out rows and columns for vertexes that are not considered as fundus vertexes, i.e., in *Fundi*
    '''
    
    if len(Dists) != len(VrtxCmpnts):
        print "Error"
        exit()
    else:
        for i in xrange(0,len(Dists)):
#            Dist = Dists[i]
#            Cmpnt = VrtxCmpnts[i]
            Num = len(VrtxCmpnts[i])
            Dist = array(Dists[i])
            for Idx, Vrtx in enumerate(VrtxCmpnts[i]):
                if not Vrtx in Fundi:
                    Dist[Idx, :] = zeros(Num)
                    Dist[:, Idx] = zeros(Num)
#            Dist = list(Dist)
            Dists[i] = [[int(D) for D in Row] for Row in list(Dist)]
        
    return Dists

def nonZeroLn(List): # activated 2011-05-25 19:14
    '''Given a 2-D list, return number of 1-D lists that contains nonzero elements
    '''
    Counter = 0
    for L in List:
        for E in L:
            if E != 0:
                Counter += 1
                break
            
    return Counter

def downsample(Special, VrtxCmpnts, NbrLst, Prob):
    '''Randomly delete vertexes on original mesh by probability Prob. 
    
    
    Parameters 
    =============
    
    Special : list of integers
        special vertexes that have to be in downsampled mesh. 
    
    VrtxCmpnts : list of lists of integers
        VrtxCmpnts[i] is a list of vertexes in the i-th connected component
        
    NbrLst : list of lists of integers
            neighbor list of vertexes
    
    Prob : float 
        The probability \in [0, 1] that a vertex is to be KEPT. 
        If it is 1, remove nothing.
        If it is 0, remove ALL. 
    
    N : list of lists of integers
        new vertexes to be left after downsampling
        
    L : list of lists of integers
        neighbor list of vertexes after downsampling
    
    Vrtx : integer
        a vertex id
    
    MinusVrtx : integer
        number of vertexes removed in mesh downsampling
    
    MinusEdge : integer
        number of edges removed in mesh downsampling
    
    Keep : list of integers
        vertexes to be kept 
    
    SpecialGroup: list of list of integers 
        each element is a list of Special vertexes in each connected component
    
    SpecialInThisGroup : list of integers
        a list of Special vertexes ID, used to represent all Special vertexes in a running component
    
    '''
    
    print "Downsampling mesh..."
    
    import random
    
    N = []
    L = [[] for i in xrange(0,len(NbrLst))]
    SpecialGroup = []
        

    for Cmpnt in VrtxCmpnts:
        for Vrtx in Cmpnt:
            L[Vrtx] = NbrLst[Vrtx]

    for CmpntID, Cmpnt in enumerate(VrtxCmpnts): # for each component
        MinusVrtx, MinusEdge = 0, 0
        NumMusthave = 0
        Keep = []
        SpecialInThisGroup = []
        for Vrtx in Cmpnt:   # for each vertex in the component
            if Vrtx in Special :
#                N[-1].append(Vrtx)
                NumMusthave += 1
                SpecialInThisGroup.append(Vrtx)
                if not Vrtx in Keep:
                    Keep.append(Vrtx)
            elif not Vrtx in Keep: # The purpose of not Vrtx in Keep is to avoid removing vertexes that should not be removed, such as Musthave.
                if Prob <= random.random():  # REMOVE Vrtx
                    MinusVrtx += 1
                    for Nbr in NbrLst[Vrtx]:
                        if Nbr in Cmpnt:
                            # step 1: put Vrtx's neighbors, which are ALSO in Cmpnt, into N and thus delete Vrtx
                            if not Nbr in Keep:
                                Keep.append(Nbr) # yield more vertexes and edges removed
                            # step 2: delete edges ending at Vrtx
                            if Vrtx in L[Nbr]:
        #                        N[-1].append(Nbr) # yield less vertexes and edges removed
        #                        Left.append(Nbr)  # yield less vertexes and edges removed
                                L[Nbr].remove(Vrtx)
                                MinusEdge += 1
                    L[Vrtx] = [] # Vrtx has no neighbor now
                else: # KEEP this vertex
                    Keep.append(Vrtx)
                    
        N.append(Keep)
        SpecialGroup.append(SpecialInThisGroup)
#        print "\t component", CmpntID+1, ":", NumMusthave, "Specials. ", "Vtx #:", len(Cmpnt), "-", MinusVrtx,  "=>", len(Keep)

    return N, L, SpecialGroup

class Prim:  # modified from the code without license @http://hurring.com/scott/code/python/mst_prim/v0.1/mst_prim.py 
    INFINITY = 2**8  # this is large enough for current problem size
    vertices = 0
    
    def __init__(self, A, r):
        """
        Prepare the inputs for mst_prime
        """
        self.vertices = A[0].__len__();
        self.nonzeros = nonZeroLn(A)   # a new member, activated Forrest 2011-05-25 18:56
        self.init_adjacency(A)
        self.remove_route(A, r)
        self.degree = [0 for i in xrange(0, len(A))] # a new member , activated Forrest 2011-05-27 20:31
        self.tree_nbr= [[] for i in xrange(0, len(A))] # record nbrs of each node, Forrest 2011-09-24
    
    def mst_prim(self, A, w, i, path, degree, tree_nbr):
        """
        'A' is the adjacency matrix
        'w' is the list of all connected vertices (in order of discovery)
        'path' is a list of tuples showing (from, to)
        i : the ID of the connected component   # Forrest 2011-05-26 00:31
        """
       
        # Stop when we've added all nodes to the path
#        if (w.__len__() == self.vertices):  # old line.  But if some nodes are not connected, it goes into infinite recursion. Deactivated Forrest 2011-05-25 19:39
        if (w.__len__() == self.nonzeros):  # new way, activated Forrest 2011-05-25 19:42
            return (A, w, path, degree, tree_nbr)
        # Find minimum path coming OUT of the known vertexes            
        (vfrom, vto, vcost) = self.find_min(A, w)
        
        # increase the degreee for vertexes vfrom and vto
        degree[vfrom] += 1
        degree[vto] += 1
        
        # update tree_nbr list for vfrom and vto Forrest 2011-09-24 10:55
        tree_nbr[vfrom].append(vto)
        tree_nbr[vto].append(vfrom)
        
        # Mark down this vertex as being a part of the MST path
        w.append(vto)
        #path.append((vfrom,vto,vcost, i)) # commented, Forrest 2011-09-24
        path.append((vfrom, vto))
        
        self.remove_route(A, vto)
        
        return self.mst_prim(A, w, i, path, degree, tree_nbr)
    
    
    def init_adjacency(self, A):
        """
        Initialize adjacency list - set 0 = INFINITY
        """
        for i in range(0, self.vertices):
            for j in range(0, self.vertices):
                if A[i][j] == 0:
                    A[i][j] = 2**8
        
    def remove_route(self, A, v):
        """
        Once we've added a node to our path, set all routes
        to this node equal to INFINITY - to prevent loops
        """
        for i in range(0, self.vertices):
            A[i][v] = self.INFINITY
    
    def find_min(self, A, w):
        """
        Find the cheapest connection we can possibly make,
        given the partially-built MST 'w'
        'vfrom' vertex to connect from
        'vto' vertex to connect to
        'vcost' cost of connection
        """
        vcost = self.INFINITY
        vto = vfrom = -1
        for v in w:
            # Get array offset of minimum  of this vertex
            i = argmin(A[v])
            if A[v][i] < vcost:
                vcost = A[v][i]
                vto = i
                vfrom = v
        return (vfrom, vto, vcost)

    # The end of Class Prim

def otherPaths(V0, APair, Degree, Pairs):
    '''Given a vertex *V0* and a path *APair*, find all other paths in *Pairs* that visits V0
    
    V0 : integer
        a vertex ID
        
    APair : a 4-tuple
    
    Degree : integer
        degree of V0
        
    Pairs : list of 4-tuples 
     
    '''
    Result = []
    Terminals = []
    for Pair in Pairs:
        if Pair[0] == V0 and Pair != APair:
            Result.append(Pair)
            Terminals.append(Pair[1])
            if len(Result) == Degree:
                return Result, Terminals
        elif Pair[1] == V0 and Pair != APair:
            Result.append(Pair)
            Terminals.append(Pair[0])
            if len(Result) == Degree:
                return Result, Terminals
            
#    print "found a standalone edge"  # This might be a potential problem. Need to check later. Forrest 2011-05-28 15:07 
    return Result, Terminals

def mst(Adjs, VrtxCmpnts, SpecialGroup, NbrLst):
    '''Using Prim algorithm to connect nodes (including fundus vertexes) within each connected component
    
    Parameters
    ==========
    
    VrtxCmpnts : list of lists of integers
        VrtxCmpnts[i] is a list of vertexes in the i-th connected component
    
    SpecialGroup : list of lists of integers
        Special[i] is a list of Special vertexes in the i-th connected component
        Special[i] can be empty 
    
    NbrLst : list of list of integers
            neighbor list of vertexes
    
    Path : list of lists (2-tuple) of integers
        Each element of *Path* is a list of the two terminals of each pair of connected fundus vertexes
        Vertexes indexed LOCALLY, i.e., they are referred by their ID in currenct component  
    
    Adj : list of lists of integers
        adjacency matrix of fundus vertexes in one connected component
        
    Adjs : list of Adj's
        adjacency matrixes of all fundus vertexes on one hemisphere
        
    W : list of integers
        vertexes taht are already connected in Prim algorithm
        It has no use.  
        
    Segs : list of lists (2-tuple) of integers
        Vertexes indexed GLOBALLY, i.e., they are referred by their ID in the entire hemisphere
        
    Degree : list of integers 
        Degrees of nodes in a component
        
    TreeNbr : list of list of integers
        Each element is a list of neighbors of a node. All LOCAL IDs.
        
    Terminal : list of integers
        Local IDs of nodes that are terminals. 
        
    Branching : list of integers
        Local IDs of nodes that are branching nodes. 
        
    
     
      
    '''
    Segs = []
    
    if len(Adjs) != len(VrtxCmpnts):
        print "Error, Adjs is not as long as VrtxCmpnts"
        exit()
    else:
        print "Connecting fundus vertexes in", len(VrtxCmpnts), "connected components."
        for i in xrange(0, len(Adjs)):  # For each component in the hemisphere
            print "\t MST on component",i+1, ",",
            if len(SpecialGroup[i]) < 2 :  # This compnent has no more than two vertexes to be connected
                print "\t Skipped. Too few vertexes (all kinds). "
            else:
                Root = VrtxCmpnts[i].index(SpecialGroup[i][0])  # always start MST from a special vertex 
#                Adj = Adjs[i]              # avoid creating new variable to speed up
#                Cmpnt = VrtxCmpnts[i]     # avoid creating new variable to speed up
                Num = len(Adjs[i])
                if Num > 1 : 
                    M = Prim(Adjs[i], Root)  # starting from the Root 
                    (Adj, W, Path, Degree, TreeNbr) = M.mst_prim(Adjs[i], [Root], i, [], M.degree, M.tree_nbr) # starting from the Root
                    #print len(Path)  # debugging line, Forrest 2011-05-25 20:33
#                    Seg = [[VrtxCmpnts[i][Idx] for Idx in Pair] for Pair in Path]  # The Idx is LOCAL (i.e., within the connected component) index of a vertex.
                    print len(Path), "links found."
                    
                    # pruning the MST Forrest 2011-09-24
                    Terminal, Branching =[], []
                    for Vrtx in xrange(0,Num):
                        if Degree[Vrtx] ==1:
                            Terminal.append(Vrtx)
                        elif Degree[Vrtx] > 2:
                            Branching.append(Vrtx)
                    
                    Special = [VrtxCmpnts[i].index(Vrtx) for Vrtx in SpecialGroup[i]] # converting global ID to local ID
                    Path = prune(Path, Degree, TreeNbr, Terminal, Branching, Special, VrtxCmpnts[i])
                    
                    for Pair in Path:                    
#                        (Src, Dst, Cost, CID) = Pair # commented, Forrest 2011-09-24
                        (Src, Dst) = Pair # Forrest 2011-09-24
#                        Segs.append( [VrtxCmpnts[i][Src], VrtxCmpnts[i][Dst], Cost, i] )
                        # pruning heuristic 1
#                        if Degree[Src] > 2 and Degree[Dst] ==1 and Cost < 4:
                            #Degree[Src] -= 1
#                        elif Degree[Dst] > 2 and Degree[Src] ==1  and Cost < 4:
                            #Degree[Dst] -= 1
                        # end of pruning heuristic 1
                        # pruning heuristic 2
#                        elif Degree[Src] ==1:
#                            OtherPairs, OtherTerminals = otherPaths(Dst, Pair, Degree[Dst], Path)
#                            Judge = [Cost < 1/3*OtherPair[2] for OtherPair in OtherPairs]
#                            Good = True
#                            for J in Judge:
#                                if J:  # if smaller do this. 
#                                    Good = False
#                                    Degree[Dst] -= 1
#                                    break
#                            if Good : # is no such a thing, do this:
#                                Segs.append( [VrtxCmpnts[i][Src], VrtxCmpnts[i][Dst], Cost, i] ) # here local indexes are converted into global indexes
#                        elif Degree[Dst] ==1:
#                            OtherPairs, OtherTerminals = otherPaths(Src, Pair, Degree[Src], Path)
#                            Judge = [Cost < 1/3*OtherPair[2] for OtherPair in OtherPairs]
#                            Good = True
#                            for J in Judge:
#                                if J:  # if smaller, do this
#                                    Good = False
#                                    Degree[Src] -= 1
#                                    break
#                            if Good: # is no such a thing, do this:
#                                Segs.append( [VrtxCmpnts[i][Src], VrtxCmpnts[i][Dst], Cost, i] ) # here local indexes are converted into global indexes
#                        # end of pruning heuristic 2
#                        # pruning heuristic 3
#                        elif Degree[Src] ==1:
#                            Results, Terminals = otherPaths(Src, Pair, Degree[Dst], Path)
#                            Good = True 
#                            for V in Terminals:
#                                if Degree[V] >1 :
#                                    Degree[Dst] -= 1
#                                    Good = False
#                                    break
#                            if Good: # is no such a thing, do this:
#                                Segs.append( [VrtxCmpnts[i][Src], VrtxCmpnts[i][Dst], Cost, i] ) # here local indexes are converted into global indexes
#                        elif Degree[Dst] ==1:
#                            Results, Terminals = otherPaths(Dst, Pair, Degree[Dst], Path)
#                            Good = True 
#                            for V in Terminals:
#                                if Degree[V] >1 :
#                                    Degree[Src] -= 1
#                                    Good = False
#                                    break
#                            if Good: # is no such a thing, do this:
#                                Segs.append( [VrtxCmpnts[i][Src], VrtxCmpnts[i][Dst], Cost, i] ) # here local indexes are converted into global indexes
                        # end of pruning heuristic 3
                        
#                        else: # no need to prune and thus put this segment into Segs directly                            
                            #Path = bfsReach([ VrtxCmpnts[i][Src], VrtxCmpnts[i][Dst], Cost, i], VrtxCmpnts[i], NbrLst) # deactivated forrest 2011-07-17
                        Path = [ VrtxCmpnts[i][Src], VrtxCmpnts[i][Dst] ] # Forrest 2011-07-17, do NOT map links onto the mesh
                        Segs.append(Path)
#                    Segs += Seg
#        print "number of segments: ", len(Segs)  
        return Segs

def getFundi(CurvFile, SurfFile, ToVTK=True, SurfFile2=''):
    '''includes rawFundi and fineFundi two functions, as well as file IO functions 
    
    This is a general framework for feature extraction
    
    Input
    ===== 
        CurvFile: string
                file name of the curvature-format file
                
        SurfFile: string
                file name of the surface-format file
                                
        ToVTK: Boolean
                Whether map the Sulci extraction result onto a surface and save in VTK format
                Default = True
 
    Example
    ========

    libvtx.getFeature('basin', CurvFile, SurfFile, ToVTK)
    '''
  
    print "Connecting clouds into fundus curves."
  
    Curvature = fileio.readCurv(CurvFile)
    CurvDisp = []
    for x in Curvature:
        if x > 0:
            CurvDisp.append(x)
        else:
            CurvDisp.append(0)
    
    Vrtx, Fc = fileio.readSurf(SurfFile)
    
    NbrLst = libbasin.vrtxNbrLst(len(Vrtx), Fc, SurfFile)
    
    Strict, Candidate = rawFundi(len(Vrtx), Curvature, NbrLst, Skip2 = False, Skip3 = False)  # turn Skip2 and Skip3 both on only when drawing curvature contours 

    VrtxCmpntFile = CurvFile + '.cmpnt.vrtx'  # need to run libbasin first to get components
    VrtxCmpnt = fileio.loadCmpnt(VrtxCmpntFile) # need to run libbasin first to get components
    Candidate
#    print len(Strict+Candidate)

# test mesh downsampling  ---

#    LeftVrtx, NewNbrLst = downsample([], VrtxCmpnt, NbrLst, Prob=0.5)
#    
#    Links = []
#    
#    for Center, Ring in enumerate(NewNbrLst): # Ring means
#        if Ring != []:
#            for Nbr in Ring: 
#                Links.append(str(Center) + " " + str(Nbr) + "\n")
#    
#    ReducedMesh = open(CurvFile + '.mesh.reduced.vtk', 'w')
#    libvtk.writeHeader(ReducedMesh)
#    libvtk.writePoint(ReducedMesh, Vrtx)
#    libvtk.writeSeg(ReducedMesh, Links)
#    ReducedMesh.close()
#    
#    exit(0)

# end of test mesh downsampling ---


# output clouds ---------------
    if Strict !=[]:
        FStrict = CurvFile + '.strict'
        fileio.writeList(FStrict, Strict)
        
    FCandidate = CurvFile + '.candidate'
    fileio.writeList(FCandidate, Candidate)
    
    if ToVTK:
        VTKFile = FCandidate + "." + SurfFile[-1*SurfFile[::-1].find('.'):] + '.vtk'
        libvtk.vrtxLst2VTK(VTKFile, SurfFile, FCandidate)
        if Strict != []:
            VTKFile = FStrict + "." + SurfFile[-1*SurfFile[::-1].find('.'):] + '.vtk'
            libvtk.vrtxLst2VTK(VTKFile, SurfFile, FStrict)
        if SurfFile2 != '':
            VTKFile = FCandidate + "." + SurfFile2[-1*SurfFile2[::-1].find('.'):] + '.vtk'
            libvtk.vrtxLst2VTK(VTKFile, SurfFile2, FCandidate)
            if Strict != []:  
                VTKFile = FStrict + "." + SurfFile2[-1*SurfFile2[::-1].find('.'):] + '.vtk'
                libvtk.vrtxLst2VTK(VTKFile, SurfFile2, FStrict)
# end of clouds output  --------- 

# get and output fundus curves -----------    

    Nodes = Strict # can be Candidates, Strict, Pits, Skeleton, etc.  

    Segs, NodeColor = lineUp(Nodes, NbrLst, VrtxCmpnt, Vrtx, CurvFile, Curvature)

    print "Saving fundus curved obtained via stringing up fundus vertexes into VTK files"

    FSeg = CurvFile + '.fundi.from.clouds'
    fileio.writeFundiSeg(FSeg, Segs)
    
    if ToVTK:
        VTKFile = FSeg + "." + SurfFile[-1*SurfFile[::-1].find('.'):] + '.vtk'
        libvtk.seg2VTK(VTKFile, SurfFile, FSeg, LUT=[CurvDisp], LUTname=['Curvature'])
        if SurfFile2 != '':
            VTKFile = FSeg + "." + SurfFile2[-1*SurfFile2[::-1].find('.'):] + '.vtk'
            libvtk.seg2VTK(VTKFile, SurfFile2, FSeg, LUT=[CurvDisp], LUTname=['Curvature'])
    
# end of get and output fundus curves ----    


# below are functions not in active use ---

def dfs4loop(Segs):
    '''The DFS search for check whether a graph is loop free.
    
    Step 1: Extract the graph as nodes and their childern from Segs
    Step 2: Run DSF on the graph
    Step 3: Report error if a visited node is visited again . Exit if DFS has visited all nodes. 
    
    
    Parameters
    ==========
    Segs    : list of 2-tuples of integers
        Each element is a 2-tuple of integers representing two vertexes. 
        
    Nodes    : list of integers
        Each element is a vertexes.
        
    Nbrs    : list of list of integers
        Each element is a list of Nbrs of of a node. Nodes are ordered by their vertex id. 
        Hence, Childern[i] is childern of the node of the i-th least vertex id. 
     
    '''
    
    pass
    
    

def loopfree():
    '''Check whether a graph is loop free 
    
    Do DFS to see whether you have traversed the 
    
    '''
  
    pass

def graphReduce(Cmpnts, Coordinates, NbrLst): 
    '''Reduce the edges in fundus face strip in one connected component for using Noah's minimal spanning tree algorithm later. 
        For each node, only a few edges left. A edge is kept by various criteria. Now we say the top 3 shortest edges.
        
        
    Parameters
    ===========
    
    NbrLst : list of list of integers
        neighbor list of vertexes 
        
    Coordinates : list of list (3-tuple) of floats
        Coordinates[i] is the coordinate of the i-th vertex on the pial surface
        
    Cmpnts : list of list integers
        each element is a list of IDs of vertexes that are in a component (could have been filtered by turnCmpnt)
        
    Dist : list of floats
        Distances from a vertex to all its neighbors, ordered as neighbors in NbrLst        
    
    Return 
    =======
    
    NewNbrLst : list of list of integers
        reduced neighbor list of vertexes
    
    Edges : list of list (2-tuple) of integers
        each element represents the IDs of two vertexes that form an edge
        Edges are left to be dumped into VTK files for visual inspection on the reduction result
           
    
    '''
    
    NewNbrLst = list(NbrLst)
    Edges = []
    
    for Cmpnt in Cmpnts:
        for Vrtx in Cmpnt:  
            Nbrs = NbrLst[Vrtx]
            Dist = {} 
            for Nbr in Nbrs:
                Dist[Nbr] = ( (Coordinates[Vrtx][0] - Coordinates[Nbr][0])**2 \
                            + (Coordinates[Vrtx][1] - Coordinates[Nbr][1])**2 \
                            + (Coordinates[Vrtx][2] - Coordinates[Nbr][2])**2)**0.5 
            Top = sorted(Dist.iteritems(), key=lambda (k,v) : (v,k))[:-1]  # the -3 here is a key, if we do a positive number, then only short edges are left and then the mesh is very discrete 
            NewNbrLst[Vrtx] = [Tuple[0] for Tuple in Top]
            [Edges.append([Vrtx, Tuple[0]]) for Tuple in Top]
                
    return NewNbrLst, Edges

def fillUp():
    '''Fill up spaces between vertex clouds 
    
    Seems no use now. 04/11/2011
    
    '''
    pass


def skeletonizeVrtx(FundiList, NbrLst, CNbr = 0.3, ANbr = 0.3, MaxIter = 10):
    '''Skeletonize vertex clouds, by iteratively removing vertexes that have a neighbor in the cloud but 'exposed' to the outer world.
    
    Input
    =====
        FundiList:    list of integers
            Each element is a vertex id. The vertex is considered as part of the raw fundi, as clouds. 
            
        NbrLst:    list of integers
            Each element is the id of vertexes neighboring with a vertex.
            
        MaxIter:    integer
            The maximum number of iteration allowed in Skeletonization by removing exposed cloud vertexes.
            
        CNbr :  
    
    Notes
    ======
    
    MaxIter = 3, 4 and 5 have no big difference after lineUp
    
    '''
    
#    ToBeRemoved = [False for i in xrange(0, len(FundiList))]
    Iteration = 0
    while Iteration <= MaxIter:
        ToBeRemoved = []
        Iteration += 1
    
        for i in xrange(0, len(FundiList)):
            Vrtx = FundiList[i]
            Nbrs = NbrLst[Vrtx]
    #        HasCloudNbr, HasAirNbr = False, False  # air neighbor means a neighbor that is not in the cloud
            CloudNbr, AirNbr = [], []# neighbors in the cloud, and neighbors not in the cloud 
            for Nbr in Nbrs:
                if Nbr in FundiList:
    #                HasiNbr = True
                    CloudNbr.append(1)
                    AirNbr.append(0)
                else:
    #                if Nbr in JustRemoved:
    #                    HasAirNbr = True
    #                    continue
    #                HasAirNbr = True
                    CloudNbr.append(0)
                    AirNbr.append(1)
    #        if HasCloudNbr and HasAirNbr:
#            if sum(CloudNbr) > len(Nbrs) * CNbr and sum(AirNbr) > len(Nbrs) * ANbr:
            if sum(CloudNbr) > 2 and sum(AirNbr) > 0:   
    #            ToBeRemoved[i] = True
                ToBeRemoved.append(Vrtx)
        
        if ToBeRemoved == []:
            break
        else: 
            for Rmv in ToBeRemoved:
                FundiList.remove(Rmv)
            
    return FundiList
