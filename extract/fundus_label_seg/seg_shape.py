# segmentation based on nearest labels for fundus vertexes, and output shape measures for grouped vertexes

def load_fundus(VTKFile):
    '''Load fundus vertex from a VTK file
    '''
    import pyvtk
    Data = pyvtk.VtkData(VTKFile)
    Fundus = Data.structure.lines
    print len(Fundus), "Fundus lines loaded"
#    Fundus = Fundus[0]
    Fundus_Vertexes = []
    for Pair in Fundus:
 #       print Pair
        Fundus_Vertexes += Pair
    
    Fundus_Vertexes = list(set(Fundus_Vertexes)) # 

    print len(Fundus_Vertexes), "Fundux vertexes loaded" 
        
    return Fundus_Vertexes

def load_labels(VTKFile):
    '''Load voted labels of a mesh 
    '''
    import pyvtk
    Labels = pyvtk.VtkData(VTKFile).point_data.data[0].scalars
    return Labels

def load_nbr(PickleFile):
    '''Load the Vertex Neighbor list of a mesh
    '''
    import cPickle
    Fp = open(PickleFile, 'r')
    Nbr = cPickle.load(Fp)
    Fp.close()
    return Nbr

def label_vertex(Vrtx, NbrLst, Labels):
    '''Give a vertex two labels from its nearest labeled neighbors, using BFS
    
    Parameters
    ===========
    
    Vrtx : integer
        The ID of a Vrtx
        
    NbrLst: List of lists of integers
        the i-th element is the neighbors of vertex i 
        
    Labels: list of integers
        Labels assigned for vertexes from 21 atlases, unlabeled vertexes have value -1, and 1-35 otherwise.
        
    Visited: list of integers
        IDs of vertexes whose labels have been checked
        
    Ring: integer
        the hop from Vrtx to the search fringe
        
    Nbrs: list of integers
        all neighbors on a ring to Vrtx
    
    Returns
    ========
    
    TwoLabels: list of two integers
        the two labels from Vrtx's nearest two labeled neighbors
        
    TwoLabels: list of two floats
        the two Euclidean distances from Vrtx's nearest two labeled neighbors
     
    ''' 
    
    TwoLabels, TwoDists = [], []
    
    if Labels[Vrtx] > 0:
        TwoLabels.append(Labels[Vrtx])
        TwoDists.append(0.0)
    
    Visited = [Vrtx]
    Fringe  = NbrLst[Vrtx]
    Ring = 1
    while len(TwoLabels) < 2:
        # step 1: check labels in fringe
        for V in Fringe:
            Visited.append(V)
            if Labels[V] > 0:
                    TwoLabels.append(Labels[V])
                    TwoDists.append(Ring)
        
        
        # step 2: search fringe propagation
        Nbrs = []
        for V in Fringe:
            for Nbr in NbrLst[V]:
                if not Nbr in Visited:
                    Nbrs.append(Nbr)
                    
        Fringe = Nbrs
        Ring += 1

    if Ring > 100: # something may go wrong
            print "Taking too long to find nearest two neighbors"
            return [-1,-1], [-1,-1]

    return TwoLabels[:2], TwoDists[:2] # let's say randomly pick the first two for fast prototyping 

def label_fundi(FundiFile, NbrLstFile, LabelFile):
    '''Label many fundus vertexes
    '''
    Fundus_Vertexes = load_fundus(FundiFile)
    NbrLst = load_nbr(NbrLstFile)
    Labels = load_labels(LabelFile)
    
    Labeled, Dists = [], [] 
    for Vrtx in Fundus_Vertexes:
        TwoLbs, TwoDists = label_vertex(Vrtx, NbrLst, Labels)
        Labeled.append(TwoLbs), Dists.append(Dists) 
        break
    
    return Labeled, Dists

#FV= load_fundus('/forrest/data/MRI/MDD/50014/surf/lh.travel.depth.depth.fundi.vtk')

Labeled, Dists = label_fundi('/forrest/data/MRI/MDD/50014/surf/lh.travel.depth.depth.fundi.vtk',\
                             '/forrest/data/MRI/MDD/50014/surf/lh.vrtx.nbr',\
                             '/forrest/data/MRI/MDD/50014/surf/lh.assign.pial.vtk')



 