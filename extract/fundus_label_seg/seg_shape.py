# segmentation based on nearest labels for fundus vertexes, and output shape measures for grouped vertexes

def load_fundus(VTKFile):
    '''Load fundus vertex from a VTK file
    
    Notes 
    ==========
    
    My current fundi vtk file has the following LUTs in order:
    curvature, travel depth, thickness, sulc, fundus length
    
    '''
    import pyvtk
    Data = pyvtk.VtkData(VTKFile)
    Vertexes = Data.structure.points
    Fundus = Data.structure.lines
    LUTs = [Data.point_data.data[i].scalars for i in xrange(5)]
    
#    Depth =     VTKReader.point_data.data[1].scalars
    
#    print len(Fundus), "Fundus lines loaded"

    Fundus_Vertexes = []
    for Pair in Fundus:
        Fundus_Vertexes += Pair
    
    Fundus_Vertexes = list(set(Fundus_Vertexes)) # 

#    print len(Fundus_Vertexes), "Fundux vertexes loaded" 
        
    return Fundus_Vertexes, len(Vertexes), LUTs

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
    NbrLst = cPickle.load(Fp)
    Fp.close()
    return NbrLst

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
        the two most comment labels from Vrtx's nearest at least 2 labeled neighbors
        
    TwoLabels: list of two floats
        the two Euclidean distances from Vrtx's nearest two labeled neighbors
     
    ''' 

    from collections import Counter
    
    TwoLabels, TwoDists = [], []
    Visited = []
    
    if Labels[Vrtx] > 0:
        TwoLabels.append(Labels[Vrtx])
        TwoDists.append(0)
        Visited = [Vrtx]
    
    Fringe  = [Vrtx]
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
                if not Nbr in Visited and not Nbr in Nbrs:
                    Nbrs.append(Nbr)
                    
        Fringe = Nbrs
        Ring += 1

        if Ring > 500: # something may go wrong
            print "Taking too long to find nearest two labeled neighbors"
            print "Vrtx is:", Vrtx
            return [-1,-1], [-1,-1]
    
    Top2Labels = Counter(TwoLabels).most_common(2)  
   
    if len(Top2Labels) < 2:
        return TwoLabels[:2], TwoDists[:2]
    else: 
        return [Top2Labels[0][0], Top2Labels[1][0]], [TwoDists[TwoLabels.index(Top2Labels[0][0])], TwoDists[TwoLabels.index(Top2Labels[1][0])]] 
    
#    return TwoLabels[:2], TwoDists[:2] # let's say randomly pick the first two for fast prototyping 

def label_fundi(Fundus_Vertexes, NbrLst, Labels):
    '''Label many fundus vertexes
    
    
    Returns 
    ========
    
    Fundus_Labels: a dictionary
        key is the ID of a vertex on a fundus; value are 2-tuples,
        which are the two nearest labels for the vertex,
        
    Fundus_Dists: a dictionary 
        key is the ID of a vertex on fundu; value are 2-tuples,
        which are distances from the vertex to two nearest labeled vertexes
        
    
    '''
    Count = 0
    Fundus_Labels, Fundus_Dists = {}, {}  
    for ID, Vrtx in enumerate(Fundus_Vertexes):
        TwoLbs, TwoDists = label_vertex(Vrtx, NbrLst, Labels)
        Fundus_Labels[Vrtx] = TwoLbs
        Fundus_Dists[Vrtx]   = TwoDists
        
#        print ID,"-th fundus vertex:", Vrtx
        
        Count += 1 
        if Count > 10:
            break 
    
    return Fundus_Labels, Fundus_Dists

def seg_fundi(Fundus_Labels):
    '''Segment fundus vertexes based on their labels
    
    Parameters
    ============
    
    Fundus_Labels: a dictionary
        key is the ID of a vertex on a fundus; value are 2-tuples,
        which are the two nearest labels for the vertex,
        
    Fundus_Dist: a dictionary 
        key is the ID of a vertex on fundu; value are 2-tuples,
        which are distances from the vertex to two nearest labeled vertexes
    
    Groups: a dictionary 
        key is 4-digit label pair, value is a list of fundus vertexes that are assigned the label pair
        using two nearest labels scheme 
        
        
    Notes
    ======
    
    Label pairs are converted into 4-digit integer here. 
    
    '''
    
    Groups = {}
    for LabelPair in Fundus_Labels.values():
#        if not Pair in LabelPairs:
#            LabelPairs.append(Pair)
        Groups[LabelPair[0]*100 + LabelPair[1]] = []
    
    for Vrtx, LabelPair in Fundus_Labels.iteritems():
        Groups[LabelPair[0]*100 + LabelPair[1]].append(Vrtx)

    return Groups

def write_seg(Fundus_Labels, Fundus_Dists, Num_Vrtx, FundiFile, FundiFile2):
    '''Write the segmentation result into VTK
    
    To make things easier, this function copies everything from FundiFile and FundiFiles2 
    as the beginning part of outputting files 
    
    Parameters
    ===========
    
    Num_Vrtx: integer
        Number of vertexes on original mesh 
    
    Notes
    =======
    
    1. Please note the way we represent two labels, it's 4 digits. 
    
    2. Distances are not dumped yet. 
    
    '''

    Segment = [-1 for i in xrange(Num_Vrtx)] # -1 means this not a fundus vertex
    for Vrtx, LabelPair in Fundus_Labels.iteritems():
        Segment[Vrtx] = LabelPair[0]*100 + LabelPair[1]

    
    SegFile = FundiFile[:-4] + '.seg.vtk'  # drop .vtk suffix
    SegFP = open(SegFile, 'w')
    
    FundiFP = open(FundiFile,'r')  
    while True:
        Line = FundiFP.readline()
        if len(Line) < 1 : 
            break 
        else:
            SegFP.write(Line)
    FundiFP.close()        
    
    SegFP.write('\nSCALARS Segment int \n')
    SegFP.write('LOOKUP_TABEL default\n')
    
    for Seg in Segment: 
        SegFP.write(str(Seg)+'\n')
    
    SegFP.close()

    if FundiFile2 != '':
        SegFile = FundiFile2[:-8] + '.seg.2nd.vtk'  # drop .2nd.vtk suffix
        SegFP = open(SegFile, 'w')
        
        FundiFP = open(FundiFile2,'r')  
        while True:
            Line = FundiFP.readline()
            if len(Line) < 1 : 
                break 
            else:
                SegFP.write(Line)
        FundiFP.close()        
        
        SegFP.write('\nSCALARS Segment int \n')
        SegFP.write('LOOKUP_TABEL default\n')
        
        for Seg in Segment: 
            SegFP.write(str(Seg)+'\n')
        
        SegFP.close()

    return 0
    
def gen_shape(Groups, LUTs, FundiFile):
    '''Generate shape table for this hemisphere and write to TSV file
    
    Parameters
    ==========
    
    LUTs: list of lists of integers
        Each element is a per-vertex list of shape descriptor.
        The order is curvature, depth, thickness, sulc and fundus length  
     
    '''
    
    from numpy import mean
    
    TSV = FundiFile[:-4]+'.shape.tsv'
    FP = open(TSV, 'w')
    FP.write('label_pair \t curvature \t depth \t thickness \t convexity \t length \n')
    
    for LabelPair in Groups.keys():
        FP.write(str(LabelPair) + '\t')
        [FP.write(str(mean([LUT[j] for j in Groups[LabelPair]])) + ' \t ') for LUT in LUTs] 
        # i know the line above is horrible. - Forrest 2012-02-19
        FP.write('\n')

    FP.close()
    return 0

#FV= load_fundus('/forrest/data/MRI/MDD/50014/surf/lh.travel.depth.depth.fundi.vtk')

def seg_shape_main(FundiFile, NbrLstFile, LabelFile):
    
    Fundus_Vertexes, Num_Vrtx, LUTs = load_fundus(FundiFile)
    NbrLst = load_nbr(NbrLstFile)
    Labels = load_labels(LabelFile)
    
    Fundus_Labels, Fundus_Dists = label_fundi(Fundus_Vertexes, NbrLst, Labels)
    Groups = seg_fundi(Fundus_Labels)
    write_seg(Fundus_Labels, Fundus_Dists, Num_Vrtx, FundiFile, FundiFile[:-4]+'.2nd.vtk')
    gen_shape(Groups, LUTs, FundiFile)

seg_shape_main('/forrest/data/MRI/MDD/50014/surf/lh.travel.depth.depth.fundi.vtk',\
                             '/forrest/data/MRI/MDD/50014/surf/lh.vrtx.nbr',\
                             '/forrest/data/MRI/MDD/50014/surf/lh.assign.pial.vtk')

 