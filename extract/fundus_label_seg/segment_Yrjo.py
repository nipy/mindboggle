# segmentation based on nearest labels for fundus vertexes

def load_yrjo_pits(File):
    '''Load a column of test file and return the list
    '''
    Fp = open(File, 'r')
    Pits = [] 
    for Line in Fp.readlines():
        Pits.append(int(Line.split()[0]) - 1 )  # Yrjo's index starts from 1 
    Fp.close()

    return Pits

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
        The ID of a vertex
        
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
           
    FourLabels: list of 4 integers
        The 1st, 2nd, 3rd and 4th nearest labels
                
    Distances: list of 3 integers
        the distance from the vertex to its 2nd, 3rd, and 4th nearest labels 
     
    LabelVrtxes: list of 3 integers
        the indexes of the vertexes who offer the 2nd, 3rd and 4th nearest labels 
        
    Notes
    =======
    
    The labels from VTK file are MATLAB-indexed. I need to convert it to Python-indexed. 
    Stick with the labels on my file. 
     
    ''' 

#    from collections import Counter
    

    FourLabels = []
    Distances = []
    LabelVrtxes = []
     
    Visited = [Vrtx]
    
    Fringe  = [Vrtx]
    Ring = 0
    while len(FourLabels) < 4:
        # step 1: check labels in fringe
        for V in Fringe:
            Visited.append(V)
            if Labels[V]-1 > 0 and not Labels[V]-1 in FourLabels: 
            # added the AND condition to ensure two different labels, F 2012-02-21
                    FourLabels.append(Labels[V]-1)
                    Distances.append(Ring) 
                    LabelVrtxes.append(V)
                
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
            return -1, [], [], []

    return FourLabels, Distances, LabelVrtxes

def label_pits(Pits, NbrLst, Labels):
    '''Label Pits    
    
    Returns 
    ========
    
   Pit_Labels: a dictionary
        key is the ID of a pit; value are 4-tuples,
        the 1st, 2nd, 3rd and 4th labels
        
    Pit_Dists: a dictionary 
        key is the ID of a pit; value are 4-tuples,
        which are distances from the vertex to 1st, 2nd, 3rd and 4th nearest labeled vertexes
        
    
    '''
#    Count = 0
    Pit_Labels, Pit_Dists = {}, {}
    Pit_Vrtxes={}
    for ID, Vrtx in enumerate(Pits):
        FourLabels, Distances, LabelVrtxes = label_vertex(Vrtx, NbrLst, Labels)
        Pit_Labels[Vrtx] = FourLabels
        Pit_Dists[Vrtx]  = Distances[1:]  # because the 1st one is 0 
        Pit_Vrtxes[Vrtx] = LabelVrtxes[1:] # because the 1st one is itself
    
    return Pit_Labels, Pit_Dists, Pit_Vrtxes

def write_seg_TSV(Pit_Labels, Pit_Dists, Pit_Vrtxes, PitsFile, LabelFile):
    '''Dump the results into TSV to give to Yrjo 
    
    
    Returns
    ==========
    SeftLabels: list of integers
        Each element is the label given to the pit
        
    SecondLabels: list of integers
        Each element is the 2nd (different) nearest label to the pit
        
    Dist: list of integers
        Each element is the distance from the pit to its 2nd nearest label 
    
    '''
    
    import pyvtk
    Data = pyvtk.VtkData(LabelFile)
    Points = Data.structure.points
    
    Num_Vrtx = len(Points)
    
    Fp = open(PitsFile + '.seg.tsv','w')
    
    for Vrtx, LabelPair in Pit_Labels.iteritems():
        Line  = str(Vrtx) + '\t' 
        for i in range(4):
            Line += (str(LabelPair[i]) + '\t' )
        for i in range(3):
            Line += (str(Pit_Dists[Vrtx][i]) + '\t' )
        for i in range(3):
            Line += (str(Pit_Vrtxes[Vrtx][i]) + '\t' )
            
        Fp.write(Line + '\n')
    
    Fp.close()
  
def seg_main(PitsFile, NbrLstFile, LabelFile):
    
    Pits = load_yrjo_pits(PitsFile)
    NbrLst = load_nbr(NbrLstFile)
    Labels = load_labels(LabelFile)
    
    Fundus_Labels, Fundus_Dists, SecondVrtxes = label_pits(Pits, NbrLst, Labels)
    #write_seg_VTK(Fundus_Labels, Fundus_Dists, PitsFile, LabelFile[:-9]+'.inflated.vtk', LabelFile)
    write_seg_TSV(Fundus_Labels, Fundus_Dists, SecondVrtxes, PitsFile, LabelFile)


#seg_main('/forrest/data/MRI/Kirby/KKI2009-11/surf/KKI2009-11_LH.txt',\
#                             '/forrest/data/MRI/Kirby/KKI2009-11/surf/lh.vrtx.nbr',\
#                             '/forrest/data/MRI//Kirby/KKI2009-11/label/lh.aparcNMMjt.pial.vtk')

import sys
seg_main(sys.argv[1], sys.argv[2], sys.argv[3]) 

#def write_seg_VTK(Pit_Labels, Pit_Dists, FundiFile, LabelFile2, LabelFile):
#    '''Obsolete function. DO NOT USE!    
#    
#    Write the segmentation result into VTK, for visualization, not giving to Yrjo
#    
#    To make things easier, this function copies everything from FundiFile and FundiFiles2 
#    as the beginning part of outputting files 
#    
#    Parameters
#    ===========
#    
#    Num_Vrtx: integer
#        Number of vertexes on original mesh 
#    
#    Notes
#    =======
#    
#    Please note the way we represent two labels, it's 4 digits. 
#
#     
#    
#    '''
#
#    import pyvtk
#    Data = pyvtk.VtkData(LabelFile)
#    Points = Data.structure.points
#    
#    Num_Vrtx = len(Points)
#
#
#    Segment = [-1 for i in xrange(Num_Vrtx)] # -1 means this not a fundus vertex
#    SegmentIdx = [-1 for i in xrange(Num_Vrtx)] # -1 means this not a fundus vertex
#    
##    Pair_Index, Num_Diff_Labels = {}, 0
##    for Label_Pair in Fundus_Labels.values():
##        if not Label_Pair in Pair_Index.keys():
##            Pair_Index[(Label_Pair[0],Label_Pair[1])] = Num_Diff_Labels
##            Num_Diff_Labels += 1
#    All_Pairs = []
#    for LabelPair in Pit_Labels.values():
#        if not LabelPair[0]*100 + LabelPair[1] in All_Pairs:
#            All_Pairs.append(LabelPair[0]*100 + LabelPair[1]) 
#              
#    for Vrtx, LabelPair in Pit_Labels.iteritems():
#        Segment[Vrtx] = LabelPair[0]*100 + LabelPair[1]# for absolute labels, use this line
#        SegmentIdx[Vrtx] = All_Pairs.index(LabelPair[0]*100 + LabelPair[1]) # for better visualization, use this line
#                
#    Distance = [-1 for i in xrange(Num_Vrtx)] # -1 means this not a fundus vertex
#    for Vrtx, DistPair in Pit_Dists.iteritems():
##        Distance[Vrtx] = DistPair[0] + DistPair[1]
#        Distance[Vrtx] = DistPair[0]*100 + DistPair[1]
#    
#    SegFile = FundiFile + '.seg.vtk'  # drop .fundi.vtk suffix
##    SegFP = open(SegFile, 'w')
#
#    Pits = [Vrtx for Vrtx in Pit_Labels.keys()]
#    pyvtk.VtkData(pyvtk.PolyData(points=Points, vertices=Pits),\
#                  pyvtk.PointData(pyvtk.Scalars(Segment, name='Segment'),\
#                                  pyvtk.Scalars(SegmentIdx, name='SegIdx'),\
#                                  pyvtk.Scalars(Distance, name='To_Nearest_Two'))).\
#                  tofile(SegFile, 'ascii')
#
#    if LabelFile2 != '':
#        SegFile = FundiFile + '.seg.2nd.vtk'  # drop .2nd.vtk suffix
#    
#        Data = pyvtk.VtkData(LabelFile2)
#        Points = Data.structure.points        
#        
#        Pits = [Vrtx for Vrtx in Pit_Labels.keys()]
#    pyvtk.VtkData(pyvtk.PolyData(points=Points, vertices=Pits),\
#                  pyvtk.PointData(pyvtk.Scalars(Segment, name='Segment'),\
#                                  pyvtk.Scalars(SegmentIdx, name='SegIdx'),\
#                                  pyvtk.Scalars(Distance, name='To_Nearest_Two'))).\
#                  tofile(SegFile, 'ascii')
#
#    return 0