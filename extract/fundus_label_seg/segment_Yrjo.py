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
    
    SelfLabel: an integer
        the Label on the vertex 
        
    SecondLabel: an integer
        The 2nd nearest label to the vertex
        
    Distance: an integer
        the distance from the vertex to its 2nd nearest label 
     
    SecondVrtx: an integer
        the index of the vertex who offers the 2nd nearest label 
     
    ''' 

#    from collections import Counter
    

    SelfLabel = Labels[Vrtx]
    SecondLabel = -2
     
    Visited = [Vrtx]
    
    Fringe  = [Vrtx]
    Ring = 1
    while SecondLabel == -2:
        # step 1: check labels in fringe
        for V in Fringe:
            Visited.append(V)
            if Labels[V] > 0 and Labels[V]!= SelfLabel: 
            # added the AND condition to ensure two different labels, F 2012-02-21
                    SecondLabel = Labels[V]
                    Distance = Ring 
                    SecondVrtx = V
                    break
                
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
            return -1,-1,-1

    return SelfLabel, SecondLabel, Distance, SecondVrtx

def label_pits(Pits, NbrLst, Labels):
    '''Label Pits    
    
    Returns 
    ========
    
   Pit_Labels: a dictionary
        key is the ID of a pit; value are 2-tuples,
        which are the two nearest labels for the vertex,
        
    Pit_Dists: a dictionary 
        key is the ID of a pit; value are 2-tuples,
        which are distances from the vertex to two nearest labeled vertexes
        
    
    '''
#    Count = 0
    Pit_Labels, Pit_Dists = {}, {}
    SecondVrtxes={}
    for ID, Vrtx in enumerate(Pits):
        SelfLabel, SecondLabel, Distance, SecondVrtx = label_vertex(Vrtx, NbrLst, Labels)
        Pit_Labels[Vrtx] = [SelfLabel, SecondLabel]
        Pit_Dists[Vrtx]  = Distance
        SecondVrtxes[Vrtx] = SecondVrtx
        
#        print ID,"-th fundus vertex:", Vrtx
        
#        Count += 1 
#        if Count > 10:
#            break 
    
    return Pit_Labels, Pit_Dists, SecondVrtxes

def write_seg_TSV(Pit_Labels, Pit_Dists, SecondVrtxes, PitsFile, LabelFile):
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
        Fp.write(str(Vrtx) + '\t' + str(LabelPair[0]) + '\t' + str(LabelPair[1]) + '\t' + str(Pit_Dists[Vrtx]) + '\t' + str(SecondVrtxes[Vrtx]) + '\n')
    
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