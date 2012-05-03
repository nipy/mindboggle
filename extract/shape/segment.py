# assign two nearest labels (a label pair) to every fundus vertex
# Label pairs are aggregated to lowest number id of the regions in the aggregate
# The ending lines of this script does not traverse on a path. Loop in Shell.  

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

def label_pair_aggregate(LabelPair):
    '''Map a label pair to the lowest number id of the same region
    
    Parameters
    ============
    LabelPair: list of two integers
        A label pair to be mapped to aggregated label pair 
        
    Returns
    ==========
    
    no variable. Direct return. 

    Notes
    ======= 
        Label pair aggregation:
        precentral: [28,24]*,[3,24]*,[18,24]*
        postcentral: [22,29],[22,31]
        intraparietal: [29,31],[29,8]
        lateral occipital sulcus: [11,8]*, [11,29]*
        anterior occipital sulcus: [11,15]*,[11,9]
        circular sulcus: [35,30],[35,34],[35,12],[35,2],[35,24],[35,22],[35,31]
        cingulate sulcus: [2,14],[2,28],[2,17],[25,17]
        calcarine fissure: [13,25],[13,2]
        lateral H-shaped orbital sulcus: [12,18],[12,3]
        occipitotemporal sulcus: [7,9],[7,11]
        collateral sulcus: [7,6],[7,16],[7,13]
        interhemispheric fissure, dorsal margin: [17,28],[17,24],[17,22],[25,29],[5,29],[5,11]
       
    '''
    if LabelPair == [3, 24] or LabelPair == [18, 24]:
        return [24,28]
    elif LabelPair == [22, 31]:
        return [22, 29]
    elif LabelPair == [8, 29]:
        return [29, 31]
    elif LabelPair == [11, 29]:
        return [8, 11]
    elif LabelPair == [9, 11]:
        return [11, 15]
    elif LabelPair == [34, 35] or LabelPair == [12, 35] or LabelPair == [2, 35]\
     or LabelPair == [24, 35] or LabelPair == [22, 35] or LabelPair == [31, 35]:
        return [30, 35]
    elif LabelPair == [2, 28] or LabelPair == [2, 17] or LabelPair == [17, 25]:
        return [2, 14]
    elif LabelPair == [2, 13]:
        return [13, 25]
    elif LabelPair == [3, 12]:
        return [12, 18]
    elif LabelPair == [7, 11]:
        return [7, 9]
    elif LabelPair == [7, 16] or LabelPair == [7, 13]:
        return [7, 6]
    elif LabelPair == [17, 24] or LabelPair == [17, 22] or LabelPair == [25, 29]\
     or LabelPair == [5, 29] or LabelPair == [5, 11]:
        return [17, 28]
    else:
        return LabelPair

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
    
    TwoLabels: list of two integers
        the two most comment labels from Vrtx's nearest at least 2 labeled neighbors
        
    TwoDists: list of two floats
        the two Euclidean distances from Vrtx's nearest two labeled neighbors
        
    Notes
    ========
    
    The code used to assign label was different from what Arno expects. 
    Now it is updated to assign two nearest DIFFERENT labels to each fundus vertex. 
    - Forrest, 2012-02-21
    
    For labels on a ring, do we assign the closest or the most common label?
    - Forrest, 2012-02-21
    
    Always return the smaller label first. So (12,32) and (32,12) are the same, (12.32). 
     
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
            if Labels[V] > 0 and not Labels[V] in TwoLabels: 
            # added the AND condition to ensure two different labels, F 2012-02-21
                    TwoLabels.append(Labels[V])
                    TwoDists.append(Ring)
                    if len(TwoLabels) ==2 :
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
            return [-1, -1], [-1, -1]

    if TwoLabels[0] > TwoLabels[1]:
#        print TwoLabels, TwoLabels[::-1]
        return TwoLabels[::-1], TwoDists[::-1] # Always smaller label first 
    else:
        return TwoLabels, TwoDists
    
#    Top2Labels = Counter(TwoLabels).most_common(2)  
   
#    if len(Top2Labels) < 2:
#        return TwoLabels[:2], TwoDists[:2]
#    else: 
#        return [Top2Labels[0][0], Top2Labels[1][0]], [TwoDists[TwoLabels.index(Top2Labels[0][0])], TwoDists[TwoLabels.index(Top2Labels[1][0])]] 
    
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
#    Count = 0
    Fundus_Labels, Fundus_Dists = {}, {}  
    for ID, Vrtx in enumerate(Fundus_Vertexes):
        TwoLbs, TwoDists = label_vertex(Vrtx, NbrLst, Labels)
        Fundus_Labels[Vrtx] = label_pair_aggregate(TwoLbs)
        Fundus_Dists[Vrtx]   = TwoDists
        
#        print ID,"-th fundus vertex:", Vrtx
        
#        Count += 1 
#        if Count > 10:
#            break 
    
    return Fundus_Labels, Fundus_Dists

def write_seg(Fundus_Labels, Fundus_Dists, Num_Vrtx, FundiFile, FundiFile2=''):
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
    SegmentIdx = [-1 for i in xrange(Num_Vrtx)] # -1 means this not a fundus vertex
    
#    Pair_Index, Num_Diff_Labels = {}, 0
#    for Label_Pair in Fundus_Labels.values():
#        if not Label_Pair in Pair_Index.keys():
#            Pair_Index[(Label_Pair[0],Label_Pair[1])] = Num_Diff_Labels
#            Num_Diff_Labels += 1
    All_Pairs = []
    for LabelPair in Fundus_Labels.values():
        if not LabelPair[0]*100 + LabelPair[1] in All_Pairs:
            All_Pairs.append(LabelPair[0]*100 + LabelPair[1]) 
              
    for Vrtx, LabelPair in Fundus_Labels.iteritems():
        Segment[Vrtx] = LabelPair[0]*100 + LabelPair[1]# for absolute labels, use this line
        SegmentIdx[Vrtx] = All_Pairs.index(LabelPair[0]*100 + LabelPair[1]) # for better visualization, use this line
                    
    Distance = [-1 for i in xrange(Num_Vrtx)] # -1 means this not a fundus vertex
    for Vrtx, DistPair in Fundus_Dists.iteritems():
#        Distance[Vrtx] = DistPair[0] + DistPair[1]
        Distance[Vrtx] = DistPair[0]*100 + DistPair[1]
    
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
    SegFP.write('LOOKUP_TABLE default\n')
    
    for Seg in Segment: 
        SegFP.write(str(Seg)+'\n')
        
    SegFP.write('\nSCALARS SegIdx int \n')
    SegFP.write('LOOKUP_TABLE default\n')
    
    for Seg in SegmentIdx: 
        SegFP.write(str(Seg)+'\n')
#
#    SegFP.write('\nSCALARS To_Nearest int \n')
#    SegFP.write('LOOKUP_TABLE default\n')
#    
#    for Dist in Distance: 
#        SegFP.write(str(Dist)+'\n')
    
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
        SegFP.write('LOOKUP_TABLE default\n')
        
        for Seg in Segment: 
            SegFP.write(str(Seg)+'\n')

        SegFP.write('\nSCALARS SegIdx int \n')
        SegFP.write('LOOKUP_TABLE default\n')
    
        for Seg in SegmentIdx: 
            SegFP.write(str(Seg)+'\n')

#        SegFP.write('\nSCALARS To_Nearest int \n')
#        SegFP.write('LOOKUP_TABLE default\n')
#        
#        for Dist in Distance: 
#            SegFP.write(str(Dist)+'\n')
            
        SegFP.close()

    return 0
    
def seg_main(FundiFile, NbrLstFile, LabelFile):
    print "Segmenting fundi in", FundiFile 
    
    Fundus_Vertexes, Num_Vrtx, LUTs = load_fundus(FundiFile)
    NbrLst = load_nbr(NbrLstFile)
    Labels = load_labels(LabelFile)
    
    Fundus_Labels, Fundus_Dists = label_fundi(Fundus_Vertexes, NbrLst, Labels)
    write_seg(Fundus_Labels, Fundus_Dists, Num_Vrtx, FundiFile, FundiFile[:-4]+'.2nd.vtk')
#    write_seg(Fundus_Labels, Fundus_Dists, Num_Vrtx, FundiFile)  # do NOT write to inflated surface, unless visualization 

# a test run 
#seg_main('/forrest/data/MRI/MDD/50014/surf/lh.travel.depth.depth.fundi.vtk',\
#                             '/forrest/data/MRI/MDD/50014/surf/lh.vrtx.nbr',\
#                             '/forrest/data/MRI/MDD/50014/surf/lh.assign.pial.vtk')

# run under a path
import sys
#seg_main(sys.argv[1], sys.argv[2], sys.argv[3])

import os
Path = '/forrest/data/MRI/MDD/'
#Hemi = [sys.argv[1]]  # like lh or rh
for Hemi in ['lh','rh']:
    for DirIndx, Dir in enumerate(os.listdir(Path)):
        if len(Dir) == 5: # and DirIndx <= 30:
            print DirIndx, Dir, Hemi
            Prefix = Path+Dir+'/surf/'+Hemi
            seg_main(Prefix+'.euclidean.depth.depth.fundi.vtk', Prefix+'.vrtx.nbr', Prefix+'.assign.pial.vtk')
            
            