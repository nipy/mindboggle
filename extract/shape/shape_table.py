# this script loads segmented fundi 42 brains and filter label pairs that do not appear in all braings
# then, it computes the shape tables, which consist of mean values of shape measures of largest 
# components for each label pairs in each brain
# The load_multi_segs takes really long to run. So use ending lines to pickle the results

def load_multi_segs(Path, Suffix, Quit=1000):
    '''Load segmented fundi from from 52x2=104 hemisphere hemispheres     
    
    
    Parameters
    ==============
    
    Path: string
        a path. The Path is highly specified for our current file hierarchy, 
        thus    Path/subject/surf/....
   
    Suffix: string
        The suffix of the VTK file that contains label pairs assigned to fundus vertexes, 
        e.g., '.euclidean.depth.depth.fundi.seg.vtk' 
   
    LUT: all SCALARS loaded from a segmented fundi file
        The order of SCALARS is:
        1. curvature    0
        2. depth        1
        3. thickness    2
        4. convexity    3
        5. funduslength 4
        6. fundusID,    5
        7. Segment, the label pair  6 
        8. Segment indexes in the hemisphere 7   -- no need to keep, will be dropped in upper stream
        9. distances to two nearest labeled vertexes  8   -- not needed by Yrjo, will be dropped in upper stream  
    
    LabelPair: dictionary 
        Keys are vertex IDs; values are label pair as a 4-digit integer  
        
    Fundus_Vertex: list of integers
        Fundus vertexes in one hemisphere
        
    Fundus_Line: list of 2-tuples of integers
        Fundus lines in one hemisphere
    
    Returns
    ==========
    
    LabelPair_All: list of dictionaries
        Each element is a LabelPair 
        in one hemisphere
        The first 52 ones will be for left hemispheres 
        and the last 52 for right hemispheres
     
    LUT_All: list of LUT's
    
    Fundus_Vertex_All: list of lists of integers
        Each element is a *Fundus_Vertex*
        Left hemisphere: first 52. 
        
    Fundus_Line_All: list of lists of 2-tuples of integers
        Each element is a *Fundus_Line*
    
    Vertexes_All: list of list of 3-tuples of floats (obsolete)
        Each element is the coordinates of all vertexes in one hemisphere
    
    Notes
    ======
    
    Variables of suffix All are for all hemispheres, the list of all corresponding Singular variables 
    
    PointData id may need to changed once upper stream code changes, e.g., dropping fundusLength and fundusID 
    
    '''
    import os
    import pyvtk
    
    LabelPair_All, LUT_All, Fundus_Vertex_All, Fundus_Line_All = [], [], [], []
    
    Count = 0 # for debug 
    
    for Hemi in ['/surf/lh','/surf/rh']:
        for Subj in os.listdir(Path):
            if len(Subj) == 5: # for MDD subjects
                VTKFile = Path + Subj + Hemi + Suffix
                print VTKFile
                Data = pyvtk.VtkData(VTKFile)
                #Points = Data.structure.points
                Fundus_Line = Data.structure.lines
                
                LUT = [Data.point_data.data[i].scalars for i in [0,1,2,3]]
                LUT_All.append(LUT)
                
                LabelPair = {} #[-1 for i in xrange(len(Vertexes))]
                for Vrtx, Label in enumerate(Data.point_data.data[6].scalars):  #  
                    if Label != -1: # -1 * 100 + (-1) = -101
                        LabelPair[Vrtx] = Label

                Fundus_Vertex = []
                for Pair in Fundus_Line:
                    Fundus_Vertex += Pair
                Fundus_Vertex = list(set(Fundus_Vertex)) # 
                Fundus_Vertex_All.append(Fundus_Vertex)
                Fundus_Line_All.append(Fundus_Line)
                
                LabelPair_All.append(LabelPair)
                                
                Count += 1  # for fast debugging
            
            if Count == Quit:   # for fast debugging 
                return LabelPair_All, LUT_All, Fundus_Vertex_All, Fundus_Line_All
    
    return LabelPair_All, LUT_All, Fundus_Vertex_All, Fundus_Line_All

def common_pairs(LabelPair_All):
    '''Load segmented fundi  and leave only those appear in all hemispheres
    
    Parameters
    ==============

    LabelPair: dictionary 
        Keys are vertex IDs; values are label pair as a 4-digit integer  
    
     LabelPair_All: list of dictionaries
        Each element is a LabelPair in one hemisphere
        The first 52 ones will be for left hemispheres 
        and the last 52 for right hemispheres    
    
    Partial: list of integers
        Label pairs (as integers) that only exist in some hemispheres
        
    Pair: 4-digit integer
        

    Return
    =========
    
    All: list of integers
       Label pairs (as integers) that exist in all hemispheres     
    
    ''' 
    
    All  = []
    for LabelPair in LabelPair_All:
        All += LabelPair.values()
    All = list(set(All))
    
    Partial = []
    
    for Pair in All: # now find out all pairs that 
        for LabelPair in LabelPair_All:
            if not Pair in LabelPair.values() and not Pair in Partial:
                Partial.append(Pair)

    for Pair in Partial:
        All.remove(Pair)

    return All

def the_other(One, Tuple):
    '''Given one element of a tuple, return the other element. Return -1, if One is not in Tuple
    
    Parameters
    Tuple: list of two integers
    
    One: integer
     
    '''
    
    if One == Tuple[0]:
        return Tuple[1]
    elif One == Tuple[1]:
        return Tuple[0]
    else:
        return -1
    
def span(Seed, Fundus_Line, Vertex_this_Group, Vertexes):
    '''From Seed, tracking along Fundus_Lines, return all contiguous vertexes of the same Label and the length of them
    
    Parameters 
    ==============
    
    Seed: integer
        a vertex from which to start find the all contiguous vertexes of the same label 

    Fundus_Line: list of 2-tuples of integers
        Each element is the two vertexes that form a fundus line segment 
        
    Vertex_this_Group: list of integers
        All vertexes that have the same label pair as Seed

    Vertexes: list of floats
        Coordinates of all vertexes 
        
    Extended: Boolean
        Whether a new round of search has found new fundus vertexes of the same label pair
    
    Returns
    =========
    
    Spanned: list of integers
        All vertexes can be visited via Fundus_Line, from Seed
        
    Length: integer (will be float soon)
        Length of this connected component of same label pairs
    
    Notes
    ======
    
    The implementation has high time-complexity but since the problem size is small, it is okay.
    
    Current length is number of hops, not mesh length 
    
    '''
    
    from math import sqrt
    
    Spanned = [Seed]
    Extended = True
    Length = 0
    
    while Extended:
        Extended = False
        for Vrtx in Spanned:
            for Line in Fundus_Line:
                Another = the_other(Vrtx, Line)
                if Another > -1:
                    if Another in Vertex_this_Group and not Another in Spanned:
                        Spanned.append(Another)
                        Extended = True
#                        Length += 1
                        Length += sqrt( (Vertexes[Vrtx][0] - Vertexes[Another][0])**2 + \
                                        (Vertexes[Vrtx][1] - Vertexes[Another][1])**2 + \
                                        (Vertexes[Vrtx][2] - Vertexes[Another][2])**2  ) 
                    
    return Spanned, Length

def largest_groups(Common, Fundus_Line, Fundus_Vertex, LabelPair, Vertexes):
    '''Find largest contiguous groups of vertexes of all common label pairs in a hemisphere
    
    Parameters
    ===========
    
    Common: list of integers
        Label pairs that are common in all hemispheres
        
    Pair: 4-digit integer
        A label pair 
        
    FundiLines: list of 2-tuples
        Each element is a line segment on fundi within a hemisphere
        
    LabelPair: dictionary
        Keys are vertex IDs; values are label pair as a 4-digit integer  

    Vertexes: list of floats
        Coordinates of all vertexes 

    Groups_of_this_Pair: list of lists of integers
        Each element is vertexes consisting of a connected component of the label pair *Pair*
      
    Length_of_this_Pair: list of integers (will be floats soon)
        Length of different connected components of the same label pair
      
    Returns
    ========
    
    Largest_Groups: dictionary
        Key is a label pair (as 4-digit integer); value is vertexes of the largest connect component of the label pair
        
    Length_LabelPairs: list of integers (will be floats soon) 
        Length of largest connected component of each label pair. 
        This is a shape measure to be included in shape table
          
    '''
    Largest_Groups = {}
    Length_LabelPairs = {}
    
#    print "Common label pairs", len(Common)
    
    for Pair in Common:  # for every label pair that is common across hemispheres
        Vertexes_of_this_Group = [Vrtx for Vrtx in Fundus_Vertex if LabelPair[Vrtx] == Pair]
        Vertexes_of_this_Group_yet_Visited =list(Vertexes_of_this_Group)
        Groups_of_this_Pair, Length_of_this_Pair = [], []
        while len(Vertexes_of_this_Group_yet_Visited) > 0:
            Spanned, Length = span(Vertexes_of_this_Group_yet_Visited[0], Fundus_Line, Vertexes_of_this_Group, Vertexes)
            Groups_of_this_Pair.append(Spanned)
            Length_of_this_Pair.append(Length)
            
            if Spanned == Vertexes_of_this_Group: # hope this can be faster without removing vertexes from Vertexes_of_this_Group_yet_Visited 
                 break
            else:
                for Vrtx in Spanned:
                    Vertexes_of_this_Group_yet_Visited.remove(Vrtx)

        Largest_Index = Length_of_this_Pair.index(max(Length_of_this_Pair))
        Largest_Groups[Pair] = Groups_of_this_Pair[Largest_Index]
        Length_LabelPairs[Pair] = Length_of_this_Pair[Largest_Index] 
                        
#        print "label pair", Pair, "size(s):", map(len, Groups_of_this_Pair)#, Largest_Index
        print "label pair", Pair, "size(s):", Length_of_this_Pair
         
    return Largest_Groups, Length_LabelPairs

def gen_shape(Groups, LUTs, AvgFile, IdvFile,  Length):
    '''Generate two kinds of shape tables for this hemisphere and write to TSV files
    
    Parameters
    ==========
    
    LUTs: list of lists of integers
        Each element is a per-vertex list of shape descriptor.
        The order is curvature, depth, thickness, sulc and fundus length  
     
    Length: dictionary 
        Keys are label pairs (as 4-digit integer) and values are lengths of the largest component of corresponding label pairs
     
    AvgFile: string
        File path to save averages of shape measures 
        
    IdvFile: string
        File path to save shape measures of every vertex on longest fundus component of a label pair
    
    Notes
    ========
    
    We generate two kinds of shape tables here. One for averages of shape measures. One for all fundus vertexes.  
    
    '''
    
    from numpy import mean
    
    FP = open(AvgFile, 'w')
    FP.write('label_pair \t curvature \t depth \t thickness \t convexity \t length \n')
    
    for LabelPair in Groups.keys():
        FP.write(str(LabelPair) + '\t')
        [FP.write(str(mean([LUT[j] for j in Groups[LabelPair]])) + ' \t ') for LUT in LUTs] 
        # i know the line above is horrible. - Forrest 2012-02-19
        FP.write(str(Length[LabelPair]))
        FP.write('\n')

    FP.close()
    
    FP = open(IdvFile, 'w')
    FP.write('label_pair \t curvature \t depth \t thickness \t convexity \t length \n')
    
    for LabelPair, Fundus_Vertexes in Groups.iteritems():
        for Vertex in Fundus_Vertexes:
            FP.write(str(LabelPair) + '\t')
            [FP.write(str(LUT[Vertex]) + ' \t ') for LUT in LUTs] 
            FP.write(str(Length[LabelPair]))
            FP.write('\n')

    FP.close()
    
    return 0

def gen_shape_all(Path, Common, LabelPair_All, LUT_All, Fundus_Vertex_All, Fundus_Line_All):
    '''Generate shape tables for all hemispheres 
    
    Parameters
    ============
    
    LabelPair_All: list of dictionaries
        Each element is a LabelPair 
        in one hemisphere
        The first 52 ones will be for left hemispheres 
        and the last 52 for right hemispheres
     
    LUT_All: list of LUT's
    
    Fundus_Vertex_All: list of lists of integers
        Each element is a Fundus_Vertex
        Left hemisphere: first 52. 
        
    Fundus_Line_All: list of lists of 2-tuples of integers
        Each element is a Fundus_Line
        
    Lengths: dictionary 
        Keys are label pairs (as 4-digit integer) and values are lengths of the largest component of corresponding label pairs
    
    Returns
    ========
    
    None, written to file by gen_shape()
    
    '''
    Count = 0
     
    import os
    import pyvtk
    
    for Hemi in ['/surf/lh','/surf/rh']:
        for Subj in os.listdir(Path):
            if len(Subj) == 5: # for MDD subjects
                print Subj + Hemi[len(Hemi)-Hemi[::-1].find('/'):]
                AvgFile = './Shape_Travel/' + Subj + Hemi[len(Hemi)-Hemi[::-1].find('/'):] + '.avg.tsv'
                IdvFile = './Shape_Travel/' + Subj + Hemi[len(Hemi)-Hemi[::-1].find('/'):] + '.idv.tsv'
                
                # load Vertexes from here
                VTKFile = Path + Subj + Hemi + '.pial.vtk'
                print "Loading coordinates from", VTKFile
                Data = pyvtk.VtkData(VTKFile)
                Vertexes = Data.structure.points
                # end load Vertexes from here
                
                Group, Lengths = largest_groups(Common, Fundus_Line_All[Count], Fundus_Vertex_All[Count], LabelPair_All[Count], Vertexes)
                gen_shape(Group, LUT_All[Count], AvgFile, IdvFile, Lengths)
                Count += 1
                
#                if Count == 2:  # for quick debugging 
#                    return 0
                 
    return 0
# now something real 

import cPickle

#S,L,FV,FL = load_multi_segs('/forrest/data/MRI/MDD/', '.travel.depth.depth.fundi.seg.vtk', Quit= 500)
#with open('labels.pickle', 'w') as f:
#    cPickle.dump([S,L,FV,FL], f)

with open('labels.pickle') as f:
    S,L,FV,FL = cPickle.load(f)

Common = common_pairs(S)

gen_shape_all('/forrest/data/MRI/MDD/',Common, S, L, FV, FL)

#Test_Groups, TestLength = largest_groups(Common, FL[0], FV[0], S[0])
#gen_shape(Test_Groups, L[0], 'test_shape.tsv', TestLength) 
