# This scripts generates 
# components for each label pairs in each brain
# The load_multi_segs takes really long to run. So use ending lines to pickle the results

import sys

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
                AvgFile = './Shape_Convexity/' + Subj + Hemi[len(Hemi)-Hemi[::-1].find('/'):] + '.avg.tsv'
                IdvFile = './Shape_Convexity/' + Subj + Hemi[len(Hemi)-Hemi[::-1].find('/'):] + '.idv.tsv'
                
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


def input_handle():
    """Process inputs from Shell and prepare maps and features for the core algorithms

    Outputs
    ===========
        Fundi : list of integers
            IDs of vertexes that form fundi, provided by ?h.pial.fundi.vtk
        
        Sulci : list of integers
            IDs of vertexes that form sulci (Forrest's depth-thresholded sulci or Yrjo's hole-free ones)
            provided by ?h.pial.sulci.vtk

        Region_Labels : list of integers
            Each element corresponds to a mesh vertex, provided by ?h.pial.labels.max.vtk
            The labels are as those used by FreeSurfer. 
            
        Sulcus_Labels : list of integers
            Each element corresponds to a mesh vertex, provided by ?h.pial.{fundi,sulci}.segmented.vtk

        Depth : list of floats
            Each element represents the travel depth of the corresponding mesh vertex, 
            provided by ?h.pial.depth.vtk

        Curv_mean: list of floats
            Each element represents the mean curvature of the corresponding mesh vertex, 
            provided by ?h.pial.curvature.mean.vtk

        Curv_max : list of floats
            Each element represents the maximum curvature of the corresponding mesh vertex, 
            provided by ?h.pial.curvature.max.vtk
 
        Curv_min : list of floats
            Each element represents the minimum curvature of the corresponding mesh vertex, 
            provided by ?h.pial.curvature.min.vtk   

        Curv_gauss : list of floats
            Each element represents the Gaussian curvature of the corresponding mesh vertex, 
            provided by ?h.pial.curvature.gauss.vtk
       
        Convexity : list of floats
            Each element represents the convexity of the corresponding mesh vertex, 
            provided by ?h.pial.sulci.vtk
 
    """
    import getopt
    import sys
    def print_help(): # a nested function
        print "\n  Usage: python shape_table.py  [OPTIONS]\n"
        print "  Options: "
        print "  --labelSurf labelSurfVTK"
        print "    The VTK file that contains a labeled (using region/gyral labels) surface"
        print "  --segSulci segSulciVTK"
        print "    The VTK file that contains labeled (using sulcus labels) sulci"
        print "  --segFundi segFundiVTK"
        print "    The VTK file that contains labeled (using sulcus labels) fundi"
        print "  --depth DepthVTK"
        print "    the VTK file containing travel depth map\n"
        print "  --curvmin curvminVTK"
        print "    the VTK file containing minimum curvature map\n"
        print "  --curvmax curvmaxVTK"
        print "    the VTK file containing maximum curvature map\n"
        print "  --curvmean curvmeanVTK"
        print "    the VTK file containing mean curvature map\n"
        print "  --curvgauss curvgaussVTK"
        print "    the VTK file containing Gaussian curvature map\n"
        print "  --thick ThicknessFile"
        print "    the thickness file provided by FreeSurfer\n"
        print "  --convex ConvexityFile"
        print "    the convexity file provided by FreeSurfer\n"
        print "  --patchtable PatchTableFile"
        print "    the file to save shape table for label surface patches"
        print "  --funditable FundiTableFile"
        print "    the file to save shape table for segmented fundi"
        print "  --sulcitable"
        print "    the file to save shape table for segmented sulci"
        print "  Examples: "
#        print "    python shape_table.py --labelSurf lh.aparcNMMjt.vtk --segFundi lh.pial.fundi.vtk --segSulci lh.pial.sulci.vtk\
#        --depth lh.travel.vtk --curvmean lh.curvature.mean.vtk --curvmin lh.curvature.min.vtk\
#        --curvmax lh.curvature.max.vtk --curvgauss lh.curvature.gauss.vtk --thick lh.thickness --convex lh.sulc \n"
        print "    python shape_table.py --labelSurf ../KKI2009-11/lh.aparcNMMjt.pial.vtk\
        --segFundi ../KKI2009-11/lh.travel.depth.depth.fundi.vtk --segSulci ../KKI2009-11/lh.pial.sulci.seg.vtk\
        --depth ../KKI2009-11/lh.travel.vtk -- curvmean ../KKI2009-11/lh.pial.curvature.mean.vtk\
        --curvmin ../KKI2009-11/lh.pial.curvature.min.vtk --curvmax ../KKI2009-11/lh.pial.curvature.max.vtk\
        --curvgauss ../KKI2009-11/lh.pial.curvature.gauss.vtk --thick ../KKI2009-11/lh.thickness\
        --convex ../KKI2009-11/lh.sulc --segFundi ../KKI2009-11/lh.travel.depth.depth.fundi.vtk\
        --segSulci ../KKI2009-11/lh.pial.sulci.seg.vtk --depth ../KKI2009-11/lh.travel.vtk\
        --curvmean ../KKI2009-11/lh.pial.curvature.mean.vtk --curvmin ../KKI2009-11/lh.pial.curvature.min.vtk\
        --curvmax ../KKI2009-11/lh.pial.curvature.max.vtk --curvgauss ../KKI2009-11/lh.pial.curvature.gauss.vtk\
        --thick ../KKI2009-11/lh.thickness --convex ../KKI2009-11/lh.sulc\
        --patchtable ../KKI2009-11/lh.patch.tsv --funditable ../KKI2009-11/lh.fundi.tsv --sulcitable ../KKI2009-11/lh.sulci.tsv"
        
    def check_opt(opts, args): # a nested function
        '''Check whether opts and args satisfy constraints 
        '''            
        # Check whether all options have their values given. Second surface unneeded for shape table generation.
        for o,p in opts:
            if p == "":
                print "[ERROR] The option",o, "is missing parameter. Please check usage and provide it."
    #            print_help()
                sys.exit()
    
    def process_opt(opts,args): # a nested function
        '''Give option parameters to variables that represent I/O file names.
        Do this AFTER check_opt.
        '''
        
        InputColor = "\033[32m"
        OutputColor = "\033[36m"
        EndColor = "\033[0m"
        
        ThickFile, ConvexFile = '', ''
#        FundiVTK, PitsVTK, SulciVTK = '','',''
        DepthVTK, CurvMinVTK, CurvMaxVTK, CurvMeanVTK, CurvGaussVTK = "","","","",""
        LabelSurfVTK, SegSulciVTK, SegFundiVTK = "","", ""
        
        Patch_Table_File, Sulci_Table_File, Fundi_Table_File = "", "", ""
                
        for o,p in opts:
            if o=="--labelSurf":
                LabelSurfVTK  = p
                print "  [Input] labeled surface:" + InputColor + LabelSurfVTK + EndColor
            elif o=="--segSulci":
                SegSulciVTK  = p
                print "  [Input] segmented sulci:" + InputColor + SegSulciVTK + EndColor
            elif o=="--segFundi":
                SegFundiVTK  = p
                print "  [Input] segmented fundi:" + InputColor + SegFundiVTK + EndColor
            elif o=="--convex":
                ConvexFile = p
                print "  [Input] convexity file:" + InputColor + ConvexFile + EndColor
            elif o=="--thick":
                ThickFile = p
                print "  [Input] thickness file:" + InputColor + ThickFile + EndColor
            elif o=="--depth":
                DepthVTK = p
                print "  [Input] depth file:" + InputColor + DepthVTK + EndColor
            elif o=="--curvmin":
                CurvMinVTK = p
                print "  [Input] minimal curvature file:" + InputColor + CurvMinVTK + EndColor
            elif o=="--curvmax":
                CurvMaxVTK = p
                print "  [Input] maximal curvature file:" + InputColor + CurvMaxVTK + EndColor
            elif o== "--curvmean":
                CurvMeanVTK = p
                print "  [Input] mean curvature file:" + InputColor + CurvMeanVTK + EndColor
            elif o== "--curvgauss":
                CurvGaussVTK = p
                print "  [Input] Gaussian curvature file:" + InputColor + CurvGaussVTK + EndColor
            elif o== "--patchtable":
                Patch_Table_File = p
                print "  [Output] Patch table file:" + OutputColor + Patch_Table_File + EndColor
            elif o =="--sulcitable":
                Sulci_Table_File = p
                print "  [Output] sulci table file:" + OutputColor + Sulci_Table_File + EndColor
            elif o =="--funditable":
                Fundi_Table_File = p
                print "  [Output] fundi table file:" + OutputColor + Fundi_Table_File + EndColor
            else:
                print "Unrecongnized option", o 
                print_help()
        return  [ThickFile, ConvexFile,\
        DepthVTK, CurvMinVTK, CurvMaxVTK, CurvMeanVTK, CurvGaussVTK,\
        LabelSurfVTK, SegSulciVTK, SegFundiVTK], [Patch_Table_File, Sulci_Table_File, Fundi_Table_File] 
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"",\
        ["thick=","convex=","segFundi=","segSulci=","labelSurf=","depth=","curvmin=","curvmax=","curvmean=","curvgauss=","patchtable=","sulcitable=","funditable="])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        print_help()
        sys.exit(2)

    def load_LUTs(InputFileList):
        """Load all look-up tables, 4 curvatures, depth, thickness, convexity, labeled surfaces and fundi and sulci segmentation results
        """
        import io_vtk
        import io_file 

#        print InputFileList

        [ThickFile, ConvexFile, DepthVTK, CurvMinVTK, CurvMaxVTK, CurvMeanVTK, CurvGaussVTK, LabelSurfVTk, SegSulciVTK, SegFundiVTK] =\
        InputFileList 

        Thickness, Convexity, Depth, CurvMin, CurvMax, CurvMean, CurvGauss, LabelGyral, SegSulci, SegFundi=\
        [], [], [], [], [], [], [], [], [], []
                       
        # load freesurfer type maps
        if ThickFile != "":
            Thickness = io_file.readCurv(ThickFile)
        if ConvexFile != "":
            Convexity = io_file.readCurv(ConvexFile)

        LUTs = [Thickness, Convexity, Depth, CurvMin, CurvMax, CurvMean, CurvGauss, LabelGyral, SegSulci, SegFundi]
        
        for i in xrange(2,10-2):
            Points, Faces, LUTs[i] = io_vtk.load_VTK_Map(InputFileList[i])
        
        return LUTs
        
    
    def load_fundi(segFundiVTK):
        """Load fundus lines from segFundiVTK
        """
        import io_vtk
        Fundi = io_vtk.load_VTK_line(segFundiVTK)
        print "Fundi loaded from", segFundiVTK
        return Fundi
        
    def load_sulci(segSulciVTK):
        """Load sulcus vertexes from segSulciVTK 
        """
        import io_vtk
        Sulci = io_vtk.load_VTK_vertex(segSulciVTK)
        print "sulci loaded from", segSulciVTK
        return Sulci

    # main part of input_handle()
    check_opt(opts, args)
    InputFiles, OutputFiles = process_opt(opts, args)
    LUTs = load_LUTs(InputFiles)
#    print len(LUTs)
#    print InputFiles[8:9]
    
    segSulciVTK  = InputFiles[8]
    Sulci, LUTs[8] = load_sulci(segSulciVTK)
        
    segFundiVTK = InputFiles[9]
    Fundi, LUTs[9] = load_fundi(segFundiVTK)

#    
    return LUTs, Fundi, Sulci, OutputFiles
    

def generate_shape_table(LUTs, Fundi, Sulci, GyralLabels, SulcalLabels, OutputFiles):
    """Generate shape tables for labeled surface patches, segmented fundi and segmented sulci
    Thus, 3 tables.  
    
    Inputs
    =======
    
    GyralLabels : list of integers
        A list of all possible regional/gyral labels 
        
    SulcalLabels : list of integers
        A list of all possible sulcal labels
    
    """
    def label_patch_shape_table(LUTs, GyralLabels):
        """Generate the shape table for labeled surface patches
        
        Inputs
        =========
        
        Index : list of integers 
            IDs of vertexes whose labels are one Gyral label. 
        
        """
        from numpy import mean

        Measures = LUTs[:7]    
        Label = map(int, LUTs[7])
        
        Table = []
        
        for GyralLabel in GyralLabels:
            Index = [i for i in xrange(len(Label)) if Label[i] == GyralLabel]  # IDs of vertexes whose labels are the gyralLabel
            Measure_of_Row = [ [Measure[i] for i in Index ] for Measure in Measures]
            Row = map(mean, Measure_of_Row)
            Row = [GyralLabel] + Row
            Table.append(Row)
            
        return Table

    def seg_sulci_shape_table(LUTs, SulcalLabels, Sulci):
        """Generate the shape table for segmented sulci
        
        Inputs
        =========
        
        Index : list of integers 
            IDs of vertexes whose labels are one Gyral label. 
        
        """
        from numpy import mean
    
        Measures = LUTs[:7]
        SegSulci = map(int, LUTs[8])
        
        Table = []
        
        for SulcalLabel in SulcalLabels:
#            Index = [i for i in xrange(len(SegSulci)) if ((SegSulci[i] == SulcalLabel) and (i in Sulci))] 
            Index = [i for i in Sulci if SegSulci[i] == SulcalLabel]
            Measure_of_Row = [ [Measure[i] for i in Index ] for Measure in Measures]
            Row = map(mean, Measure_of_Row)
            Row = [SulcalLabel] + Row
            Table.append(Row)
#            break # for quick testing 
            
        return Table

    def seg_fundi_shape_table(LUTs, GyralLabels, Fundi):
        """Generate the shape table for segmented fundi
        
        Inputs
        =========
        
        Index : list of integers 
            IDs of vertexes whose labels are one Gyral label. 
        
        """
        from numpy import mean
    
        Measures = LUTs[:7]
        SegFundi = map(int, LUTs[9])
        Fundus_Vertexes = []
        for Lines in Fundi:
            Fundus_Vertexes += Lines
        Fundus_Vertexes = list(set(Fundus_Vertexes))
#        print Fundus_Vertexes
#        print SegFundi

#        print [i for i in xrange(len(SegFundi)) if ((SegFundi[i] != -1) and (i in Fundus_Vertexes))]
        
        Table = []
        
        for SulcalLabel in SulcalLabels:
#            Index = [i for i in xrange(len(SegFundi)) if ((SegFundi[i] == SulcalLabel) and (i in Fundus_Vertexes))]
            Index = [i for i in Fundus_Vertexes if SegFundi[i] == SulcalLabel]
#            print SulcalLabel, Index 
            Measure_of_Row = [ [Measure[i] for i in Index ] for Measure in Measures]
            Row = map(mean, Measure_of_Row)
            Row = [SulcalLabel] + Row
            Table.append(Row)
            
#            break # for quick testing 
            
        return Table

    def write_shape_table(Table, Filename):
        Fp = open(Filename, 'w')
        Header = ['label','thickness','convexity','depth','curvmin','curvmax','curvmean','curvgauss']
        Line = "\t".join(Header)
        Fp.write(Line+"\n")
        for Row in Table:
#            print Row
            Row = map(str,Row)
            Line = "\t".join(Row)
            Fp.write(Line+"\n")
        Fp.close()

    # unpack Outputfiles
    [Patch_Table_File, Sulci_Table_File, Fundi_Table_File] = OutputFiles
    
    # shape table for each labeled surface patches
    Label_Patch_Table = label_patch_shape_table(LUTs, GyralLabels)
    write_shape_table(Label_Patch_Table, Patch_Table_File)
    
    # shape table for each sulci patch of common labels
    Sulci_Table = seg_sulci_shape_table(LUTs, SulcalLabels, Sulci)
    write_shape_table(Sulci_Table, Sulci_Table_File)

    # shape table for each fundi segment of commom labels
    Fundi_Table = seg_fundi_shape_table(LUTs, SulcalLabels, Fundi)
    write_shape_table(Fundi_Table, Fundi_Table_File)
    
    return Label_Patch_Table, Sulci_Table, Fundi_Table
    


if __name__ == "__main__":
    
    GyralLabels = range(50)
    SulcalLabels = range(50)
    
    LUTs, Fundi, Sulci, OutputFiles = input_handle()
    Label_Patch_Table, Sulci_Table, Fundi_Table = generate_shape_table(LUTs, Fundi, Sulci, GyralLabels, SulcalLabels, OutputFiles)

#    # Now begins the code to generate shape table 
#    import cPickle
#    
#    #S,L,FV,FL = load_multi_segs('/forrest/data/MRI/MDD/', '.euclidean.depth.depth.fundi.seg.vtk', Quit= 500)
#    #S,L,FV,FL = load_multi_segs('/forrest/data/MRI/MDD/', '.sulc.fundi.from.pits.pial.seg.vtk', Quit= 500)
#    #with open('labels.pickle', 'w') as f:
#    #    cPickle.dump([S,L,FV,FL], f)
#    
#    with open('labels.pickle') as f:
#        S,L,FV,FL = cPickle.load(f)
#    
#    Common = common_pairs(S)
#    
#    gen_shape_all('/forrest/data/MRI/MDD/',Common, S, L, FV, FL)
#    
#    #Test_Groups, TestLength = largest_groups(Common, FL[0], FV[0], S[0])
#    #gen_shape(Test_Groups, L[0], 'test_shape.tsv', TestLength) 

