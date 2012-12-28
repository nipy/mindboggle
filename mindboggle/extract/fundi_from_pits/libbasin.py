#!/usr/bin/python
"""
Extracting features from VTK input files. 

Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

For algorithmic details, please check:     
    Forrest S. Bao, et al., Automated extraction of nested sulcus features from human brain MRI data, 
    IEEE EMBC 2012, San Diego, CA

Dependencies:
    python-vtk: vtk's official Python binding
    numpy
    io_vtk : under mindboggle/utils

   
"""
 
from numpy import mean, std, median, array, zeros, eye, flatnonzero, sign, matrix, zeros_like
import os.path
import cPickle
#import io_vtk  # Assummng io_vtk is in PYTHONPATH
from mindboggle.utils import io_vtk
import sys


#-----------------Begin function definitions------------------------------------------------------------- 

def fcNbrLst(FaceDB, Hemi):
    '''Get a neighbor list of faces, also the vertex not shared with current face
    
    Data structure:
    
    NbrLst: a list of size len(FaceDB)
    NbrLst[i]: two lists of size 3 each. 
    NbrLst[i][0] = [F0, F1, F2]: F0 is the neighbor of face i facing V0 where [V0, V1, V2] is face i. And so forth.  
    NbrLst[i][1] = [V0p, V1p, V2p]: V0p is the vertex of F0 that is not shared with face i        
    
    '''
    
    
    NbrFile = Hemi + '.fc.nbr'

    if os.path.exists(NbrFile):
        #return fileio.loadFcNbrLst(NbrFile)
        print "loading  face nbr lst from:" , NbrFile
        Fp = open(NbrFile, 'r')
        NbrLst = cPickle.load(Fp)
        Fp.close()
        return NbrLst        

    print "calculating face neighbor list"
    
    FaceNo = len(FaceDB)
    
    NbrLst = []
    [NbrLst.append([[-1,-1,-1], [-1,-1,-1]]) for i in xrange(FaceNo)]  
    
    Done =[]
    [Done.append(0) for i in xrange(FaceNo)]
    
    for i in xrange(0, FaceNo):
#    for i in xrange(0, 2600+1):
#        print i
        Face = FaceDB[i]
#        [V0, V1, V2] = Face
#        Found = 0  # if Found  == 1, no need to try other faces
        for j in xrange(i+1, FaceNo):
            AnotherFace = FaceDB[j]

            for Idx in xrange(0,2):
                ChkFc1 = Face[Idx]
                for ChkFc2 in Face[Idx+1:3]:
                    if ChkFc1 in AnotherFace:
                        if ChkFc2 in AnotherFace:
                            NbrID1 = 3 - Face.index(ChkFc1) - Face.index(ChkFc2)   # determine it's F0, F1 or F2.  
                            NbrLst[i][0][NbrID1] = j
                    
                            NbrID2 = 3 - AnotherFace.index(ChkFc1) - AnotherFace.index(ChkFc2)   # determine it's F0, F1 or F2.
                            NbrLst[j][0][NbrID2] = i
                    
#                        Vp1 = AnotherFace[NbrID2]# determine V{0,1,2}p
#                        Vp2 = Face[NbrID1]# determine V{0,1,2}p
                            NbrLst[i][1][NbrID1] = AnotherFace[NbrID2]
                            NbrLst[j][1][NbrID2] = Face[NbrID1]
                    
                            Done[i] += 1
                            Done[j] += 1
                    
            if Done[i] ==3:
                break  # all three neighbors of Face has been found 

    Fp = open(NbrFile, 'w')
    
# Commented 2011-11-27 23:54     
#    for i in xrange(0, len(FaceDB)
#        for j in NbrLst[i]:
#            Fp.write(str(j[0]) + '\t' + str(j[1]) + '\t' + str(j[2]) + '\t')
#        Fp.write('\n')
# End of Commented 2011-11-27 23:54

    cPickle.dump(NbrLst, Fp)

    Fp.close()
                    
    return NbrLst

def vrtxNbrLst(VrtxNo, FaceDB, Hemi):
    """Given the number of vertexes and the list of faces, find the neighbors of each vertex, in list formate. 
    """

    NbrFile = Hemi + '.vrtx.nbr'
        
    if os.path.exists(NbrFile):
        #return fileio.loadVrtxNbrLst(NbrFile) # change to cPickle
        print "Loading vertex nbr lst from:", NbrFile
        Fp = open(NbrFile, 'r')  # need to use cPickle
        NbrLst = cPickle.load(Fp)
        Fp.close()
        return NbrLst
       
    print "Calculating vertex neighbor list"
    
    NbrLst = [[] for i in xrange(0, VrtxNo)]
        
    for Face in FaceDB:
        [V0, V1, V2] = Face
        
        if not V1 in NbrLst[V0]:
            NbrLst[V0].append(V1)
        if not V2 in NbrLst[V0]:
            NbrLst[V0].append(V2)

        if not V0 in NbrLst[V1]:
            NbrLst[V1].append(V0)
        if not V2 in NbrLst[V1]:
            NbrLst[V1].append(V2) 

        if not V0 in NbrLst[V2]:
            NbrLst[V2].append(V1)
        if not V1 in NbrLst[V2]:
            NbrLst[V2].append(V1)
    
    Fp = open(NbrFile, 'w')  # need to use cPickle

# Commented 2011-11-27 23:54    
#    for i in xrange(0, VrtxNo):
#        [Fp.write(str(Vrtx) + '\t') for Vrtx in NbrLst[i]]
#        Fp.write('\n')    
# End of Commented 2011-11-27 23:54
    cPickle.dump(NbrLst, Fp)

    Fp.close()

    return NbrLst

def compnent(FaceDB, Basin, NbrLst, PathHeader):
    '''Get connected component, in each of all basins, represented as faces and vertex clouds
    
    Parameters
    -----------
    
    NbrLst : list
        neighbor list of faces, NOT VERTEXES
        
    PathHeader : header of the path to save component list
    
    '''
    
    FcCmpntFile = PathHeader + '.cmpnt.face' 
    VrtxCmpntFile = PathHeader + '.cmpnt.vrtx'
        
    if os.path.exists(FcCmpntFile) and os.path.exists(VrtxCmpntFile):
#        return fileio.loadCmpnt(FcCmpntFile), fileio.loadCmpnt(VrtxCmpntFile)
        print "Loading Face Components from:", FcCmpntFile
        Fp = open(FcCmpntFile, 'r')
        FcCmpnt = cPickle.load(Fp)
        Fp.close()
        
        print "Loading Vertex Components from:", VrtxCmpntFile
        Fp = open(VrtxCmpntFile, 'r')
        VrtxCmpnt = cPickle.load(Fp)
        Fp.close()
        return FcCmpnt, VrtxCmpnt
        
    print "calculating face and vertex components"
    
    Visited = [False for i in xrange(0, len(Basin))]
    
    FcCmpnt, VrtxCmpnt = [], [] 
    
    while not allTrue(Visited):
        Seed = dfsSeed(Visited, Basin)# first basin face that is not True in Visited
#        print Seed
        Visited, FcMbr, VrtxMbr = dfs(Seed, Basin, Visited, NbrLst, FaceDB)# DFS to fine all connected members from the Seed
        FcCmpnt.append(FcMbr)
        VrtxCmpnt.append(VrtxMbr)
            
#    fileio.writeCmpnt(FcCmpnt, FcCmpntFile)
#    fileio.writeCmpnt(VrtxCmpnt, VrtxCmpntFile)
    Fp = open(FcCmpntFile, 'w')
    cPickle.dump(FcCmpnt, Fp)
    Fp.close()
    Fp = open(VrtxCmpntFile, 'w')
    cPickle.dump(VrtxCmpnt, Fp)
    Fp.close()
                
    return FcCmpnt, VrtxCmpnt

def judgeFace1(FaceID, FaceDB, CurvatureDB, Threshold = 0):
    """Check whether a face satisfies the zero-order criterion
    
    If all three vertexes of a face have negative curvature, return True. O/w, False.
    
    Input
    ======
    
        FaceID: integer
                the ID of a face, indexing from 0
        
        FaceDB: list
                len(FaceDB) == number of faces in the hemisphere
                FaceDB[i]: a 1-D list of the IDs of three vertexes that consist of the i-th face
                
        CurvatureDB: list
                len(CurvatureDB) == number of vertexes in the hemisphere
                CurvatureDB[i]: integer, the curvature of the i-th vertex 
    
    """
    
    [V0, V1, V2] = FaceDB[FaceID]
##    
#    if (CurvatureDB[V0] > Threshold) and (CurvatureDB[V1] > Threshold) and (CurvatureDB[V2] > Threshold):
#        return True
#    else:
#        return False
##
    if (CurvatureDB[V0] <= Threshold) or (CurvatureDB[V1] <= Threshold) or (CurvatureDB[V2] <= Threshold):
        return False
    else:
        return True    

def basin(FaceDB, CurvatureDB, Threshold = 0):
    '''Given a list of faces and per-vertex curvature value, return a list of faces comprising basins
    '''
    Basin = []
    Left = []
    for FaceID in xrange(0, len(FaceDB)):
        if judgeFace1(FaceID, FaceDB, CurvatureDB, Threshold = Threshold):
            Basin.append(FaceID)
        else:
            Left.append(FaceID)
    
    return Basin, Left

def allTrue(List):
    '''Check whether a logical list contains non-True elements. 
    '''    
#    for Bool in List:
#        if not Bool:
#            return False
#    return True

    return all(x==True for x in List) 

def dfsSeed(Visited, Basin):
    '''Given a list of faces comprising the basins, find a face that has not been visited which will be used as the seeding point for DFS.
    '''
    for i in xrange(0, len(Visited)):
        if not Visited[i]:
            return Basin[i]
             
def dfs(Seed, Basin, Visited, NbrLst, FaceDB):
    '''Return all members (faces and vertexes) of the connected component that can be found by DFS from a given seed point
    
    Parameters
    -----------
    
    NbrLst : list
        neighbor list of faces, NOT VERTEXES    
     
    '''
    Queue = [Seed]
    FcMbr = []  # members that are faces of this connected component
    VrtxMbr = [] # members that are vertex of this connected component
    while Queue != []:
#        print Queue
        Seed = Queue.pop()
        if Seed in Basin:
            if not Visited[Basin.index(Seed)]:
                Visited[Basin.index(Seed)] = True
                FcMbr.append(Seed)
                for Vrtx in FaceDB[Seed]:
                    if not (Vrtx in VrtxMbr):
                        VrtxMbr.append(Vrtx)
                Queue += NbrLst[Seed][0]
            
    return Visited, FcMbr, VrtxMbr

def pmtx(Adj):
    '''Print a matrix as shown in MATLAB stdio
    '''
    
    for j in xrange(0,25):
        print j,
    print '\n'
     
    for i in xrange(0, 25):
        print i,
        for j in xrange(0, 25):
            print Adj[i,j],
        print '\n'    

def all_same(items):
    return all(x == items[0] for x in items)

def univariate_pits(CurvDB, VrtxNbrLst, VrtxCmpnt, Thld): 
    '''Finding pits using one variable, e.g., depth. 
     
    '''
    
    print "Extracting pits"
    
#    Stack, P, Child, M, B, End, L  = [], [], {}, -1, [], {}, 10
    C = [-1 for i in xrange(0, len(VrtxNbrLst))]
    Child = {}
    End = {}
    M = -1
    B = []

    for Cmpnt in VrtxCmpnt:  # for each component
        
        Curv=dict([(i, CurvDB[i]) for i in Cmpnt])
        
        Stack = []
        for Vrtx, Cvtr in sorted(Curv.iteritems(), key=lambda (k,v): (v,k)):
            Stack.append(Vrtx)
    
        Visited = []
        while len(Stack) >0:
            Skip_This_Vrtx = False # updated Forrest 2012-02-12, skip vertexes whose neighbors are not in the component to denoise
            Vrtx = Stack.pop()
            
            WetNbr = []
            NbrCmpnt = []
            for Nbr in list(set(VrtxNbrLst[Vrtx])):
                if not Nbr in Cmpnt:
                    Skip_This_Vrtx = True 
                if Nbr in Visited:  # This condition maybe replaced by If C[Vrtx] ==-1
                    WetNbr.append(Nbr)
                    if C[Nbr] != -1: 
                        NbrCmpnt.append(C[Nbr])
            
            if Skip_This_Vrtx :
                continue
            
            Visited.append(Vrtx)
            
            if len(WetNbr) == 1: # if the vertex has one neighbor that is already wet     
                [Nbr] = WetNbr
                if End[C[Nbr]]:
                    C[Vrtx] = Child[C[Nbr]]
                else:
                    C[Vrtx] = C[Nbr]
    #                print C[Nbr], "==>", C[V] 
            elif len(WetNbr) >1 and all_same(NbrCmpnt): # if the vertex has more than one neighbors which are in the same component
                if End[NbrCmpnt[0]]:
                    C[Vrtx] = Child[NbrCmpnt[0]]
                else:
                    C[Vrtx] = NbrCmpnt[0]
            elif len(WetNbr) >1 and not all_same(NbrCmpnt):  # if the vertex has more than one neighbors which are NOT in the same component
                M += 1
                C[Vrtx] = M
                for Nbr in WetNbr:
                    Child[C[Nbr]] = M
                    End[C[Nbr]] = True
                End[M] = False
#            elif : # the vertex's neighbor are not fully in the component 
            else:
                M += 1
                if CurvDB[Vrtx] > Thld:
                    B.append(Vrtx)
                End[M] = False
                C[Vrtx] = M
    return B, C, Child

def clouchoux(MCurv, GCurv):
    '''Judge whether a vertex is a pit in Clouchoux's definition
    
    Parameters
    ===========
    
    MCurv: float
        mean curvature of a vertex
        H in Clouchoux's paper 
    
    GCurv: float
        mean curvature of a vertex 
        K in Clouchoux's paper
      
    Returns 
    ========
    
    True if this is a pit. False, otherwise. 
    
    Notes
    =========
    
    (Since Joachim's code updates all the time, this settings has to be updated accordingly)
    
    In Clochoux's paper, the following definitions are used:
        H > 0, K > 0: pit, in Clouchoux's paper
        H < 0, K > 0: peak, in Clouchoux's paper    
    
    If features are computed by ComputePricipalCurvature(), 
    use this settings to get proper pits:
        H > 3, K < 0 (curvatures not normalized)
        H > 0.2, K < 0 (curvatures normalized)
        
    
    '''

#    if (MCurv > 3) and (GCurv < 0):
    if  (MCurv > 0.2) and (GCurv < 0):
        return True
    else:
        return False
    
def clouchoux_pits(Vertexes, MCurv, GCurv):
    '''Extract pits using Clouchoux's definition
    '''
    
    Pits = []

    for i in xrange(len(Vertexes)):
        if clouchoux(MCurv[i], GCurv[i]):
            Pits.append(i)
    
    print  len(Pits), "Pits found"
    
    return Pits

def getBasin_and_Pits(Maps, Mesh, SulciVTK, PitsVTK, SulciThld = 0, PitsThld = 0, Quick=False, Clouchoux=False, SulciMap='depth'):
    '''Extracting basin and pits (either local minimum approach or Clouchoux's)
    
    Parameters
    =============
        Maps: dictionary
            Keys are map names, e.g., depth or curvatures.
            Values are per-vertex maps, e.g., curvature map.
                        
        Mesh: 2-tuple of lists
            the first list has coordinates of vertexes while the second defines triangles on the mesh
            This is a mandatory surface, normally a non-inflated surface. 

        SulciThld: float
            the value to threshold the surface to separate sulci and gyri
            
        PitsThld: float
            vertexes deeper than this value can be considered as pits
            
        Quick: Boolean
            If true, extract sulci only (no component ID, only thresholding), skipping pits and later fundi. 
            
        Clouchoux: Boolean
            If true, extract pits using Clouchoux's definition. O/w, local minimum approach.
            
        SulciMap: string
            The map to be used to get sulci
            by default, 'depth'

    '''
    
    def write_surface_with_LUTs(File, Points, Faces, Maps):
        """Like write_scalars in io_vtk but no writing of vertices 
        """
        print "writing sulci into VTK file:", File

        Fp = open(File,'w')
        io_vtk.write_vtk_header(Fp)
        io_vtk.write_vtk_points(Fp, Points)
        io_vtk.write_vtk_faces(Fp, Faces)
        if len(Maps) > 0:
        # Make sure that LUTs is a list of lists
            Count = 0
            for LUT_name, LUT in Maps.iteritems():
                if Count == 0 :
                    io_vtk.write_vtk_scalars(Fp, LUT, LUT_name)
                else:
                    io_vtk.write_vtk_scalars(Fp, LUT, LUT_name, begin_scalars=False)
                Count += 1
        Fp.close()
        return None
        
    def write_pits_without_LUTs(File, Points, Indexes):
        """Like write_scalars in io_vtk but no writing of vertices 
        """
        print "writing pits into VTK file:", File
        Fp = open(File,'w')
        io_vtk.write_vtk_header(Fp)
        io_vtk.write_vtk_points(Fp, Points)
        io_vtk.write_vtk_vertices(Fp, Indexes)
        Fp.close()
        return None
    
    print "\t thresholding the surface using threshold = ", SulciThld
    
    [Vertexes, Faces] = Mesh
        
    MapBasin = Maps[SulciMap]   
    Basin, Gyri = basin(Faces, Maps[SulciMap], Threshold = SulciThld)

    if not Quick:
        LastSlash = len(SulciVTK) - SulciVTK[::-1].find('/')
        Hemi =  SulciVTK[:SulciVTK[LastSlash:].find('.')+LastSlash]# path up to which hemisphere, e.g., /home/data/lh
        VrtxNbr = vrtxNbrLst(len(Vertexes), Faces, Hemi)
        FcNbr   = fcNbrLst(Faces, Hemi)
        FcCmpnt, VrtxCmpnt = compnent(Faces, Basin, FcNbr, ".".join([Hemi, SulciMap, str(SulciThld)])) 
        CmpntLUT = [-1 for i in xrange(len(MapBasin))]
        for CmpntID, Cmpnt in enumerate(VrtxCmpnt):
            for Vrtx in Cmpnt:
                CmpntLUT[Vrtx] = CmpntID
                
        Maps['CmpntID'] = CmpntLUT

        if Clouchoux:
            Pits = clouchoux_pits(Vertexes, Maps['meancurv'], Maps['gausscurv'])
        else: # local minimum approach  
            MapPits =  Maps[SulciMap]  # Users will get the option to select pits extraction map in the future.
            Pits, Parts, Child = univariate_pits(MapPits, VrtxNbr, VrtxCmpnt, PitsThld)
            Maps['hierarchy'] = Parts
        
    else:
        print "\t\t Thresholding the surface to get sulci only."
    
    Faces = [map(int,i) for i in Faces]# this is a temporal fix. It won't cause precision problem because sys.maxint is 10^18.
    Vertexes = map(list, Vertexes)    
    write_surface_with_LUTs(SulciVTK, Vertexes, [Faces[i] for i in Basin], Maps)    

    if Quick:
        sys.exit()
        
    write_pits_without_LUTs(PitsVTK, Vertexes, Pits)

    # output tree hierarchies of basal components
#    print "writing hierarchies of basal components"
#    WetFile = PrefixExtract + '.pits.hier'
#    WetP = open(WetFile,'w')
#    for LowComp, HighComp in Child.iteritems():
#        WetP.write(str(LowComp) + '\t' + str(HighComp) + '\n')
#    WetP.close()
    # end of output tree hierarchies of basal components
    
# End of Get pits Forrest 2011-05-30 10:16

# a monolithic code output each component
#    Dic = {}
#    for CID, Cmpnt in enumerate(FcCmpnt):
#        Dic[CID] = len(Cmpnt)
#    
#    #Dic = sorted(Dic.iteritems(), key= lambda (k,v,) : (v,k))
#    Counter = 1
#    for CID, Size in sorted(Dic.iteritems(), key=lambda (k,v): (v,k)):  
##        print Size
#        Rank = len(FcCmpnt) - Counter +1
#        Fp = open(BasinFile + '.' + SurfFile[-1*SurfFile[::-1].find('.'):] + '.' + str(Rank) +'-th.vtk','w')
#        Vertex, Face = fileio.readSurf(SurfFile)
#        FundiList = FcCmpnt[CID]
#        libvtk.wrtFcFtr(Fp, Vertex, Face, FundiList)
#        Fp.close()
#        Counter += 1
    
# a monolithic code output each component
           
#---------------End of function definitions---------------------------------------------------------------  