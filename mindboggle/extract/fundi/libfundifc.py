# All functions to extract fundi as face strips 
# Forrest Sheng Bao (GPL) 2011
# This file is called libfundifc.py before 2011-10-08

import fileio, libvtk, libbasin
from numpy import std, mean, median
import os

#def judgeFace1(FaceID, FaceDB, CurvatureDB, CenterDB = [], Threshold = 0):
def judgeFace1(FaceID, CenterDB, Threshold = 0):
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
    
# if there is CenterDb
#    if CenterDB != []:
    if CenterDB[FaceID] > Threshold:
        return True
    else:
        return False   
    
# if there is NO CenterDb --    
    
#    [V0, V1, V2] = FaceDB[FaceID]
##    if ((CurvatureDB[V0] > 0) and (CurvatureDB[V1] > 0)) \
##    or ((CurvatureDB[V1] > 0) and (CurvatureDB[V2] > 0)) \
##    or ((CurvatureDB[V0] > 0) and (CurvatureDB[V2] > 0)):    
#    if (CurvatureDB[V0] > Threshold) and (CurvatureDB[V1] > Threshold) and (CurvatureDB[V2] > Threshold):
#        return True
#    else:
#        return False

def getCenter(FaceDB, CurvDB):
    '''Compute the average per-vertex value of a face 
    '''
    
    Center = []
    
    for i in xrange(0,len(FaceDB)):
        [V0, V1, V2] = FaceDB[i]
        Center.append(mean([CurvDB[V0], CurvDB[V1], CurvDB[V2]]))

    return Center
   
#def thirdVertex(FaceID, V0, V1, FaceDB):
#    """Given the two vertexes' ID and a face ID, find the ID of the third vertex in the face
#    
#    Test
#    =====
#        import curvature as curv
#        >>> curv.thirdVertex(0, 1, 3, [[1,2,3]])
#        2
#    
#    """
#    Face = FaceDB[FaceID]
#    V0, V1 = Face.index(V0) , Face.index(V1)  # find the index of V0 and V1 in FaceDB[FaceID], e.g, 1 and 3
#    V2 = 3 - V0 - V1  # index of V2 in FaceDB[FaceID]
#    return Face[V2]

def neighborFaceVertex(FaceID, V0, V1, NbrLst, FaceDB):
    """ Given two vertexes and a face, find the face and the third vertex of the face that co-edges the two vertexes
        
    Test
    ======
        import curvature as curv
        >>> curv.neighborFaceVertex(0, 1, 3, [[1,2,3],[1,7,4],[3,1,5]])
        (2, 5)
    
    """ 
    
    Nbrs = NbrLst[FaceID]
    Face = FaceDB[FaceID]
    Idx = 3 - Face.index(V0) - Face.index(V1) 
    
    return Nbrs[0][Idx], Nbrs[1][Idx]
    
    
def judgeFace2(FaceID, FaceDB, NbrLst, CurvatureDB, Lb=-0.04, Ub=0.12, Level = 1):
    """Check whether a face satisfies the first-order criterion
    
    Input
    ======
        FaceID: integer
            the ID of a face, indexing from 0
        
        FaceDB: list
            len(FaceDB) == number of faces in the hemisphere
            FaceDB[i]: a 1-D list of the IDs of three vertexes that consist of the i-th face
        
        NbrLst: list
            len(FaceDB) == number of faces in the hemisphere
            The neighbor list
                
        CurvatureDB: list
            len(CurvatureDB) == number of vertexes in the hemisphere
            CurvatureDB[i]: integer, the curvature of the i-th vertex 
                
        Level: integer
            There are many ways to judge whether a face can satisfy the 1st-order derivative test, e.g., comparing each neighbor
            with the face or comparing the average value of all neighbors with the face. Level allows us to choose among different
            ways. 
            
            Level = 1: A face has to satisfy a condition with each of its neighbors/directions
            Level = 2: A face has to satisfy a condition with the average of all its neighbors/directions 
    
    """
    [V0, V1, V2] = FaceDB[FaceID]
    
    #===========================================================================
    # The vertex against V0 and V1 is V3. The one against V0 and V2 is V4. 
    # So the one against V1 and V2 is V5. 
    #===========================================================================
    
    Tmp, V3 = neighborFaceVertex(FaceID, V0, V1, NbrLst, FaceDB)
    Tmp, V4 = neighborFaceVertex(FaceID, V0, V2, NbrLst, FaceDB)
    Tmp, V5 = neighborFaceVertex(FaceID, V1, V2, NbrLst, FaceDB)
    
    if Level == 4:
        if ((Ub > CurvatureDB[V3] - CurvatureDB[V2] >= Lb) and (Ub > CurvatureDB[V5] - CurvatureDB[V0] >= Lb)) \
        or ((Ub > CurvatureDB[V5] - CurvatureDB[V0] >= Lb) and (Ub > CurvatureDB[V4] - CurvatureDB[V1] >= Lb))\
        or ((Ub > CurvatureDB[V4] - CurvatureDB[V1] >= Lb) and (Ub > CurvatureDB[V3] - CurvatureDB[V2] >= Lb)):
            return True
        else:
            return False
#        Diff=[CurvatureDB[V3] - CurvatureDB[V2], CurvatureDB[V5] - CurvatureDB[V0], CurvatureDB[V4] - CurvatureDB[V1]]
#        Judge = [Ub > D > Lb for D in Diff ] 
#        FalseCount = 0
#        for J in Judge:
#            if not J:
#                FalseCount += 1
#            
#        if FalseCount <= 0.1 * 3:
#            return True
#        else:
#            return False
    elif Level == 1:
        if (Ub > CurvatureDB[V3] - CurvatureDB[V2] >= Lb) and (Ub > CurvatureDB[V5] - CurvatureDB[V0] >= Lb) and (Ub > CurvatureDB[V4] - CurvatureDB[V1] >= Lb):   
            return True
        else:
            return False
    elif Level == 2: 
        if Ub > mean([CurvatureDB[V3] - CurvatureDB[V2], CurvatureDB[V5] - CurvatureDB[V0], CurvatureDB[V4] - CurvatureDB[V1]]) > Lb:
#        if ((Ub > CurvatureDB[V3] - CurvatureDB[V2] >= Lb) or (Ub > CurvatureDB[V5] - CurvatureDB[V0] >= Lb))\
#         or ((Ub > CurvatureDB[V4] - CurvatureDB[V1] >= Lb)):
            return True
        else:
            return False
    elif Level == 3:
        if ((Ub > CurvatureDB[V3] - CurvatureDB[V2] >= Lb) or (Ub > CurvatureDB[V5] - CurvatureDB[V0] >= Lb))\
         or ((Ub > CurvatureDB[V4] - CurvatureDB[V1] >= Lb)):
            return True
        else:
            return False

    
def judgeFace3(FaceID, FaceDB, NbrLst, Curv, Threshold = 0.01, Level = 1):
    """Check whether a face satisfies the second-order criterion
        
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
                
        Level: integer
            There are many ways to judge whether a face can satisfy the 2nd-order derivative test, e.g., comparing each neighbor
            with the face or comparing the average value of all neighbors with the face. Level allows us to choose among different
            ways. 
            
            Level = 1: A face has to satisfy a condition with each of its neighbors/directions
            Level = 2: A face has to satisfy a condition with the average of all its neighbors/directions 
        
    """
        
    [V0, V1, V2] = FaceDB[FaceID]      
        
    FaceID2, V4 = neighborFaceVertex(FaceID, V0, V1, NbrLst, FaceDB)
    FaceID1, V3 = neighborFaceVertex(FaceID, V0, V2, NbrLst, FaceDB)
    FaceID3, V5 = neighborFaceVertex(FaceID, V1, V2, NbrLst, FaceDB)
    
    # The second-hop neighbors for V0-V2 are V6 and V7.
    Tmp, V6 = neighborFaceVertex(FaceID1, V0, V3, NbrLst, FaceDB)
#    if Tmp == -1:
#        return False
    Tmp, V7 = neighborFaceVertex(FaceID1, V2, V3, NbrLst, FaceDB)
#    if Tmp == -1:
#        return False
    
    # The second-hop neighbors for V1-V2 are V8 and V9.
    Tmp, V9 = neighborFaceVertex(FaceID3, V1, V5, NbrLst, FaceDB)
#    if Tmp == -1:
#        return False
    Tmp, V8 = neighborFaceVertex(FaceID3, V2, V5, NbrLst, FaceDB)
#    if Tmp == -1:
#        return False
    
    # The second-hop neighbors for V0-V1 are V10 and V11.
    Tmp, V10 = neighborFaceVertex(FaceID2, V1, V4, NbrLst, FaceDB)
#    if Tmp == -1:
#        return False
    Tmp, V11 = neighborFaceVertex(FaceID2, V0, V4, NbrLst, FaceDB)
#    if Tmp == -1:
#        return False
    
    D1 = (Curv[V6] - Curv[V2] + Curv[V7] - Curv[V0]) /2  - (Curv[V3] -Curv[V1])
    D2 = (Curv[V8] - Curv[V1] + Curv[V9] - Curv[V2]) /2  - (Curv[V5] -Curv[V0])
    D3 = (Curv[V11] - Curv[V1] + Curv[V10] - Curv[V0]) /2  - (Curv[V4] -Curv[V2])
    
    if Level == 1:    
        if D1 < Threshold and D2 < Threshold and D3 < Threshold:
            return True
        else:
            return False
    elif Level == 2:
        if mean([D1,D2,D3]) < Threshold :
            return True
        else:
            return False
    elif Level ==3:  # Level 3 is currently not active 
        if D1 < Threshold or D2 < Threshold or D3 < Threshold:
            return True
        else:
            return False
    elif Level == 4:
        if ( D1 < 0 and D2 < 0 ) \
        or ( D1 < 0 and D3 < 0 ) \
        or ( D2 < 0 and D3 < 0 ):
            return True
        else:
            return False

def faceFundi(FaceDB, CurvDB, NbrLst, CenterDb, Ratio=0.1, Skip2 = True, Skip3 = True):
    """Output strict fundus faces and candidate fundus faces 
    
    """
    
    print "Extracting fundi as face strips"
    
    Strict, Candidate = [], [] # to store IDs for strict and candidate faces 

    t1Lb, t2Lb, t2Ub, t3Lb = mean(CurvDB) + 0.5*std(CurvDB), -1*std(CurvDB)*Ratio, 1*std(CurvDB)*Ratio, -1*std(CurvDB)*Ratio
#    t1Lb, t2Lb, t2Ub, t3Lb = 0, -0.02, 0.05, 0.0045
    
    # empirical result: 
    # 1. We need to tolerate negative differences, other, the extracted fundi have holes, like a fish net. 
    # 2. t1Lb = 0.5std or 1mean+1std is good.
    # 3. if t1Lb, t2Lb, t2Ub, t3Lb = 0, -0.02, 0.05, 0.0045, set L2@4 and L3@1
        
    for i in xrange(0, len(FaceDB)):
#        print i
#        if judgeFace1(i, FaceDB, CurvDB, Threshold = t1Lb):
        if judgeFace1(i, CenterDb, Threshold = t1Lb):
            if Skip2:
                Candidate.append(i)
            else:
                if judgeFace2(i, FaceDB, NbrLst, CurvDB, Lb=t2Lb, Ub = t2Ub, Level = 2):
                    if Skip3: 
                        Candidate.append(i)
                    else:
                        if judgeFace3(i, FaceDB, NbrLst, CurvDB, Threshold = t3Lb, Level = 2):
                            Strict.append(i)
                        else:
                            Candidate.append(i)

                
    return Candidate, Strict

def skeletonFc(FundiList, NbrLst, MaxIter = 1):
    '''Skeletonize face strips, by iteratively removing vertexes that have a neighbor in the cloud but 'exposed' to the outer world.
    
    Input
    =====
        FundiList:    list of integers
            Each element is a face id. The vertex is considered as part of raw fundi, as face strips 
            
        NbrLst:    list of integers
            Each element contains ids of faces neighboring with a face
            
        MaxIter:    integer
            The maximum number of iteration allowed in Skeletonization by removing exposed faces. 
    
    '''
    
#    ToBeRemoved = [False for i in xrange(0, len(FundiList))]
    Iteration = 0
    while Iteration < MaxIter:
        ToBeRemoved = []
        Iteration += 1
    
        for i in xrange(0, len(FundiList)):
            Face = FundiList[i]
            Nbrs = NbrLst[Face][0]   
    #        HasCloudNbr, HasAirNbr = False, False  # air neighbor means a neighbor that is not in the cloud
            CloudNbr, AirNbr = [], []# neighbors in the cloud, and neighbors not in the cloud 
            for Nbr in Nbrs:
                if Nbr in FundiList:
    #                HasiNbr = True
                    CloudNbr.append(1)
                    AirNbr.append(0)
                else:
                    CloudNbr.append(0)
                    AirNbr.append(1)
            if sum(CloudNbr) > 1 and AirNbr > 1:   
                ToBeRemoved.append(Face)
        
        if ToBeRemoved == []:
            break
        else: 
            for Rmv in ToBeRemoved:
                FundiList.remove(Rmv)
            
    return FundiList

def fc2Vrtx(FaceDB, FaceList):
    '''Convert a list of faces into a list of vertexes that consist of those faces
    '''
    
    VrtxLst = []
    
    for Face in FaceList:  # obtain Face ID
        for Vrtx in FaceDB[Face]: # obtain 3 vertexes consisting the Face
            VrtxLst.append(Vrtx)
             
    return list(set(VrtxLst))  # The reason i use set() here is to remove repeated elements 

def fc2VrtxCmpnt(FaceDB, FaceList, FaceCmpnts):
    '''Convert a list of faces into a component-wise list of vertexes that consist of faces 
    '''
    print "Coverting fundus face strips into fundus vertex clouds"
    VrtxLst = [[] for i in xrange(0,len(FaceCmpnts))]
    
    for CmpntID, Cmpnt in enumerate(FaceCmpnts):
        for Face in Cmpnt:
            if Face in FaceList:  # obtain Face ID
                for Vrtx in FaceDB[Face]: # obtain 3 vertexes consisting the Face
                    VrtxLst[CmpntID].append(Vrtx)
             
    return [list(set(VrtxLst[i])) for i in xrange(0, len(FaceCmpnts))]

def myDT(VrtxClouds, VrtxNbrLst):
    '''Run my own distance transform, a node having neighbors not in the VrtxNbrLst has a distance 0 to the ``air'', o/w, 1  + min(distances of neighbors)
    
    Parameters
    =============
    
    VrtxClouds    : list of lists of integers
        each element is a list of vertexes that are fundi in one connected component 
        
    VrtxNbrLst    : list of list of integers
        each element is a list of ids of vertexes that are neighbors of a vertexes
    
    '''
    print "Performing distance transform on fundus vertex clouds obtained from fundus face strips"
    Map = []
    Large = 100 #  A large enough number
       
    for Cloud in VrtxClouds:
#        print " a new DT"
        Dist = [Large for i in xrange(0, len(Cloud))]

        if Dist == []:  # a component may have no fundus vertex 
            Map.append([])
            continue

        Queue = list(Cloud)  # elements in Queue are vertexes whose DT value has NOT been assigned yet.

        # pass 1: initialize vertexes that have DT value 0 
        while Queue != [] :
            Vrtx = Queue.pop()
            Nbrs = VrtxNbrLst[Vrtx]
            for Nbr in Nbrs:
                if not Nbr in Cloud: # a neighbor of Vrtx is an ``air'' vertex
                    Dist[Cloud.index(Vrtx)] = 0
                    break
        # iterative passes: calculate DT values for ``non-air'' vertexes  

# Commented 2011-05-01 01:34
#        Queue = list(Cloud)  # elements in Queue are vertexes whose DT value has NOT been assigned yet.
#        while Queue != []:
#            Vrtx = Queue.pop()
#            Nbrs = VrtxNbrLst[Vrtx]
#            if Dist[Cloud.index(Vrtx)] == Large: 
#                NbrCost = [Dist[Cloud.index(Nbr)] for Nbr in Nbrs]
#                LowestCost = min(NbrCost)
#                if LowestCost < Large:
#                    Dist[Cloud.index(Vrtx)] = 1 + LowestCost
#                else:  # no neighbor of Vrtx has yet been assigned 
#                    Queue.insert(0, Vrtx)  # put it back to be processed later
# End of commented 2011-05-01 01:34
# Block replaced by the block below
# New block, using gradually increasing DT value
        DT = 0
 #       print type(Dist), len(Dist)
        NotDone = True
        NewQueue = list(Cloud)
        while NotDone: # as long as there are unassigned vertex
            NotDone = False
            Queue = list(NewQueue)  # elements in Queue are vertexes whose DT value has NOT been assigned yet.
            NewQueue = [] 
#            print len(Queue)
            while Queue != []:
                Vrtx = Queue.pop()
                Nbrs = VrtxNbrLst[Vrtx]
                if Dist[Cloud.index(Vrtx)] == Large:
                    NbrCost = [Dist[Cloud.index(Nbr)] for Nbr in Nbrs]
                    LowestCost = min(NbrCost)
                    if LowestCost == DT:
                        Dist[Cloud.index(Vrtx)] = 1 + LowestCost
                    else:
                        NotDone = True
                        NewQueue.append(Vrtx)
            DT += 1
#        print "DT:", DT
# End of new block
        Map.append(Dist)
    return Map

# The function distFilt is obsolete on 2011-05-03 23:50. Part of it was replaced by dtCut() 
#def distFilt(VrtxClouds, Maps):
#    '''Remove vertexes according to their distance to the boundary of the fundus clouds
#    This function is current NOT in use. 2011-04-30 23:21
#    
#    Remove vertexes that are near the edge of the vertex clouds (thus, far from the medians of clouds). 
#    This function work component by component (two vertexes are in the same component if they are in the same basin). 
#    
#    Parameters
#    ===========
#    
#    VrtxClouds    : list of list of integers
#        each element is a list of ids of vertexes comprising fundi
#        
#        
#    Maps    : list of list of integers
#        each element is a list of distances of vertexes to the boundary of fundus clouds
#
#    Threshold    : integer
#        vertexes whose distance value below OR EQUAL TO Threshold will be removed. 
#     
#    '''
#    
#    for i in xrange(0, len(Maps)):  # i indicates the index of a component. The dimension of VrtxClouds is the same as the one of Map. 
#        Map, Cloud = Maps[i], VrtxClouds[i]
#        if len(Map) > 0: # we may have components that have no fundus vertex at all. 
#            ToBeRemoved = []  # vertexes that should be removed
#            if median(Map) > 2 :
#                Threshold = max(Map) - 2   
#            else:
#                Threshold = 1
##            Threshold = median(Map)  
#            for j in xrange(0, len(Map)):
#                if Map[j] <= Threshold:
#                    ToBeRemoved.append(Cloud[j])
#            for Vrtx in ToBeRemoved:
#                Cloud.remove(Vrtx) 
#    return VrtxClouds
# End of obsolete function distFilt() 2011-05-03 23:50 

def dtCut(Map):
    '''Determine the cut-off distance value for filtering vertexes in a component
    '''
# option 1: Using median and max     
#    if median(Map) >= 2 :
#        return max(Map) - 2
#    else:
#        return 1
# End of option 1

# option 2: Using top 10% of vertexes
#    Map.sort()
#    return Map[int(len(Map)*0.9)]
# End of option 2
    
# option 3: Using median
    return mean(Map) +2
# End of option 3

def covered(Seed, Radius, VrtxNbrLst):
    '''Use BFS to find all vertexes that is in a circle which centers at Seed with radius Radius
    
    Parameters 
    =============
    Radius    : integer 
        it is actually the depth of BFS
        
    Covered    : list of integers
        ids of vertexes that are within the range of Radius from vertex Seed
        
    Front : list of integers
        ids of vertexes that are in the front of BFS searching. It increases in each iteration from [].        
        
    Queue : list of integers
        ids of vertexes that are to be processed for this iteration of BFS. It decreases in each iteration to [].
    '''
    
    Covered, Queue  = [Seed], [Seed]
    
    while Radius > 0 :
        Front = []
        while Queue != []:
            Front += VrtxNbrLst[Queue.pop()] 
        Queue = list(set(Front))  # prepare queue for next iteration
        Covered += Queue 
        Radius -= 1
        
    return list(set(Covered))

def Skl(VrtxClouds, Maps, VrtxNbrLst):
    '''Skeletonize vertex clouds in my own algorithm from my own distance transformation
    
    This function is executed component by component (two vertexes are in the same component if they are in the same basin).
    The distance value of a vertex is
    the radius (Hamiltonian distance) of a circle that centers at the vertex and is tangent to the cloud, i.e., if the radius increase by 1,
    some vertexes inside the circle will not be fundus vertexes. 
    A *circle* from a vertex is a circle which centers at the vertex and of the radius the distance value of the vertex. 
    The skeletonization here is to find a set of vertexes, such that circles from those vertexes 
    can cover an entire (or most of) cloud within one component, for each component.  Well, necessary linking is needed. 
    
    For each component, sort distance values (obtained via myDT) of all vertexes. Label all vertexes as unvisited. 
    Begin from vertex X of the largest distance value. Call the vertex a seed. Label all vertexes that is within the circle from vertex X as visited and seeded from the vertex X. 
    Repeat this process on a vertex of the next largest distance value, until all vertexes are labeled visited and have recorded their seeds. 
    
    We can do line up later based on this result.       
    
    Parameters
    ===========
    
    VrtxClouds    : list of list of integers
        each element is a list of ids of vertexes comprising a fundus cloud in a component 
        
    VrtxNbrLst    : list of list of integers
        each element is a list of ids of vertexes that are neighbors of a vertexes
            
    Maps    : list of list of integers
        each element is a list of distances of vertexes to the boundary of fundus clouds    
    
    Seed    : list of integers 
        a list of seed vertexes that form the skeleton of a fundus cloud in one component
        
    Seeds    : list of list of integers
         each element is a Seed for one component      
    
    Skeleton    : list of list of integers
        each element is a list of ids of vertexes comprising the skeleton of a fundus cloud in a component 
    
    Dist    : a dictionary of integer keys and values 
        keys are vertex ids and values are their distance to the boundary of the cloud
    
    '''
    
    print "Skeletonizing DT map of fundus vertex clouds obtained from fundus face strips"
    
    Skeleton = []
    
    for i in xrange(0, len(Maps)):  # i indicates the index of a component. The dimension of VrtxClouds is the same as the one of Map.
#        print "Skeletonizing component:", i
        Map, Cloud = Maps[i], VrtxClouds[i]
        if len(Map) > 0: # we may have components that have no fundus vertex at all. 
            Seeds = [] # a list of  skeleton vertexes
         
            Dist = dict( [ (Cloud[j], Map[j]) for j in xrange(0, len(Map))] )
            Dist = sorted(Dist.iteritems(), key= lambda (k,v,) : (v,k))  

            CutoffDT = dtCut(Map)# 2011-05-03 23:48
            
            UnVisited = list(Cloud)
            while Dist != [] and UnVisited != []:
                Seed, Radius = Dist.pop()
                if Radius > CutoffDT:
                    Covered = covered(Seed, Radius, VrtxNbrLst)
                    HasNew = False
                    for Vtx in Covered:
                        if Vtx in UnVisited:
                            UnVisited.remove(Vtx)
                            HasNew = True
                    if HasNew:
                        Seeds.append(Seed)  # put Seed into Seeds if Seed covers vertexes that have not been covered
            if len(Seeds) >1:  # Added 2011-05-04 00:17 
                Skeleton.append(Seeds) # Added 2011-05-04 00:17
 
    return Skeleton

def faceToCurve(FaceDB, FaceList, FaceCmpnts, VrtxNbrLst):
    '''Extract fundus (curve) segments from fundus (curve) strips
    
    Parameters
    ===========
    
    FaceDB    : list of 3-tuple lists of integers
        each element of FaceDB are the 3 vertexes that consist of a face on the pial surface
        
    FaceList    : list of integers
        each element is the id of a face that is a fundus face 
        
    FaceCmpnts    : list of lists of integers 
        Each element is a list of ids of faces that are in the same basin component, and connected 
        
    
     
    '''
    
    Clouds = fc2VrtxCmpnt(FaceDB, FaceList, FaceCmpnts) # step 1: Convert faces in to Vrtx, clustered by components
    Map = myDT(Clouds, VrtxNbrLst)  # step 2: my own distance transformation
#    Clouds = distFilt(Clouds, Map) # step 3: filter before we do skeletonization
    Clouds = Skl(Clouds, Map, VrtxNbrLst)
    
#    return Map
    return Clouds

def isNbr(Face0, Face1, NbrDB): # copied from original fundi/link.py 2011-05-05 23:13
    '''Determining whether two faces, Face0 and Face1 are neighbors. Depending on my implementation, NbrDB or FaceDB may be provided. 
    '''
    if Face0 == Face1:
        return False   # a quick and easy way such that you don't have to check this when calling it in other functions.
    
    if Face1 in NbrDB[Face0]:
        return True 
    else:
        return False

def isolate(CandidateList, StrictList, NbrDB): # copied from original fundi/link.py 2011-05-05 23:13
    '''Denoise on candidate and strict set. 
    
    For candidate set: 
    Given a list of candidate faces and a table of neighbors of faces, 
    first remove faces who co-edge with NO strict face, 
    then remove those who co-edge with NO more than two neighbor in the candidate/strict set
    
    For strict set:
    If a strict set has more than two  
    
    '''
    NewCan, NewStc = [], []
    List = CandidateList + StrictList
    for Candidate in CandidateList:
        [F0, F1, F2] = NbrDB[Candidate][0]
        if F0 in StrictList or F1 in StrictList or F2 in StrictList:  # the face must first neighbor with a strict face
            if (F0 in List and F1 in List) or (F0 in List and F2 in List) or (F1 in List and F2 in List):
                NewCan.append(Candidate)
    List = NewCan + StrictList
    for Strict in StrictList:
        [F0, F1, F2] = NbrDB[Strict][0]
        if (F0 in List and F1 in List) or (F0 in List and F2 in List) or (F1 in List and F2 in List):
                NewStc.append(Strict)
    return NewCan, NewStc

def up1(Candidate, Strict, NbrDB): # copied from original fundi/link.py 2011-05-05 23:13
    '''Upgrading Candidate faces that share two out of three edges with strict faces into Strict sets
        ''' 
    Upgrade = [] # to store faces to be upgraded  
    Denoise = [] # strict faces to be removed as noises
    for Can in Candidate:
        [F0, F1, F2] = NbrDB[Can][0]
        if (F0 in Strict and F1 in Strict) or (F0 in Strict and F2 in Strict) or (F1 in Strict and F2 in Strict):
            Upgrade.append(Can)
#                print i, Can, NbrDB[Can]
    for S in Strict:
        [F0, F1, F2] = NbrDB[S][0]
        if not (F0 in Strict or F1 in Strict or F2 in Strict):
            Denoise.append(S) 
#    print len(Upgrade)     
    Strict += Upgrade
    for R in Upgrade:
        Candidate.remove(R)
    for R in Denoise:
        Strict.remove(R)
        
    return Candidate, Strict#, Upgrade