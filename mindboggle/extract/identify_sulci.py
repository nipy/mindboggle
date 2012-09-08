# This script loads a labeled surface, and identify sulci within folds. 
#
#      
#
#Sulcus label pair lists:  (sulcus label lists are the unions of each of these)
#
#label_pairs = [[[28,12]],
#
#[[28,3]],
#
#[[3,18]],
#
#[[28,24], [3,24], [18,24]],
#
#[[24, 22]],
#
#[[22,29], [22,31]],
#
#[[29,31], [29,8]],
#
#[[31,8]],
#
#[[31,30]],
#
#[[11,8], [11,29]],
#
#[[11,15], [11,9]],
#
#[[30,15]],
#
#[[15,9]],
#
#[[35,30], [35,34], [35,12], [35,2], [35,24], [35,22], [35,31]],
#
#[[30,34]],
#
#[[2,14], [2,28], [2,17], [25,17]],
#
#[[17,28]],
#
#[[5,25]],
#
#[[13,25], [13,2]],
#
#[[28,14]],
#
#[[2,4]],
#
#[[12,18], [12,3]],
#
#[[12,14]],
#
#[[7,9], [7,11]],
#
#[[7,6], [7,16], [7,13]]]
#
#for pairs in label_pairs:
#
#   upairs = np.unique(np.asarray(pairs))
#
#   if len(upairs) > 2:
#
#       print(upairs)
#
#three_label_lists = [
#
#[22,29,31],
#
#[8,29,31],
#
#[8,11,29],
#
#[9,11,15],
#
#[2,13,25],
#
#[3,12,18],
#
#[7,9,11]]
#
#four_label_lists = [
#
#[3,18,24,28],
#
#[6,7,13,16]]
#
#five_label_lists = [
#
#[2,14,17,25,28]]
#
#eight_label_list = [
#[2,12,22,24,30,31,34,35]]


def identify(Labels, Folds, Sulcus_Label_Pair_Lists, Sulcus_Index, Neighbor):
    """The main function
    Parameters
    ------------

    Labels:  labels for all vertices (list of integers, 
             [# vertexes x 1],  indexed from 1)
    Folds: Lists of vertexes in folds, list of lists of integers
    Label_Pairs: list of label pairs for each sulcus, all contained within a list
    Sulcus_Index: sulcus indexes assigned  or to be assigned to all vertexes [# vertexes x 1], indexed from 1

    Returns
    ----------
    Sulcus_Index: sulcus indexes for all vertices (a list of integers, [# vertexes x 1],
                  indexed from 1, with zeros for non-sulcus vertices)
                  
    Notes
    -------
    
    A fold is a group of connected vertices that are deeper than a depth threshold.

    The ''label list'' is a list of all labels that help to define all sulci in the DKT protocol.
    
    A ''sulcus label list'' is the list of labels used to define a single sulcus.
    A sulcus label list is a subset of the label list.

    A ''sulcus label pair list'' contains some pairs of labels from a sulcus label list, 
    where each pair defines a boundary between two labeled regions. 
    No two sulcus label pair lists share a label pair.

    A ''sulcus ID'' or ''sulcus index'' identifies a sulcus.  
    It is the index within the list of sulcus label pair lists for a given sulcus (starting from 1).

    Here is the algorithm. 

    1. Assign -1 to all fold (i.e., non-zero) vertices with labels that are not in the label list.

    For each fold (vertices with the same non-zero label):
    
        2. If the set of labels in the fold is the same as in one of the sulcus label lists, 
           assign the sulcus index for that list to all vertices in the fold.
    
        3. If the set of labels in the fold is a subset of labels in only one of the sulcus label lists, 
           assign the sulcus index for that list to all vertices in the fold.
        
        4. If the labels in the fold are a subset of the labels in more than one of the sulcus label lists 
           (e.g., the labels in the fold are [1,2,3] and we have two label sets [1,2,3,4] and [1,2,3,5] 
           in sulcus label lists):
        
          (a) Select the sulcus label list, if any, which is the only one that has two of the fold labels 
          as one of its pairs.
        
          (b) If there is more than one such label list, skip to #6.
        
        5. If the labels in the fold are a superset of the labels for only one of the sulcus label lists (e.g., the labels in the fold are 
        [1,2,3,4] and there is only one sulcus label list which is its subset [1,2,3]):
        
          (a) Assign the sulcus index for that list to just the vertices in the fold with labels in its list 
          (i.e., vertices with label 1, 2 or 3).
        
          (b) Treat the remaining vertices (i.e., vertices of label 4) as a new fold and start from #1.
        
        6. If the labels in the fold are a superset of the labels in more than one of the
           sulcus label lists (e.g., labels in the fold are [1,2,3,4] and there are two sulcus 
           label lists [1,2,3] and [2,3,4]), or if the assignment is ambiguous (as in 4b), 
           perform segmentation:
        
          (a) Find all label boundaries within the fold
               (vertices with two different labels in its neighborhood of connected vertices).
        
          (b) Find label boundary vertices whose two different labels are not a sulcus label pair.
        
          (c) Treat each non-sulcus boundary vertex as a seed, and propagate an assignment
               of -2 until the termination criterion is met for every new seed:
                   -- the nearest label to the seed (that is in the label list or zero,
                      and is different than the  seed's label) is not one of the labels
                      of the original seeds.
        
           (d) Reassign each -2 vertex to its nearest label in the label list.
           (e) Return to #1.
                  
    """
    
    # Nested function definitions 
    def check_list_subset(Labels, Sulcus_Label_Lists):
        """Check if a list of labels, *Labels*, is a subset or subsets of a member list of *Sulcus_Label_Lists*
        
        Example
        ---------
        >>> check_list_subset([1,2],[[1,2],[3,4]])
        >>> [0]
        >>> check_list_subset([1,2],[[1,2,3],[1,2,5]])
        >>> [0, 1]
        >>> check_list_subset([1,2],[[2,3],[1,2,5]])
        >>> [1]

        
        
        """
#        for i, Pair_List in enumerate(Sulcus_Label_Pair_Lists):
#            if Pair in Pair_List:
#                return i
#        return -1
        Labels = set(Labels) 
        Subset_Index = []
        for Id, Pair_List in enumerate(Sulcus_Label_Lists):
            if Labels.issubset(set(Pair_List)):
                Subset_Index.append(Id)
                
        return Subset_Index
                
    def check_subset_pair(Fold_Labels, Sulcus_Label_Pair_Lists, Selected_Index=[]):
        """Check how many lists in sulcus label pair lists have its pair as subset of Fold_Labels. 
        
        Parameters
        -----------
        
        Selected_Index : Index (0-index) of Sulcus_Label_Pair_Lists that we should check, list of integers
        
        Note
        -----
        
        A fold can have more than 2 labels. 
          
        Example
        ----------
        
        >>> check_subset_pair([1,2,3],[[[1,2]]])
         [0]
        
        >>> check_subset_pair([1,2,3],[[[1,2],[3,4]], [[1,3]], [[1,5],[4,3]]])
         [0, 1]

         
        """
        
        if Selected_Index==[]:
            Selected_Index = range(len(Sulcus_Label_Pair_Lists))
        
        Fold_Labels = set(Fold_Labels)
        Subset_Pair_Index = []
        for i in Selected_Index:
            Row = Sulcus_Label_Pair_Lists[i]
            for Pair in Row:
                if set(Pair).issubset(Fold_Labels):
                    Subset_Index.append(i)
        
        return Subset_Pair_Index 

    def check_list_superset(Labels, Sulcus_Label_Lists):
        """Check if a list of labels, *Labels*, is a subset or subsets of a member list of *Sulcus_Label_Lists*
        
        Example
        ---------
        >>> check_list_superset([1,2,3,4],[[1,2],[2,3,4]])
         [0, 1]

        
        
        """
#        for i, Pair_List in enumerate(Sulcus_Label_Pair_Lists):
#            if Pair in Pair_List:
#                return i
#        return -1
        Labels = set(Labels)
        Superset_Index = []
        for Id, Pair_List in enumerate(Sulcus_Label_Lists):
            if set(Pair_List).issubset(Labels):
                Superset_Index.append(Id)
                
        return Superset_Index
    
    def detect_boundary(Labels, Fold, Neighbor):
        """Detect the label boundary in a fold
        
        Parameters
        -------------
        Labels:  labels for all vertices (list of integers, 
             [# vertexes in a hemisphere x 1],  indexed from 1)
        Fold: a list of integers, representing the vertexes in this fold 
        Neighbor: neighbors of each vertex in a hemisphere
                 (list of lists of integers, [# vertexes in a hemisphere x 3], indexed from 0) 
        
        Returns
        --------
        Boundary: a list of integers, representing vertexes at label boundaries. 
        Labels_of_Boundary_Vertexes: a dictionary, whose keys are vertexes IDs and values are lists of integers
        
        """
        Boundary_Vertexes = []
        Labels_of_Boundary_Vertexes = {}
        for Vertex in Fold:
            Labels_of_Neighbors = [Labels[A_Neighbor] for A_Neighbor in Neighbor[Vertex] if A_Neighbor in Fold]
            if not all(x == Labels_of_Neighbors[0] for x in Labels_of_Neighbors): 
                Boundary_Vertexes.append(Vertex)
                Labels_of_Boundary_Vertexes[Vertex] = Labels_of_Neighbors 
            
        return Boundary_Vertexes, Labels_of_Boundary_Vertexes
    
    def only_one_pair(A_Label, A_Label_Pair_List):
        """Check whether A_Label only appear in one pair in A_Label_Pair_List
        
        if it does, return such a pair. O/w, an empty list. 
        """
        Appeared = False
        for Pair in A_Label_Pair_List:
            if A_Label in Pair:
                if Appeared: 
                    return False, []
                else:
                    Appeared = True
        return True, Pair
    
    # End of nested function definitions 
        
    import numpy as np
    # Prepared data structures 
    Sulcus_Label_Lists = []
    for Row in Sulcus_Label_Pair_Lists:
        Sulcus_Label_Lists.append(list(np.unique(np.asarray(Row))))
        
    All_Protocol_Pairs = [] 
    for Row in Sulcus_Label_Pair_Lists:
        for Pair in Row:
            if not Row in All_Protocol_Pairs:
                All_Protocol_Pairs.append(Pair)
     
    Label_List = []
    for Row in Sulcus_Label_Lists:
        Label_List += Row
    Label_List = list(set(Label_List))
    # End of prepare data structures 
    
    for Fold in Folds:          
        Labels_in_this_Fold = set([Labels[Vertex] for Vertex in Fold])
        
        # Case 1: If the fold has only one label, let sulcus index of all vertexes in this fold be -1
        if len(Labels_in_this_Fold) == 1:
            continue  # Sulcus_Index is initialized with -1's
        
        # Case 2: If the set of labels in the fold is the same as in one 
        # of the sulcus label lists, assign the sulcus index for that 
        # list to all vertices in the fold.
        elif Labels_in_this_Fold in Sulcus_Label_Lists:
            Index = Sulcus_Label_Lists.index(Labels_in_this_Fold) + 1 # The plus 1 is because sulcus index begins from 1. 
            for Vertex in Fold:
                Sulcus_Index[Vertex] = Index
                
        else: # Cases 3-6
            Subset_Index = check_list_subset(Labels_in_this_Fold, Sulcus_Label_Lists)
            if len(Subset_Index) == 1:
                # Case 3: 3. If the set of labels in the fold is a subset of labels 
                # in only one of the sulcus label lists, assign the sulcus index 
                # for that list to all vertices in the fold.
                for Vertex in Fold:
                    Sulcus_Index[Vertex] = Subset_Index[0] + 1
            elif len(Subset_Index) > 1: 
                # Case 4: If the labels in the fold are a subset of the labels in more 
                # than one of the sulcus label lists:
                # (a) Select the sulcus label list, if any, which is the only one that has two of 
                # the fold labels as one of its pairs. 
                Subset_Pair_Index = check_subset_pair(Labels_in_this_Fold, Sulcus_Label_Pair_Lists, Selected_Index=Subset_Index)
                if len(Subset_Pair_Index) ==1:
                    for Vertex in Fold: 
                        Sulcus_Index[Vertex] = Subset_Pair_Index[0] + 1
            else: # Cases 5-6
                Superset_Index = check_list_superset(Labels_in_this_Fold, Sulcus_Label_Lists)
                # Case 5. If the labels in the fold are a superset of the labels for only one 
                # of the sulcus label lists (e.g., the labels in the fold are [1,2,3,4] 
                # and there is only one sulcus label list which is its subset [1,2,3]):
                # (a) Assign the sulcus index for that list to just the vertices in the 
                # fold with labels in its list (i.e., vertices with label 1, 2 or 3).
                # (b) Treat the remaining vertices (i.e., vertices of label 4) as a new fold and start from #1.
                if len(Superset_Index) == 1: # Case 5
                    SubFolds = [[]]
                    for Vertex in Fold:
                        if Labels[Vertex] in Sulcus_Label_Lists[Superset_Index[0]]: #subcase a
                            Sulcus_Index[Vertex] = Superset_Index[0] + 1 
                        else: 
                            SubFolds[0].append(Vertex)                    
                            Sulcus_Index = identify(Labels, SubFolds, Sulcus_Label_Pair_Lists, Sulcus_Index)
                else: # Cases 6 to be done. 
                    print "case 6!"
                    # step (a): Detect label boundary vertexes
                    Boundary_Vertexes, Labels_of_Boundary_Vertexes = detect_boundary(Labels, Fold, Neighbor)
                    
                    # step (b): Find label boundary vertices whose two different labels (A) are 
                    # and (B) are not a sulcus label pair (in the protocol).  
                    A, B = [], []  # The A and B as we defined in the pseudocode 
                    for Vertex in Boundary_Vertexes:
                        Label_Pair_at_this_Boundary_Vertex = Labels_of_Boundary_Vertexes[Vertex]
                        if Label_Pair_at_this_Boundary_Vertex in All_Protocol_Pairs:
                            if not Label_Pair_at_this_Boundary_Vertex in A:
                                A.append(Label_Pair_at_this_Boundary_Vertex)
                        else:
                            if not Label_Pair_at_this_Boundary_Vertex in B:
                                B.append(Label_Pair_at_this_Boundary_Vertex)
                                
                    # step (c): For each label in A:
                    # If there is only one label pair P in the fold containing that label,
                    # assign all vertices with that label the sulcus index of the sulcus pair list 
                    # (corr. to the boundary) that contains that label. Otherwise, skip to next label.
                    All_Labels_in_A = []
                    for Pair in A:
                        if not Pair[0] in All_Labels_in_A:
                            A.append(Pair[0])
                        if not Pair[1] in All_Labels_in_A:
                            A.append(Pair[1])
                    
                    for Label in All_Labels_in_A:
                        Appear_in_one_Pair_Only, Such_Pair = only_one_pair(Label, A)
                        if Appear_in_one_Pair_Only:
                            if Such_Pair == []:
                                print "Error"
                                exit()
                            for RowID, Row in enumerate(Sulcus_Label_Pair_Lists):
                                if Such_Pair in Row:
                                    Assign_Index = RowID + 1
                        # Assign the proper index to all vertexes in this fold bearing this label
                        for Vertex in Fold:
                            if Labels[Vertex] == Label:
                                Sulcus_Index[Vertex] = Assign_Index       
                        
                    # step (d): For each remaining vertex (with label in A or B):
                    # Assign the vertex the sulcus index of the nearest A boundary.
                    
                    # To be done
                    
    return Sulcus_Index 
    

if __name__ == "__main__":
    from mindboggle.utils import io_vtk
    import sys
    Points, Faces, Labels = io_vtk.load_scalar(sys.argv[1], return_arrays=0)
    Points, Faces, Fold_IDs = io_vtk.load_scalar(sys.argv[2], return_arrays=0)  
    # the second surface has to be inflated. 
    
    Sulcus_Label_Pair_List = [[[28,12]],\
                   [[28,3]],\
                   [[3,18]],\
                   [[28,24], [3,24], [18,24]],\
                   [[24, 22]],\
                   [[22,29], [22,31]],\
                   [[29,31], [29,8]],\
                   [[31,8]],\
                   [[31,30]],\
                   [[11,8], [11,29]],\
                   [[11,15], [11,9]],\
                   [[30,15]],\
                   [[15,9]],\
                   [[35,30], [35,34], [35,12], [35,2], [35,24], [35,22], [35,31]],\
                   [[30,34]],\
                   [[2,14], [2,28], [2,17], [25,17]],\
                   [[17,28]],\
                   [[5,25]],\
                   [[13,25], [13,2]],\
                   [[28,14]],\
                   [[2,4]],\
                   [[12,18], [12,3]],\
                   [[12,14]],\
                   [[7,9], [7,11]],\
                   [[7,6], [7,16], [7,13]]]
    
    Unique_Fold_ID = set(Fold_IDs) # we later traverse each fold
    
    Folds = []
    for i in Unique_Fold_ID: # on every fold
        Fold = [Vertex for Vertex in xrange(len(Labels)) if Fold_IDs[Vertex] == i]
        Folds.append(Fold)
        
    Sulcus_Index = [-1 for i in xrange(len(Labels))]  # Initially -1 for all vertexes
    # Since we do not touch gyral vertexes and vertexes whose labels are not in label list, 
    # or vertexes having only one label, their sulcus indexes will remain -1. 
    
    Sulcus_Index = identify(Labels, Folds, Sulcus_Label_Pair_List, Sulcus_Index, Neighbor)
    # The Neighbor should be loaded from the 3rd file
    
    # finally, write points, faces and sulcus_index back to a new vtk file 
    