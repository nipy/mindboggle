#!/usr/bin/python
"""
Identify sulci from surface files with fold IDs and labels.


Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def generate_neighbor_lists(n_vertices, faces):
    """
    Generate the neighbor lists for all vertices.
    
    Parameters
    ------------
    n_vertices : Number of vertices in this hemisphere, an integer
    faces : faces in this hemisphere, [# faces x 3], list of lists of integers
     
    Returns
    --------
    neighbor_lists : neighbor lists for vertices, [#vertices x n],
            list of lists of integers;
            n is the number of neighbors for a vertex. n varies. 
     
    """
       
    print "Calculating vertex neighbor lists..."
    
    neighbor_lists = [[] for i in xrange(n_vertices)]
        
    for face in faces:
        [V0, V1, V2] = face
        
        if not V1 in neighbor_lists[V0]:
            neighbor_lists[V0].append(V1)
        if not V2 in neighbor_lists[V0]:
            neighbor_lists[V0].append(V2)

        if not V0 in neighbor_lists[V1]:
            neighbor_lists[V1].append(V0)
        if not V2 in neighbor_lists[V1]:
            neighbor_lists[V1].append(V2) 

        if not V0 in neighbor_lists[V2]:
            neighbor_lists[V2].append(V1)
        if not V1 in neighbor_lists[V2]:
            neighbor_lists[V2].append(V1)
    
    print "  neighbor list calculation done."
    
    return neighbor_lists

def identify(labels, folds, protocol_label_pair_lists, sulcus_indices,
             neighbor_lists, points):
    """
    Main function for identifying sulci from folds and labels.

    Definitions
    -----------

    A fold is a group of connected vertices that are deeper than a depth threshold.

    The ''label list'' is a list of all labels that help to define all sulci
    in the DKT protocol.

    A ''sulcus label list'' is the list of labels used to define a single sulcus.
    A sulcus label list is a subset of the label list.

    A ''sulcus label pair list'' contains some pairs of labels from a sulcus
    label list, where each pair defines a boundary between two labeled regions.
    No two sulcus label pair lists share a label pair.

    A ''sulcus ID'' or ''sulcus index'' identifies a sulcus. It is the index
    within the list of sulcus label pair lists for a given sulcus (starting from 1).

    Algorithm
    ---------

    1. First, assign -1 to all fold (i.e., non-zero) vertices with labels that are
       not in the label list.

    For each fold (vertices with the same non-zero label):

        2. If the set of labels in the fold is the same as in one of the sulcus
           label lists, assign the sulcus index for that list to all vertices
           in the fold.

        3. If the set of labels in the fold is a subset of labels in only one
           of the sulcus label lists, assign the sulcus index for that list
           to all vertices in the fold.

        4. If the labels in the fold are a subset of the labels in more than one
           of the sulcus label lists (e.g., the labels in the fold are [1,2,3]
           and we have two label sets [1,2,3,4] and [1,2,3,5] in sulcus label lists):

           (a) Select the sulcus label list, if any, which is the only one
               that has two of the fold labels as one of its pairs.

           (b) If there is more than one such label list, skip to #6.

        5. If the labels in the fold are a superset of the labels for only one
           of the sulcus label lists (ex: the labels in the fold are [1,2,3,4]
           and there is only one sulcus label list which is its subset [1,2,3]):

           (a) Assign the sulcus index for that list to just the vertices
               in the fold with labels in its list (vertices with label 1, 2 or 3).

           (b) Treat the remaining vertices (vertices of label 4) as a new fold
               and start from #2.

        6. If the labels in the fold are a superset of the labels in more than 
           one of the sulcus label lists (e.g., labels in the fold are [1,2,3,4] 
           and there are two sulcus label lists [1,2,3] and [2,3,4]), or if the 
           assignment is ambiguous (as in 4b), perform segmentation:

           (a) Find all label boundaries within the fold (vertices with 
               two different labels in its neighborhood of connected vertices).

           (b) Find label boundary vertices whose two different labels 
               (A) are and (B) are not a sulcus label pair (in the protocol). 
               ''Two different labels'' means a label pair.

           (c) For each label in A:
               If there is only one label pair P in the fold containing that label,
               assign all vertices with that label the sulcus index of the 
               sulcus pair list (corr. to the boundary) that contains that label. 
               Otherwise, skip to next label.

           (d) For each remaining vertex (with label in A or B):
               Assign the vertex the sulcus index of the nearest A boundary.

    Parameters
    ----------
    labels : labels for all vertices (list of integers, 
        [# vertices x 1],  indexed from 1)
    folds : lists of vertices in folds, list of lists of integers
    label_pairs : list of label pairs for each sulcus, all contained within a list
    sulcus_indices : sulcus indices assigned or to be assigned to all vertices
        [# vertices x 1], indexed from 1
    points : coordinates of all vertices on the surface [# vertices x 3]
    neighbor_lists : neighbor lists of all vertices, [#vertices x n],
        list of lists of integers.
        n is the number of neighbors for a vertex. n varies.

    Returns
    -------
    sulcus_indices : sulcus indices for all vertices (a list of integers,
        [# vertices x 1], indexed from 1, with zeros for non-sulcus vertices)
                  
    """
    import numpy as np

    #---------------------------------------------------------------------------
    # Nested function definitions
    #---------------------------------------------------------------------------
    def find_superset_lists(labels, sulcus_label_lists):
        """
        Find lists within *sulcus_label_lists* that are supersets 
        of a list of *labels*.
        
        Example
        -------
        >>> find_superset_lists([1,2],[[1,2],[3,4]])
        >>> [0]
        >>> find_superset_lists([1,2],[[1,2,3],[1,2,5]])
        >>> [0, 1]
        >>> find_superset_lists([1,2],[[2,3],[1,2,5]])
        >>> [1]
        
        """
        labels = set(labels)
        subset_indices = []
        for Id, pair_list in enumerate(sulcus_label_lists):
            if labels.issubset(set(pair_list)):
                subset_indices.append(Id)
                
        return subset_indices
                
    def find_subset_pairs(fold_labels, protocol_label_pair_lists, selected_indices=[]):
        """
        Find label pairs within *protocol_label_pair_lists* 
        that are subsets of *fold_labels*.

        Parameters
        ----------
        selected_indices : Indices (0-index list of integers) of 
                           protocol_label_pair_lists we should check
        
        Example
        -------
        >>> find_subset_pairs([1,2,3],[[[1,2]]])
         [0]
        >>> find_subset_pairs([1,2,3],[[[1,2],[3,4]], [[1,3]], [[1,5],[4,3]]])
         [0, 1]
         
        """
        
        if selected_indices==[]:
            selected_indices = range(len(protocol_label_pair_lists))
        
        fold_labels = set(fold_labels)
        subset_indices = []
        for i in selected_indices:
            row = protocol_label_pair_lists[i]
            for pair in row:
                if set(pair).issubset(fold_labels):
                    subset_indices.append(i)
        
        return subset_indices

    def find_subset_lists(labels, sulcus_label_lists):
        """
        Find lists within *sulcus_label_lists* that are subsets of *labels*.
        
        Example
        ---------
        >>> find_subset_lists([1,2,3,4],[[1,2],[2,3,4]])
         [0, 1]     
        
        """
        labels = set(labels)
        superset_index = []
        for Id, pair_list in enumerate(sulcus_label_lists):
            if set(pair_list).issubset(labels):
                superset_index.append(Id)
                
        return superset_index
    
    def detect_boundaries(labels, fold, neighbor_lists):
        """
        Detect the label boundaries in a fold.
        
        Parameters
        ----------
        labels :  labels for all vertices (list of integers, 
             [# vertices in a hemisphere x 1],  indexed from 1)
        fold : a list of integers, representing the vertices in this fold 
        neighbor_lists : neighbor lists of all vertices, [#vertices x n],
            list of lists of integers.
            n is the number of neighbors for a vertex. n varies.

        Returns
        -------
        boundary_vertices : a list of integers, representing vertices at
                            label boundaries
        boundary_labels : a dictionary whose keys are vertex IDs
                          and values are lists of integers
        
        """
        boundary_vertices = []
        boundary_labels = {}
        for vertex in fold:
            labels_of_neighbors = [labels[a_neighbor]
                                   for a_neighbor in neighbor_lists[vertex]
                                   if a_neighbor in fold]
            if not all(x == labels_of_neighbors[0] for x in labels_of_neighbors): 
                boundary_vertices.append(vertex)
                boundary_labels[vertex] = labels_of_neighbors 
            
        return boundary_vertices, boundary_labels
    
    def only_one_pair(a_label, a_label_pair_list):
        """
        Check whether a_label only appears in one pair in a_label_pair_list.
        
        If it does, return such a pair. Otherwise, an empty list.
        """
        appeared = False
        for pair in a_label_pair_list:
            if a_label in pair:
                if appeared: 
                    return False, []
                else:
                    appeared = True
        return True, pair
    
    def euclidean_distance(point1, point2):
        """
        Estimate the Euclidean distance between two points.
        """
        from numpy import sqrt
        return sqrt((point1[0] - point2[0]) ** 2 + \
                    (point1[1] - point2[1]) ** 2 + \
                    (point1[2] - point2[2]) ** 2)

    #---------------------------------------------------------------------------
    # Prepare data structures
    #---------------------------------------------------------------------------
    # Prepare list of sulcus label lists (labels within a sulcus)
    sulcus_label_lists = []
    for row in protocol_label_pair_lists:
        sulcus_label_lists.append(list(np.unique(np.asarray(row))))

    # Prepare list of all unique labels in the labeling protocol
    label_list = []
    [label_list.append(x)
     for lst in sulcus_label_lists
     for x in lst
     if x not in label_list]

    # Prepare list of all unique sorted label pairs in the labeling protocol
    protocol_pairs = []
    [protocol_pairs.append(list(set(x)))
     for lst in protocol_label_pair_lists
     for x in lst
     if list(set(x)) not in protocol_pairs]

    #---------------------------------------------------------------------------
    # Loop through folds
    #---------------------------------------------------------------------------
    for fold in folds:
        labels_in_this_fold = list(set([labels[vertex] for vertex in fold]))
        
        # Case 1: If the fold has only one label, assign sulcus index of -1
        # to all vertices in this fold
        if len(labels_in_this_fold) == 1:
            continue  # sulcus_indices already initialized with -1 values
        
        # Case 2: If the set of labels in the fold is the same as in one 
        # of the sulcus label lists, assign the sulcus index for that 
        # list to all vertices in the fold.
        elif labels_in_this_fold in sulcus_label_lists:
            # Add 1 because the sulcus index begins from 1
            Index = sulcus_label_lists.index(labels_in_this_fold) + 1
            for vertex in fold:
                sulcus_indices[vertex] = Index
                
        else: # Cases 3-6
            superset_list_indices = find_superset_lists(labels_in_this_fold,
                                                        sulcus_label_lists)
            if len(superset_list_indices) == 1:
                # Case 3: If the set of labels in the fold is a subset of labels
                # in only one of the sulcus label lists, assign the sulcus index 
                # for that list to all vertices in the fold.
                for vertex in fold:
                    sulcus_indices[vertex] = superset_list_indices[0] + 1
            elif len(superset_list_indices) > 1:
                # Case 4: If the labels in the fold are a subset of the labels
                # in more than one of the sulcus label lists:
                # select the sulcus label list, if any, which is the only
                # one that has two of the fold labels as one of its pairs.
                subset_pair_indices = find_subset_pairs(labels_in_this_fold,
                                          protocol_label_pair_lists,
                                          selected_indices=superset_list_indices)
                if len(subset_pair_indices) == 1:
                    for vertex in fold: 
                        sulcus_indices[vertex] = subset_pair_indices[0] + 1
                else:
                    pass
            else: # Cases 5-6
                subset_list_indices = find_subset_lists(labels_in_this_fold,
                                                        sulcus_label_lists)
                # Case 5: If the labels in the fold are a superset of the labels
                # for only one of the sulcus label lists
                # (ex: The labels in the fold are [1,2,3,4] and there is
                #      only one sulcus label list which is its subset [1,2,3]):
                # (a) Assign the sulcus index for that list to just the vertices
                #     in the fold with labels in its list (1, 2, or 3).
                # (b) Treat the remaining vertices (vertices of label 4)
                #     as a new fold and start from #1.
                if len(subset_list_indices) == 1: # Case 5
                    unassigned = [[]]
                    for vertex in fold:
                        if labels[vertex] in sulcus_label_lists[subset_list_indices[0]]:
                            # subcase a
                            sulcus_indices[vertex] = subset_list_indices[0] + 1
                        else:
                            # subcase b
                            unassigned[0].append(vertex)                    
                            sulcus_indices = identify(labels, unassigned,
                                protocol_label_pair_lists,
                                sulcus_indices, neighbor_lists, points)
                else: # Case 6
                    print "  case 6..."
                    # step (a): Detect label boundary vertices
                    boundary_vertices, boundary_labels = detect_boundaries(labels,
                        fold, neighbor_lists)
                    
                    # step (b): Find label boundary vertices whose two different
                    # labels (A) and (B) are not a label pair in the protocol.
                    A, B = [], []
                    for vertex in boundary_vertices:
                        vertex_labels = boundary_labels[vertex]
                        if vertex_labels not in protocol_pairs:
                            if not vertex_labels in A:
                                A.append(vertex_labels)

                    # step (c): For each label in A:
                    # If there is only one label pair P in the fold containing
                    # that label, assign all vertices with that label the sulcus
                    # index of the sulcus pair list (corr. to the boundary)
                    # that contains that label. Otherwise, skip to next label.
                    All_labels_in_A = []
                    for pair in A:
                        if not pair[0] in All_labels_in_A:
                            A.append(pair[0])
                        if not pair[1] in All_labels_in_A:
                            A.append(pair[1])
                            
                    # We create a new data structure to store vertices
                    # that are assigned sulcus indices
                    vertices_Assigned_sulcus_indices = []
                    
                    for Label in All_labels_in_A:
                        Appear_in_one_pair_Only, Such_pair = only_one_pair(Label, A)
                        if Appear_in_one_pair_Only:
                            if Such_pair == []:
                                print "Error"
                                exit()
                            for rowID, row in enumerate(protocol_label_pair_lists):
                                if Such_pair in row:
                                    Assign_Index = rowID + 1
                        # Assign the proper index to all vertices in this fold bearing this label
                        for vertex in fold:
                            if labels[vertex] == Label:
                                sulcus_indices[vertex] = Assign_Index
                                vertices_Assigned_sulcus_indices.append(vertex)      
                        
                    # step (d): For each remaining vertex (with label in A or B):
                    # Assign the vertex the sulcus index of the nearest A boundary.
                    
                    # Only do (d) on vertices in this fold but not in vertices_Assigned_sulcus_indices
                    # Please note that vertices_Assigned_sulcus_indices will NOT be updated from now
                    for vertex in fold:
                        if vertex not in vertices_Assigned_sulcus_indices:
                            # Estimate the distance from this vertex to all label boundary vertices
                            # And find the closest boundary vertex in Euclidean distance
                            Shortest_Distance = 1000000
                            Closest_boundary_vertex = -1
                            
                            for boundary_vertex in boundary_vertices: # Note the plural form
                                Distance = euclidean_distance(points[boundary_vertex], points[vertex]) 
                                if  Distance < Shortest_Distance:
                                    Closest_boundary_vertex = boundary_vertex
                                    Shortest_Distance = Distance
                                     
                            # Assign the sulcus index of the closest 
                            sulcus_indices[vertex] = sulcus_indices[Closest_boundary_vertex]
                    
    return sulcus_indices
    

if __name__ == "__main__":

    import sys
    from mindboggle.utils import io_vtk

    # Load labels and folds (the second surface has to be inflated).
    points, faces, labels = io_vtk.load_scalar(sys.argv[1], return_arrays=0)
    points, faces, fold_IDs = io_vtk.load_scalar(sys.argv[2], return_arrays=0)

    # Calculate neighbor lists for all vertices
    neighbor_lists = generate_neighbor_lists(len(points), faces)

    protocol_label_pair_lists = [[[28,12]],\
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
    
    # Create a list of lists of folds
    unique_fold_IDs = set(fold_IDs) # we later traverse each fold
    folds = []
    for i in unique_fold_IDs: # on every fold
        fold = [vertex for vertex in xrange(len(labels)) if fold_IDs[vertex] == i]
        folds.append(fold)
    
    # Assign a sulcus index to each fold vertex.
    # Initialize all fold vertices as -1.
    sulcus_indices = [-1 for i in xrange(len(labels))]
    # Since we do not touch gyral vertices and vertices whose labels
    # are not in label list, or vertices having only one label,
    # their sulcus indices will remain -1.
    sulcus_indices = identify(labels, folds, sulcus_label_pair_list, \
        sulcus_indices, neighbor_lists, points)

    # Finally, write points, faces and sulcus_indices to a new vtk file
    io_vtk.write_scalars(sys.argv[3], points, range(len(points)), faces, \
        LUTs=[sulcus_indices], LUT_names=['Sulcus indices'])
    