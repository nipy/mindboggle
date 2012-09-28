#!/usr/bin/python
"""
Identify sulci from surface files with fold IDs and labels.


Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import sys
#import numpy as np
#from utils.io_vtk import load_scalar, write_scalars
#from info.sulcus_boundaries import sulcus_boundaries
#from utils.mesh_operations import segment, detect_boundaries, \
#    find_neighbors, compute_distance

verbose = 1


def identify(labels, folds, label_pair_lists, sulcus_IDs,
             neighbor_lists, points):
    """
    Identify sulcus folds in a brain surface according to a labeling protocol
    that includes a list of label pairs defining each sulcus.

    Parameters
    ----------
    labels : list of integers, indexed from 1
        labels for all vertices
    folds : list of lists of integers
        each list contains indices to vertices of a fold
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus
    sulcus_IDs : list of integers, indexed from 1
        sulcus IDs assigned or to be assigned to all vertices
    points : list of lists of three floats
        coordinates of all vertices on the surface
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    sulcus_IDs : list of integers, indexed from 1
        sulcus IDs for all vertices, with -1s for non-sulcus vertices

    Definitions
    -----------

    A fold is a group of connected vertices that are deeper than a depth threshold.

    The ''label list'' is a list of all labels that help to define all sulci
    in the DKT protocol.

    The ''label pairs list'' contains all pairs of labels from the protocol.

    A ''sulcus label list'' is the list of labels used to define a single sulcus.
    A sulcus label list is a subset of the label list.

    A ''sulcus label pair list'' contains some pairs of labels from a sulcus
    label list, where each pair defines a boundary between two labeled regions.
    No two sulcus label pair lists share a label pair.

    A ''sulcus ID'' uniquely identifies a sulcus.
    It is the index to a particular sulcus label list (or sulcus label pair list)
    in its parent list (indices start from 1).

    Algorithm
    ---------

    First, remove all fold vertices with labels that are not in the label list
    by assigning -1 to them.

    For each fold (vertices with the same non-zero label):

        Case 1:

        If the fold has only one label, remove the fold by assigning a
        sulcus ID of -1 to all vertices in this fold.

        Case 2:

        If the set of labels in the fold is the same as in one of the protocol's
        sulci, assign the index to the corresponding sulcus label list
        to all vertices in the fold.

        Case 3:

        If the set of labels in the fold is a subset of labels in
        only one of the protocol's sulcus label lists, assign the index
        for that list to all vertices in the fold.
        (Ex: the labels in the fold are [1,2,3] and there is only one
        sulcus label list containing 1, 2, and 3: [1,2,3,4])

        Case 4:

        If the set of labels in the fold is a subset of labels in more than
        one of the protocol's sulcus label lists, find corresponding sulcus
        label pair lists that contain a label pair with two of the fold labels.

        (a) If there is only one such pair, assign the index for that list
            to all vertices in the fold.

            Ex: the labels in the fold are [1,2,3] and there are
            sulcus label lists [1,2,3,4] and [1,2,3,5] corresponding to
            sulcus label pair lists [[1,4],[[2,4],[3,4]] and [[1,2],[2,3],[3,5]];
            the index to the latter list would be chosen because fold labels
            1 and 2 are in [1,2] and 2 and 3 are in [2,3].

        (b) If there is more than one such label list, the fold is unresolved.

        (c) If there is no such label list, the fold is unresolved.

        Case 5:

        If the set of labels in the fold are a superset of the labels for only
        one of the protocol's sulcus label lists (Ex: the labels in the fold are
        [1,2,3,4] and there is only one list which is its subset, [1,2,3]):

        (5.1) Assign the index for that list to just the vertices in the fold
              with labels in its list (vertices with label 1, 2 or 3).

        (5.2) Treat the remaining vertices (vertices with label 4) as a new fold
              and start again from #1.

        Case 6:

        If the labels in the fold are a superset of the labels in more than
        one of the sulcus label lists (Ex: labels in the fold are [1,2,3,4]
        and there are two sulcus label lists [1,2,3] and [2,3,4]), segment:

        (6.1) Find all label boundaries within the fold (vertices with two
              different labels in its neighborhood of connected vertices:
              [1,2],[1,3],[2,3],[2,4],[3,4] in the example below).

              Example
              -------

                |       |
              1 |   2   | 4
                +-------+
              1 |   3   | 4
                |       |

        (6.2) Find label boundary label pairs that are not members
              of the label pairs list ([2,3]).

              The example fold contains two sulci and four labels:
              Sulcus 1: [[1,2], [1,3]]
              Sulcus 2: [[2,4], [3,4]]
              and a non-sulcus label boundary: [[2,3]]
              We wish to segment this fold into two separate sulci.

                |       |
              1 |   2   | 4
                xxxxxxxxx
              1 |   3   | 4
                |       |

        Case 6a:

        Some fold vertices have labels independent of non-protocol
        boundary labels.

        (6.4) Segment groups of connected subfold vertices into separate subfolds.

                |       |
              1 |       | 4
                +       +
              1 |       | 4
                |       |

        (6.5) Assign each vertex with a nonpair label to the closest subfold.

                |   |   |
              1 | 2 | 2 | 4
                + - - - +
              1 | 3 | 3 | 4
                |   |   |

        (6.6) Treat each subfold as a new fold and start again from #1.

        Case 6b:

        All fold vertices have labels in common with non-protocol boundaries,
        and the fold is unresolved.

        Case 7:

        If the fold labels are neither a superset nor a subset of any of the
        sulcus label lists, the fold is unresolved.

    """
    import numpy as np
    from utils.mesh_operations import segment, detect_boundaries, \
        compute_distance

    #---------------------------------------------------------------------------
    # Nested function definitions
    #---------------------------------------------------------------------------
    def find_superset_lists(labels, sulcus_label_lists):
        """
        Find lists within *sulcus_label_lists* that are supersets
        of a list of *labels*.

        Parameters
        ----------
        labels : list of integers, indexed from 1
            some labels
        sulcus_label_lists : list of lists of integers
            each list contains the labels defining a sulcus

        Returns
        -------
        subset_indices : list of integers, indexed from zero
            indices to sulcus_label_lists

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

    def find_subset_pairs(fold_labels, label_pair_lists, selected_indices=[]):
        """
        Find label pairs within *label_pair_lists*
        that are subsets of *fold_labels*.

        Parameters
        ----------
        fold_labels : list of integers, indexed from 1
            labels for all vertices in a fold
        label_pair_lists : list of sublists of subsublists of integers
            each subsublist contains a pair of labels, and the sublist of these
            label pairs represents the label boundaries defining a sulcus
        selected_indices : list of integers, starting from 0
            indices to label_pair_lists

        Returns
        -------
        subset_indices : list of integers, indexed from zero
            indices to label pairs in label_pair_lists

        Example
        -------
        >>> find_subset_pairs([1,2,3],[[[1,2]]])
         [0]
        >>> find_subset_pairs([1,2,3],[[[1,2],[3,4]], [[1,3]], [[1,5],[4,3]]])
         [0, 1]

        """

        if selected_indices==[]:
            selected_indices = range(len(label_pair_lists))

        fold_labels = set(fold_labels)
        subset_indices = []
        for i in selected_indices:
            row = label_pair_lists[i]
            for pair in row:
                if set(pair).issubset(fold_labels):
                    subset_indices.append(i)

        return subset_indices

    def find_subset_lists(labels, sulcus_label_lists):
        """
        Find lists within *sulcus_label_lists* that are subsets of *labels*.

        Parameters
        ----------
        labels : list of integers, indexed from 1
            labels for all vertices
        sulcus_label_lists : list of lists of integers
            each list contains the labels defining a sulcus

        Returns
        -------
        superset_index : list of integers, indexed from zero
            indices to lists within sulcus_label_lists

        Example
        -------
        >>> find_subset_lists([1,2,3,4],[[1,2],[2,3,4]])
         [0, 1]

        """

        labels = set(labels)
        superset_index = []
        for Id, pair_list in enumerate(sulcus_label_lists):
            if set(pair_list).issubset(labels):
                superset_index.append(Id)

        return superset_index

    #---------------------------------------------------------------------------
    # Prepare data structures
    #---------------------------------------------------------------------------
    # Prepare list of sulcus label lists (labels within a sulcus)
    sulcus_label_lists = []
    for row in label_pair_lists:
        sulcus_label_lists.append(list(np.unique(np.asarray(row))))
    if verbose:
        print('\n  sulcus_label_lists:  {0}'.format(sulcus_label_lists))

    # Prepare list of all unique labels in the labeling protocol
    label_list = []
    [label_list.append(x)
     for lst in sulcus_label_lists
     for x in lst
     if x not in label_list]

    # Prepare list of all unique sorted label pairs in the labeling protocol
    protocol_pairs = []
    [protocol_pairs.append(list(set(x)))
     for lst in label_pair_lists
     for x in lst
     if list(set(x)) not in protocol_pairs]

    #---------------------------------------------------------------------------
    # Loop through folds
    #---------------------------------------------------------------------------
    for ifold, fold in enumerate(folds):
        fold_labels = [labels[vertex] for vertex in fold]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels)]
        if verbose:
            print('\n  Fold {0} of {1}\n  labels in fold: {2}'.format(
                ifold + 1, max(fold_IDs), unique_fold_labels))

        # Case 1:
        # If the fold has only one label, remove the fold by assigning a
        # sulcus ID of -1 to all vertices in this fold.
        if len(unique_fold_labels) <= 1:
            print "  Case 1: fold has only one label  -- UNASSIGNED"
            continue  # sulcus_IDs already initialized with -1 values

        # Case 2:
        # If the set of labels in the fold is the same as in one of the
        # protocol's sulci, assign the index to the corresponding sulcus
        # label list to all vertices in the fold.
        elif unique_fold_labels in sulcus_label_lists:
            print "  Case 2: fold labels match one sulcus label list"
            # Add 1 because the sulcus ID begins from 1
            Index = sulcus_label_lists.index(unique_fold_labels) + 1
            if verbose:
                print('    {0}'.format(sulcus_label_lists[Index-1]))
            for vertex in fold:
                sulcus_IDs[vertex] = Index

        # Cases 3-6:
        else:
            # Find lists within the sulcus label lists that are supersets
            # of the list of labels in this fold.
            superset_list_indices = find_superset_lists(unique_fold_labels,
                                                        sulcus_label_lists)
            # Cases 3-4:
            if len(superset_list_indices):

                # Case 3:
                # If the set of labels in the fold is a subset of labels in
                # only one of the protocol's sulcus label lists, assign the
                # index for that list to all vertices in the fold.
                if len(superset_list_indices) == 1:
                    print "  Case 3: fold labels a subset of only one sulcus label list"
                    if verbose:
                        print('    {0}'.format(sulcus_label_lists[superset_list_indices[0]]))
                    for vertex in fold:
                        sulcus_IDs[vertex] = superset_list_indices[0] + 1

                # Case 4:
                # If the set of labels in the fold is a subset of labels in
                # more than one of the protocol's sulcus label lists, find
                # corresponding sulcus label pair lists that contain a label
                # pair with two of the fold labels.
                elif len(superset_list_indices) > 1:
                    print "  Case 4: fold labels a subset of multiple sulcus label lists"
                    if verbose:
                        print('Sulcus label superset lists:\n')
                        for index in superset_list_indices:
                            print('    {0}'.format(sulcus_label_lists[index]))

                    subset_pair_indices = find_subset_pairs(unique_fold_labels,
                        label_pair_lists, selected_indices=superset_list_indices)
                    if len(subset_pair_indices):
                        if verbose:
                            print('Sulcus label subset pairs:\n')
                            for index in subset_pair_indices:
                                print('      {0}'.format(sulcus_label_lists[index]))

                        # (a) If there is only one such pair, assign the index
                        #     for that list to all vertices in the fold.
                        if len(subset_pair_indices) == 1:
                            print "  Case 4a: only one sulcus label pair list " \
                                  "contains a label pair with two fold labels"
                            if verbose:
                                print('    {0}'.format(
                                    sulcus_label_lists[subset_pair_indices[0]]))
                            for vertex in fold:
                                sulcus_IDs[vertex] = subset_pair_indices[0] + 1
                        # (b) If there is more than one such label list,
                        #     the fold is unresolved.
                        else:
                            print "  Case 4b: multiple sulcus label pair " \
                                  "lists contain a label pair with two of " \
                                  "the fold labels -- UNASSIGNED"
                    else:
                        # (c) If there is no such label list, the fold is unresolved.
                        print "  Case 4c: no sulcus label pair lists " \
                              "contain a label pair with two of the fold " \
                              "labels -- UNASSIGNED"

            # Cases 5-7:
            else:
                # Find lists within sulcus label lists that are subsets
                # of labels in this fold.
                subset_list_indices = find_subset_lists(unique_fold_labels,
                                                        sulcus_label_lists)
                # Cases 5-6:
                if len(subset_list_indices):

                    # Case 5:
                    # If the set of labels in the fold are a superset of the
                    # labels for only one of the protocol's sulcus label lists:
                    if len(subset_list_indices) == 1:
                        print "  Case 5: fold labels a superset of only one " + \
                              "sulcus label list"
                        if verbose:
                            print('    {0}'.format(
                                sulcus_label_lists[subset_list_indices[0]]))

                        unassigned = [[]]
                        for vertex in fold:
                            # (5.1) Assign the index for that list to just the
                            #       vertices in the fold with labels in its list.
                            if labels[vertex] in sulcus_label_lists[subset_list_indices[0]]:
                                sulcus_IDs[vertex] = subset_list_indices[0] + 1
                            # (5.2) Treat the remaining vertices as a new fold
                            #       and start again from #1.
                            else:
                                unassigned[0].append(vertex)
                        if len(unassigned):
                            sulcus_IDs = identify(labels, unassigned,
                                label_pair_lists, sulcus_IDs, neighbor_lists,
                                points)

                    # Case 6:
                    # If the labels in the fold are a superset of the labels
                    # in more than one of the sulcus label lists, segment:
                    else:
                        print "  Case 6: fold labels have multiple subset " + \
                              "sulcus label lists"
                        if verbose:
                            for index in subset_list_indices:
                                print('    {0}'.format(sulcus_label_lists[index]))

                        # (6.1) Find all label boundaries within the fold
                        foo, boundary_label_pairs = detect_boundaries(labels,
                            fold, neighbor_lists)

                        # (6.2) Find label boundary label pairs in the fold
                        #       that are not members of the protocol's label
                        #       pairs list and return the unique list of
                        #       'nonpair labels'.
                        nonpair_labels = []
                        for vertex_labels in boundary_label_pairs:
                            if np.unique(vertex_labels).tolist() not in protocol_pairs:
                                [nonpair_labels.append(x) for x in vertex_labels
                                 if x not in nonpair_labels]

                        # (6.3) Find fold vertices whose labels are nonpair labels
                        #       and consider the remaining to be 'subfold' vertices.
                        nonpair_vertices = [i for i, x in enumerate(fold_labels)
                                            if x in nonpair_labels]
                        subfold_vertices = [i for i, x in enumerate(fold_labels)
                                            if x not in nonpair_labels]

                        if len(subfold_vertices):
                            print("  Case 6a: fold vertices have labels " + \
                                  "independent of non-protocol boundary labels")
                            subfold_points = [points[x] for x in subfold_vertices]

                            # (6.4) Segment groups of connected subfold vertices,
                            #       if any, into separate subfolds.
                            subfolds, n_subfolds, foo = segment(subfold_vertices,
                                neighbor_lists, 50)

                            # (6.5) Assign each vertex with a nonpair label
                            #       to the closest subfold.
                            # Create a list of lists of subfolds
                            subfold_lists = []
                            for i in range(n_subfolds):
                                subfold_list = [x for x in subfold_vertices
                                                if subfolds[x] == i + 1]
                                subfold_lists.append(subfold_list)
                            # For each nonpair vertex find the closest subfold vertex
                            for nonpair_vertex in nonpair_vertices:
                                foo, index = compute_distance(points[nonpair_vertex],
                                                              subfold_points)
                                # Assign the subfold index of the closest vertex
                                if subfolds[index]:
                                    subfold_lists[subfolds[index] - 1].append(index)

                            # (6.6) Treat each subfold as a new fold and
                            #       start again from #1.
                            sulcus_IDs = identify(labels, subfold_lists,
                                label_pair_lists, sulcus_IDs, neighbor_lists, points)

                        else:
                            print "  Case 6b: all fold vertices have " + \
                                  "labels in common with non-protocol " + \
                                  "boundaries -- UNASSIGNED"

                # Case 7: If the fold labels are neither a superset nor subset
                #         of any sulcus label list, the fold is unresolved.
                else:
                    print "  Case 7: fold labels neither a superset nor " + \
                          "subset of any sulcus label list -- UNASSIGNED"

    return sulcus_IDs


if __name__ == "__main__":

    import sys
    from utils.io_vtk import load_scalar, rewrite_scalars
    from info.sulcus_boundaries import sulcus_boundaries
    from utils.mesh_operations import find_neighbors

    # Load labels and folds (the second surface has to be inflated).
    points, faces, labels, n_vertices = load_scalar(sys.argv[1], return_arrays=0)
    points, faces, fold_IDs, n_vertices = load_scalar(sys.argv[2], return_arrays=0)

    # Calculate neighbor lists for all vertices
    neighbor_lists = find_neighbors(faces, len(points))

    # Prepare list of all unique sorted label pairs in the labeling protocol
    label_pair_lists = sulcus_boundaries()

    # Create a list of lists of folds
    unique_fold_IDs = set(fold_IDs)
    folds = []
    for id in unique_fold_IDs:
        if id > 0:
            fold = [i for i,x in enumerate(fold_IDs) if x == id]
            folds.append(fold)

    # Assign a sulcus ID to each fold vertex.
    # Initialize all fold vertices as -1.
    sulcus_IDs = [-1 for i in xrange(n_vertices)]
    # Since we do not touch gyral vertices and vertices whose labels
    # are not in label list, or vertices having only one label,
    # their sulcus IDs will remain -1.
    sulcus_IDs = identify(labels, folds, label_pair_lists,
        sulcus_IDs, neighbor_lists, points)

    # Finally, write points, faces and sulcus_IDs to a new vtk file
    rewrite_scalars(sys.argv[1], sys.argv[3], sulcus_IDs, filter_scalars=[])