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
             neighbor_lists, points, sulcus_names=''):
    """
    Identify sulcus folds in a brain surface according to a labeling protocol
    that includes a list of label pairs defining each sulcus.

    Parameters
    ----------
    labels : list of integers
        labels for all vertices
    folds : list of lists of integers
        each list contains indices to vertices of a fold
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus
    sulcus_IDs : list of integers
        sulcus IDs assigned or to be assigned to all vertices
    points : list of lists of three floats
        coordinates of all vertices on the surface
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    sulcus_IDs : list of integers
        sulcus IDs for all vertices, with -1s for non-sulcus vertices

    Definitions
    -----------

    A fold is a group of connected vertices that are deeper than a depth threshold.

    The ''label list'' is a list of all labels that help to define all sulci
    in the DKT protocol.

    The ''label pair list'' contains all pairs of labels from the protocol.

    A ''sulcus label list'' is the list of labels used to define a single sulcus.
    A sulcus label list is a subset of the label list.

    A ''sulcus label pair list'' contains pairs of labels, where each pair
    defines a boundary between two labeled regions.
    No two sulcus label pair lists share a label pair.

    A ''sulcus ID'' uniquely identifies a sulcus.
    It is the index to a particular sulcus label list (or sulcus label pair list)
    in its parent list.

    Algorithm
    ---------

    First, remove all fold vertices with labels that are not in the label list
    by assigning -1 to them.

    For each fold (vertices with the same non-zero label):

        Case 1:  Unassigned -- isolated fold (only one label)

        If the fold has only one label, remove the fold by assigning a
        sulcus ID of -1 to all vertices in this fold.

        Case 2:  fold labels a perfect match with one sulcus

        If the set of labels in the fold is the same as in one of the protocol's
        sulci, assign the index to the corresponding sulcus label list
        to all vertices in the fold.

        Case 3:  fold labels a subset of only one sulcus

        If the set of labels in the fold is a subset of labels in
        only one of the protocol's sulcus label lists, assign the index
        for that list to all vertices in the fold.
        (Ex: the labels in the fold are [1,2,3] and there is only one
        sulcus label list containing 1, 2, and 3: [1,2,3,4])

        Case 4:  AMBIGUOUS -- fold labels a subset of more than one sulcus

        If the set of labels in the fold is a subset of labels in more than
        one of the protocol's sulcus label lists, find corresponding sulcus
        label pair lists that share a label pair with the fold label pairs.

        (a) fold label pair in only one sulcus

            If there is only one such pair, assign the index for that list
            to all vertices in the fold.

            Ex: the labels in the fold are [1,2,3] with label pairs [1,2] and
            [1,3], and there are sulcus label lists [1,2,3,4] and [1,2,3,5]
            for sulcus label pair lists [[1,4],[[2,4],[3,4]] and
            [[1,2],[2,3],[3,5]]; the index to the latter list would be chosen
            because the fold shares the label pair [1,2].

        (b) UNRESOLVED -- fold label pair in more than one sulcus

            If there is more than one such label list, the fold is unresolved.

        (c) UNRESOLVED -- no shared label pair

            If there is no such pair, the fold is unresolved.

        Case 5:  PARTIAL ASSIGNMENT -- fold labels a superset of only one sulcus

        If the set of labels in the fold are a superset of the labels for only
        one of the protocol's sulcus label lists (Ex: the labels in the fold are
        [1,2,3,4] and there is only one list which is its subset, [1,2,3]):

        (5.1) Assign the index for that list to just the vertices in the fold
              with labels in its list (vertices with label 1, 2 or 3).

        (5.2) Treat the remaining vertices (vertices with label 4) as a new fold
              and start again from #1.

        Case 6:  AMBIGUOUS -- fold labels a superset of more than one sulcus

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

        Case 7:  AMBIGUOUS -- fold labels neither a subset nor superset of any sulcus

        If the fold labels are neither a superset nor a subset of any of the
        sulcus label lists, the fold is unresolved.

    """
    import numpy as np
    from utils.mesh_operations import segment, detect_boundaries, \
        compute_distance

    #---------------------------------------------------------------------------
    # Nested function definitions
    #---------------------------------------------------------------------------
    def find_superset_subset_lists(labels, label_lists):
        """
        Find *label_lists* that are supersets or subsets of *labels*.

        Parameters
        ----------
        labels : list of integers
            label numbers
        label_lists : list of lists of integers
            each list contains label numbers

        Returns
        -------
        superset_indices : list of integers
            indices to label_lists that are a superset of labels
        subset_indices : list of integers
            indices to label_lists that are a subset of labels

        Example
        -------
        >>> find_superset_subset_lists([1,2],[[1,2],[3,4]])
        >>> [0]
        >>> find_superset_subset_lists([1,2],[[1,2,3],[1,2,5]])
        >>> [0, 1]
        >>> find_superset_subset_lists([1,2],[[2,3],[1,2,5]])
        >>> [1]

        """

        labels = set(labels)
        superset_indices = []
        subset_indices = []
        for Id, label_list in enumerate(label_lists):
            if labels.issubset(set(label_list)):
                superset_indices.append(Id)
            if set(label_list).issubset(labels):
                subset_indices.append(Id)

        return superset_indices, subset_indices

    def find_pairs(label_pairs, label_pair_lists):
        """
        Find sublists within *label_pair_lists* that contain
        at least one label pair within *label_pairs*.

        Parameters
        ----------
        label_pairs : list of sublists of (pairs of) integers
            each sublist contains a pair of labels
        label_pair_lists : list of sublists of subsublists of (pairs of) integers
            this list contains sublists like label_pairs

        Returns
        -------
        indices : list of integers
            indices to the sublists of label_pair_lists that contain
            at least one label pair in *label_pairs*.

        Example
        -------
        >>> find_pairs([[2,1]], [[[0,2],[1,2]]])
         [1]
        >>> find_pairs([[1,2],[3,1]], [[[1,4],[1,2]], [[1,3]], [[1,5],[2,3]]])
         [0,1]

        """

        indices = []
        for label_pair in label_pairs:
            for index, label_pair_list in enumerate(label_pair_lists):
                if list(set(label_pair)) in label_pair_list:
                    if index not in indices:
                        indices.append(index)

        return indices

    #---------------------------------------------------------------------------
    # Prepare data structures
    #---------------------------------------------------------------------------
    # Prepare list of sulcus label lists (labels within a sulcus)
    sulcus_label_lists = []
    for row in label_pair_lists:
        sulcus_label_lists.append(list(np.unique(np.asarray(row))))

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

        # List the labels in this fold
        fold_labels = [labels[vertex] for vertex in fold]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels)]
        if verbose:
            print("\n  Fold {0} vs. protocol:".format(ifold + 1))

        # Case 1: NO PAIR -- only one label in fold
        if len(unique_fold_labels) <= 1:
            print("  [1] NO PAIR -- only one label in fold: {0}".format(
                unique_fold_labels[0]))
            # Ignore: sulcus_IDs already initialized with -1 values
            continue

        # Case 2: perfect match between sulcus and fold labels
        elif unique_fold_labels in sulcus_label_lists:
            print("  [2] perfect match: {0} ({1})".
                format(unique_fold_labels,
                sulcus_names[sulcus_label_lists.index(unique_fold_labels)]))
            index = sulcus_label_lists.index(unique_fold_labels)
            # Assign the ID of the matching sulcus to all vertices in the fold
            for vertex in fold:
                sulcus_IDs[vertex] = index

        # Cases 3-7:
        else:
            # Find all label boundary pairs within the fold
            foo, fold_pairs, unique_fold_pairs = detect_boundaries(labels,
                fold, neighbor_lists)

            if verbose:
                print("  Label pairs in fold: {0}".format(unique_fold_pairs))

            # Find indices to label pair lists that contain fold label pairs
            list_indices = find_pairs(unique_fold_pairs, label_pair_lists)

            # Case 3: DISJOINT -- no label pair in common between sulci and fold
            if not len(list_indices):
                print("  [3] DISJOINT -- no sulcus label pair")

            # Cases 4-7:
            # There is at least one label pair in common between sulci and fold
            else:

                # Find sulcus label lists that contain or are contained by
                # the fold labels
                superset_indices, subset_indices = find_superset_subset_lists(
                    unique_fold_labels, sulcus_label_lists)

                # Cases 4-5:
                # If at least one sulcus label list contains fold labels
                if len(superset_indices):

                    # Case 4: only one sulcus contains the fold labels
                    if len(superset_indices) == 1:
                        print("  [4] one match -- sulcus contains "
                              "the fold labels: {0} ({1})".
                            format(sulcus_label_lists[superset_indices[0]],
                                   sulcus_names[superset_indices[0]]))
                        # Assign ID of the matching sulcus to all vertices in the fold
                        for vertex in fold:
                            sulcus_IDs[vertex] = superset_indices[0]

                    # Case 5: UNRESOLVED -- more than one sulcus contains the fold labels
                    else:
                        print("  [5] UNRESOLVED -- "
                              "more than one sulcus contains the fold labels")
                        if verbose:
                            print("Sulcus label superset lists:")
                            for index in superset_indices:
                                print("    {0} ({1})".format(
                                    sulcus_label_lists[index],
                                    sulcus_names[index]))

                # Cases 6-7:
                # If at least one sulcus label list within the fold labels
                elif len(subset_indices):

                    # Case 6: PARTIAL ASSIGNMENT -- one sulcus within fold labels
                    if len(subset_indices) == 1:
                        print("  [6] PARTIAL ASSIGNMENT -- "
                            "sulcus within fold labels: {0} ({1})".format(
                            sulcus_label_lists[superset_indices[0]],
                            sulcus_names[superset_indices[0]]))

                        # Assign sulcus index to fold vertices that have
                        # labels in the corresponding sulcus label list
                        unassigned = [[]]
                        for vertex in fold:
                            if labels[vertex] in sulcus_label_lists[subset_indices[0]]:
                                sulcus_IDs[vertex] = subset_indices[0]
                            else:
                                unassigned[0].append(vertex)

                        # Treat remaining vertices as a new fold and start again
                        if len(unassigned):
                            sulcus_IDs = identify(labels, unassigned,
                                label_pair_lists, sulcus_IDs, neighbor_lists,
                                points, sulcus_names)

                    # Case 7: AMBIGUOUS -- more than one sulcus within fold labels
                    else:
                        print("  [7] AMBIGUOUS -- "
                              "more than one sulcus within fold labels")
                        if verbose:
                            for index in subset_indices:
                                print("    {0} ({1})".format(
                                    label_pair_lists[index],
                                    sulcus_names[sulcus_label_lists.index(
                                        sulcus_label_lists[index])]))

                        # Find label boundary pairs in the fold that are not
                        # members of the protocol's label pairs list and
                        # return the unique list of 'nonpair labels'
                        nonpair_labels = []
                        for vertex_labels in fold_pairs:
                            if np.unique(vertex_labels).tolist() not in protocol_pairs:
                                [nonpair_labels.append(x) for x in vertex_labels
                                 if x not in nonpair_labels]

                        # Find fold vertices whose labels are nonpair labels
                        # and consider the remaining to be 'subfold' vertices
                        nonpair_vertices = [i for i, x in enumerate(fold_labels)
                                            if x in nonpair_labels]
                        subfold_vertices = [i for i, x in enumerate(fold_labels)
                                            if x not in nonpair_labels]

                        if len(subfold_vertices):
                            print("  [7a] fold vertices have labels "
                                  "independent of non-protocol boundary labels")
                            subfold_points = [points[x] for x in subfold_vertices]

                            # Segment groups of connected subfold vertices,
                            # if any, into separate subfolds
                            subfolds, n_subfolds, foo = segment(subfold_vertices,
                                neighbor_lists, 50)

                            # Assign each vertex with a nonpair label
                            # to the closest subfold
                            # Create a list of lists of subfolds
                            subfold_lists = []
                            for i in range(n_subfolds):
                                subfold_list = [x for x in subfold_vertices
                                                if subfolds[x] == i]
                                subfold_lists.append(subfold_list)
                            # For each nonpair vertex find the closest subfold vertex
                            for nonpair_vertex in nonpair_vertices:
                                foo, index = compute_distance(points[nonpair_vertex],
                                                              subfold_points)
                                # Assign the subfold index of the closest vertex
                                if subfolds[index]:
                                    subfold_lists[subfolds[index] - 1].append(index)

                            # Treat each subfold as a new fold and start again
                            sulcus_IDs = identify(labels, subfold_lists,
                                label_pair_lists, sulcus_IDs, neighbor_lists,
                                points, sulcus_names)

    return sulcus_IDs


if __name__ == "__main__":

    """
    run extract/identify_sulci.py /Applications/freesurfer/subjects/MMRR-21-1/label/lh.labels.DKT25.manual.vtk /Users/arno/Desktop/output/results/features/_hemi_lh_subject_MMRR-21-1/sulci.vtk /Users/arno/Documents/Projects/Mindboggle/mindboggle/mindboggle/info/sulcus_names.txt test_identify.vtk
    """
    import sys
    from utils.io_vtk import load_scalar, rewrite_scalars
    from utils.io_file import read_columns
    from info.sulcus_boundaries import sulcus_boundaries
    from utils.mesh_operations import find_neighbors

    # Arguments
    labels_file = sys.argv[1]
    folds_file = sys.argv[2]
    sulcus_names_file = sys.argv[3]
    vtk_file = sys.argv[4]

    # Load labels and folds (the second surface has to be inflated).
    points, faces, labels, n_vertices = load_scalar(labels_file, return_arrays=0)
    points, faces, fold_IDs, n_vertices = load_scalar(folds_file, return_arrays=0)
    sulcus_names = read_columns(sulcus_names_file, n_columns=1)[0]

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
        sulcus_IDs, neighbor_lists, points, sulcus_names)

    # Finally, write points, faces and sulcus_IDs to a new vtk file
    rewrite_scalars(labels_file, vtk_file, sulcus_IDs, filter_scalars=sulcus_IDs)
