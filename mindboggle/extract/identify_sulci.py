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

    A fold is a group of connected vertices deeper than a depth threshold.

    A ''label list'' is the list of labels used to define a single sulcus.

    A ''label pair list'' contains pairs of labels, where each pair
    defines a boundary between two labeled regions.
    No two label pair lists share a label pair.

    A ''sulcus ID'' uniquely identifies a sulcus.
    It is the index to a sulcus label list (or sulcus label pair list).

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

    #---------------------------------------------------------------------------
    # Prepare data structures
    #---------------------------------------------------------------------------
    # Prepare list of sulcus label lists (labels within a sulcus)
    label_lists = []
    for row in label_pair_lists:
        label_lists.append(list(np.unique(np.asarray(row))))

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
        unassigned = []
        seed_lists = []

        # List the labels in this fold
        fold_labels = [labels[vertex] for vertex in fold]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels)]

        # Case 0: NO LABEL
        if not len(unique_fold_labels):
            print("  Fold {0}: NO LABEL")
            # Ignore: sulcus_IDs already initialized with -1 values
            continue

        # Case 1: NO PAIR -- only one label in fold
        elif len(unique_fold_labels) == 1:
            print("  Fold {0}: LONE LABEL ({1})".
                  format(ifold, unique_fold_labels[0]))
            # Ignore: sulcus_IDs already initialized with -1 values
            continue

        # Case 2: perfect match between sulcus and fold labels
        elif unique_fold_labels in label_lists:
            print("  Fold {0}: perfect match: {1} ({2})".
                format(ifold,
                sulcus_names[label_lists.index(unique_fold_labels)],
                  ', '.join([str(x) for x in unique_fold_labels])))
            index = label_lists.index(unique_fold_labels)
            # Assign ID of the matching sulcus to all fold vertices
            for vertex in fold:
                sulcus_IDs[vertex] = index

        # Cases 3-7:
        # There is at least one fold label but no perfect match
        else:
            # Find all label boundary pairs within the fold
            indices_fold_pairs, fold_pairs, unique_fold_pairs = detect_boundaries(
                labels, fold, neighbor_lists)
            if verbose:
                print("  Fold labels: {0}".
                      format(', '.join([str(x) for x in unique_fold_labels])))
                #print("  Fold label pairs: {0}".
                #      format(', '.join([str(x) for x in unique_fold_pairs])))

            # Find fold label pairs in the protocol
            unique_fold_pairs = [x for x in unique_fold_pairs
                                 if np.unique(x).tolist() in protocol_pairs]
            if verbose:
                print("  Fold label pairs in protocol: {0}".
                      format(', '.join([str(x) for x in unique_fold_pairs])))

            # Case 3: NO MATCHING PAIR -- no label pair in common
            #                             between sulci and fold
            if not len(unique_fold_pairs):
                print("  Fold {0}: NO MATCHING PAIR".format(ifold))

            # Cases 4-7:
            # There is at least one label pair in common between sulci and fold
            else:

                # Find sulcus label lists that contain or are contained by
                # the fold labels
                superset_indices, subset_indices = find_superset_subset_lists(
                    unique_fold_labels, label_lists)

                # Cases 4-5:
                # If fold labels contained by one or more sulcus label lists
                if len(superset_indices):

                    # Case 4: one match -- fold labels in one sulcus
                    if len(superset_indices) == 1:
                        print("  Fold {0}: one match -- "
                              "fold labels in one sulcus: {1} ({2})".
                            format(ifold, sulcus_names[superset_indices[0]],
                            ', '.join([str(x)
                            for x in label_lists[superset_indices[0]]])))
                        # Assign ID of matching sulcus to all fold vertices
                        for vertex in fold:
                            sulcus_IDs[vertex] = superset_indices[0]

                    # Case 5: UNRESOLVED -- fold labels in more than one sulcus
                    # NOTE: Could check to see which sulcus shares a label pair
                    else:
                        print("  Fold {0}: UNRESOLVED -- "
                              "fold labels in more than one sulcus".
                              format(ifold))
                        if verbose:
                            print("           Sulci containing fold labels:")
                            for index in superset_indices:
                                print("           {0} ({1})".format(
                                    label_lists[index],
                                    sulcus_names[index]))

                # Cases 6-7:
                # If the fold labels contain at least one sulcus label list
                elif len(subset_indices):

                    # Case 6: PARTIAL MATCH -- fold labels contain one sulcus
                    if len(subset_indices) == 1:
                        supind = label_lists[superset_indices[0]]
                        print("  Fold {0}: PARTIAL MATCH -- "
                            "fold labels contain one sulcus: {1} ({2})".format(
                            ifold, sulcus_names[superset_indices[0]],
                            ', '.join([str(x) for x in supind])))

                        # Assign a sulcus ID to fold vertices that have
                        # labels in the corresponding sulcus label list,
                        # and save unassigned vertices
                        unassigned = []
                        for vertex in fold:
                            if labels[vertex] in label_lists[subset_indices[0]]:
                                sulcus_IDs[vertex] = subset_indices[0]
                            else:
                                unassigned.append(vertex)

                    # Case 7: AMBIGUOUS -- fold labels contain more than one sulcus
                    else:
                        print("  Fold {0}: AMBIGUOUS -- "
                              "fold labels contain more than one sulcus".
                              format(ifold))
                        #if verbose:
                        #    for index in subset_indices:
                        #        print("           {0} ({1})".format(
                        #            ', '.join([str(x)
                        #                       for x in label_pair_lists[index]]),
                        #            sulcus_names[label_lists.index(
                        #                label_lists[index])]))

                        # Find label boundary pairs in the fold whose labels
                        # are and are not shared by any other label pairs
                        # in the fold, and store the sulcus IDs for these pairs
                        unique_pairs = []
                        IDs_unique_pairs = []
                        remainder_pairs = []
                        IDs_remainder_pairs = []
                        for pair in unique_fold_pairs:
                            labels_in_pairs = [x for sublst in unique_fold_pairs
                                               for x in sublst]
                            if len([x for x in labels_in_pairs if x in pair]) == 2:
                                unique_pairs.append(pair)
                                IDs_unique_pairs.extend(
                                    [i for i,x in enumerate(label_pair_lists)
                                     if np.sort(pair).tolist() in x])
                            else:
                                remainder_pairs.append(pair)
                                IDs_remainder_pairs.extend([i
                                    for i,x in enumerate(label_pair_lists)
                                    if np.sort(pair).tolist() in x])

                        # Assign a sulcus ID to fold vertices that have unique
                        # label pair labels, and store unassigned vertices
                        if len(unique_pairs):
                            print("           Assign sulcus IDs to vertices "
                                  "with labels in only one fold pair")
                            if verbose:
                                for index, ID_unique_pair \
                                  in enumerate(IDs_unique_pairs):
                                    print("           {0} ({1})".format(
                                        unique_pairs[index],
                                        sulcus_names[ID_unique_pair]))
                            unassigned = []
                            for vertex in fold:
                                index = [i for i,x in enumerate(unique_pairs)
                                         if labels[vertex] in x]
                                if len(index):
                                    sulcus_IDs[vertex] = IDs_unique_pairs[index[0]]
                                else:
                                    unassigned.append(vertex)
                        else:
                            unassigned = fold[:]

                        # Segment remaining vertices into connected sets
                        # of vertices using label boundaries as seeds,
                        # and assign a sulcus ID to each resulting segment
                        if len(unassigned) and len(remainder_pairs):
                            print("           Assign sulcus IDs to vertices "
                                  "with labels in only one fold pair")
                            if verbose:
                                for index, ID_remainder_pair \
                                  in enumerate(IDs_remainder_pairs):
                                    print("           {0} ({1})".format(
                                        remainder_pairs[index],
                                        sulcus_names[ID_remainder_pair]))
                            unassigned = []
                            for vertex in fold:
                                index = [i for i,x in enumerate(unique_pairs)
                                         if labels[vertex] in x]
                                if len(index):
                                    sulcus_IDs[vertex] = IDs_remainder_pairs[index[0]]
                                else:
                                    unassigned.append(vertex)


                            if len(unassigned):
                                seed_lists = [[] for x in remainder_pairs]
                                for i, remainder_pair in enumerate(remainder_pairs):
                            seed_lists[i] = [ for i,x
                                               in enumerate(indices_fold_pairs)
                                               if fold_pairs[i]]]

                        else:
                            print("           Every fold pair has a duplicate label")
                            print("           Segment into separate label-pair regions")


            # Segment fold into separate label-pair regions
            subfolds, n_subfolds, foo = segment_multiple_seeds(seeds_lists,
                neighbor_lists, 50)


        # If there are remaining unassigned vertices
        if len(unassigned):
            print("           Segment {0} unassigned vertices "
                  "into separate subfolds and start again...".
                format(len(unassigned)))

            # Segment groups of connected unassigned vertices into separate subfolds
            subfolds, n_subfolds = segment(unassigned, seed_lists, neighbor_lists, 50)

            # Create a list of lists of subfolds
            subfold_lists = []
            for i in range(n_subfolds):
                subfold_list = [x for x in unassigned
                                if subfolds[x] == i]
                subfold_lists.append(subfold_list)

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
    fid = open(sulcus_names_file, 'r')
    sulcus_names = fid.readlines()
    sulcus_names = [x.strip('\n') for x in sulcus_names]

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
