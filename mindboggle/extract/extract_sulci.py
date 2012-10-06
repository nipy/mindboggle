#!/usr/bin/python
"""
Function to extract sulcus folds from a labeled triangular surface mesh,
using a sulcus labeling protocol.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import numpy as np
#from time import time
#from mindboggle.measure.measure_functions import compute_percentile
#from mindboggle.utils.mesh_operations import find_neighbors, detect_boundaries,\
#    segment, fill_holes
#from mindboggle.utils.io_vtk import load_scalar

#==============
# Extract folds
#==============
def extract_sulci(label_pair_lists, label_file, depth_file,
                  fraction_folds, min_sulcus_size, do_fill_holes=False):
    """
    Use depth and a sulcus labeling protocol to extract sulcus folds
    from an anatomically labeled triangular surface mesh and fill holes
    resulting from shallower areas within a fold.

    Parameters
    ----------
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus
    label_file : str
        surface mesh file in VTK format with faces and label scalar values
    depth_file : str
        surface mesh file in VTK format with faces and depth scalar values
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices
    fraction_folds : float
        fraction of surface mesh considered folds
    min_sulcus_size : int
        minimum sulcus size (number of vertices)
    do_fill_holes : Boolean
        segment and fill holes?

    Returns
    -------
    sulcus_IDs : array of integers
        sulcus IDs for all vertices, with -1s for non-sulcus vertices
    n_sulci :  int
        number of sulcus folds

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

    First, extract folds by a depth threshold, and remove all non-fold
    vertices by assigning -1 to them.

    Construct seed lists of label boundary vertices, segment into sets
    of vertices connected to label boundary seeds that have the same labels,
    and assign a sulcus ID to each segment.

    Example
    -------
    >>> import os
    >>> from mindboggle.extract.extract_sulci import extract_sulci
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> label_pair_lists = sulcus_boundaries()
    >>> fraction_folds = 0.10  # low to speed up
    >>> min_sulcus_size = 50
    >>> sulcus_IDs, n_sulci = extract_sulci(label_pair_lists, label_file,
    >>>     depth_file, fraction_folds, min_sulcus_size, do_fill_holes=False)
    >>> # Write results to vtk file and view with mayavi2:
    >>> from mindboggle.utils.io_vtk import rewrite_scalars, load_scalar
    >>> rewrite_scalars(label_file, 'test_extract_sulci.vtk',
    >>>                 sulcus_IDs, sulcus_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_sulci.vtk &')
    >>> # Write and view manual labels restricted to sulci:
    >>> points, faces, labels, n_vertices = load_scalar(label_file, True)
    >>> rewrite_scalars(label_file, 'test_extract_sulci_labels.vtk',
    >>>                 labels, sulcus_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_sulci_labels.vtk &')

    """
    import numpy as np
    from time import time
    from mindboggle.measure.measure_functions import compute_percentile
    from mindboggle.utils.mesh_operations import find_neighbors, detect_boundaries,\
        segment, fill_holes
    from mindboggle.utils.io_vtk import load_scalar

    #---------------------------------------------------------------------------
    # Load depth and label values from VTK files; find neighboring vertices
    #---------------------------------------------------------------------------
    points, faces, depths, n_vertices = load_scalar(depth_file, return_arrays=True)
    points, faces, labels, n_vertices = load_scalar(label_file, return_arrays=True)
    neighbor_lists = find_neighbors(faces, len(points))

    # Compute the minimum depth threshold for defining folds by determining the
    # percentile of depth values for the fraction of vertices that are not folds.
    # For example, if we consider the shallowest one-third of vertices not to be
    # folds, we compute the depth percentile, and two-thirds of vertices would
    # have at least this depth value and would be considered folds.
    min_depth = compute_percentile(np.sort(depths),
                                   1 - fraction_folds, key=lambda x:x)
    #---------------------------------------------------------------------------
    # Segment sulcus folds of a surface mesh
    #---------------------------------------------------------------------------
    print("  Segment surface mesh into sulcus folds...")
    t0 = time()

    # Array of sulcus IDs for fold vertices, initialized as -1.
    # Since we do not touch gyral vertices and vertices whose labels
    # are not in the label list, or vertices having only one label,
    # their sulcus IDs will remain -1.
    sulcus_IDs = -1 * np.ones(len(neighbor_lists))

    # Extract deep vertices
    print("    Extract vertices deeper than {0:.2f}".format(min_depth))
    indices_folds = np.where(depths > min_depth)[0]
    if len(indices_folds):

        # Find all label boundary pairs within the fold
        indices_boundaries, label_pairs, unique_label_pairs = detect_boundaries(
            indices_folds, labels, neighbor_lists)

        # Construct seed lists containing indices
        # to the vertices of each label boundary
        seed_lists = []
        for label_pair_list in label_pair_lists:
            seed_lists.append([x for i,x in enumerate(indices_boundaries)
                if np.sort(label_pairs[i]).tolist() in label_pair_list])

        # Segment into sets of vertices connected to label boundary seeds
        print("    Segment into separate label-pair regions...")
        t1 = time()
        sulcus_IDs, n_sulci = segment(indices_folds, neighbor_lists,
            seed_lists, min_sulcus_size, spread_same_labels=True,
            labels=labels, label_pair_lists=label_pair_lists)
        print("    ...Segmented {0} sulcus folds in {1:.2f} seconds".
              format(n_sulci, time() - t1))

    else:
        n_sulci = 0

    #---------------------------------------------------------------------------
    # Fill holes in folds
    #---------------------------------------------------------------------------
    if n_sulci > 0 and do_fill_holes:

        # Find fold vertices that have not yet been segmented
        # (because they weren't sufficiently deep)
        t2 = time()
        vertices_to_segment = [i for i,x in enumerate(sulcus_IDs) if x==-1]

        # Segment holes in the folds
        print("    Segment holes...")
        holes, n_holes = segment(vertices_to_segment, neighbor_lists,
            seed_lists=[], min_region_size=1,
            spread_same_labels=False, labels=[], label_pair_lists=[])

        # If there are any holes
        if n_holes > 0:

            # Ignore the largest hole (the background) and renumber holes
            max_hole_size = 0
            max_hole_index = 0
            for ihole in range(n_holes):
                I = np.where(holes == ihole)
                if len(I) > max_hole_size:
                    max_hole_size = len(I)
                    max_hole_index = ihole
            holes[holes == max_hole_index] = -1
            if max_hole_index < n_holes:
                holes[holes > max_hole_index] -= 1
            n_holes -= 1
            print('    ...{0} holes segmented ({1:.2f} seconds)'.
                  format(n_holes, time() - t2))

            # Fill holes
            t3 = time()
            sulcus_IDs = fill_holes(sulcus_IDs, holes, n_holes, neighbor_lists)
            print("    Filled holes ({0:.2f} seconds)".format(time() - t3))

    print("  ...Extracted and filled {0} sulcus folds in {1:.2f} seconds".
          format(n_sulci, time() - t0))

    return sulcus_IDs, n_sulci
