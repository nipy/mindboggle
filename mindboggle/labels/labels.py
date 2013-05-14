#!/usr/bin/env python
"""
Functions for labels.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#-----------------------------------------------------------------------------
# Extract borders between labels
#-----------------------------------------------------------------------------
def extract_borders(indices, labels, neighbor_lists,
                    ignore_values=[], return_label_pairs=False):
    """
    Detect the label boundaries in a collection of vertices such as a region.

    Label boundaries are the set of all vertices
    whose neighbors do not share the same label.

    Parameters
    ----------
    indices : list of integers
        indices to (a subset of) vertices
    labels : numpy array of integers
        label numbers for all vertices (-1 for unlabeled vertices)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    ignore_values : list of integers
        integers to ignore (e.g., background)

    Returns
    -------
    boundary_indices : list of integers
        indices to label boundary vertices
    boundary_label_pairs : list of lists of sorted pairs of integers
        sorted label pairs
    unique_boundary_label_pairs : list of sorted pairs of integers
        unique, sorted label pairs

    Examples
    --------
    >>> # Small example:
    >>> from mindboggle.labels.labels import extract_borders
    >>> from mindboggle.utils.plots import plot_vtk
    >>> indices = [0,1,2,4,5,8,9]
    >>> labels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, -1, -1]
    >>> neighbor_lists = [[1,2,3], [1,2], [2,3], [2], [4,7], [3,2,3]]
    >>> extract_borders(indices, labels, neighbor_lists, [], True)
        ([1, 2, 4, 5],
         [[20, 30], [30, 40], [50, 80], [30, 40]],
         [[20, 30], [30, 40], [50, 80]])
    >>> # Real example -- extract sulcus label boundaries:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.labels import extract_borders
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(labels_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> indices_boundaries, label_pairs, foo = extract_borders(range(npoints),
    >>>     labels, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> IDs = -1 * np.ones(npoints)
    >>> IDs[indices_boundaries] = 1
    >>> rewrite_scalars(labels_file, 'extract_borders.vtk',
    >>>                 IDs, 'boundaries', IDs)
    >>> plot_vtk('extract_borders.vtk')

    """
    import numpy as np

    # Make sure arguments are numpy arrays:
    if not isinstance(labels, np.ndarray):
        labels = np.array(labels)

    # Construct an array of labels corresponding to the neighbor lists:
    L = np.array([list(set(labels[lst])) for lst in neighbor_lists])

    # Find indices to sets of two labels:
    boundary_indices = [indices[y] for y in
                        [i for i,x in enumerate(L[indices])
                         if len(set(x)) == 2]]

    if return_label_pairs:
        boundary_label_pairs = [np.sort(L[indices[j]]).tolist() for j in
                                [i for i,x in enumerate(L[indices])
                                 if len(set(x)) == 2]]
    else:
        boundary_label_pairs = []

    if ignore_values:
        Ikeep = [i for i,x in enumerate(boundary_label_pairs)
                 if not len(frozenset(x).intersection(ignore_values))]
        boundary_indices = [x for i,x in enumerate(boundary_indices)
                            if i in Ikeep]
        if return_label_pairs:
            boundary_label_pairs = [x for i,x in enumerate(boundary_label_pairs)
                                    if i in Ikeep]

    if return_label_pairs:
        unique_boundary_label_pairs = []
        for pair in boundary_label_pairs:
            if pair not in unique_boundary_label_pairs:
                unique_boundary_label_pairs.append(pair)
    else:
        unique_boundary_label_pairs = []

    return boundary_indices, boundary_label_pairs, unique_boundary_label_pairs

#-----------------------------------------------------------------------------
# Extract border values from a second surface
#-----------------------------------------------------------------------------
def extract_borders_2nd_surface(labels_file, mask_file='', values_file=''):
    """
    Extract borders (between labels) on a surface.
    Options: Mask out values; extract border values on a second surface.

    Parameters
    ----------
    labels_file : string
        file name for surface mesh with labels
    mask_file : string
        file name for surface mesh with mask (>-1) values
    values_file : string
        file name for surface mesh with values to extract along borders

    Returns
    -------
    border_file : string
        file name for surface mesh with label borders (-1 background values)
    border_values : numpy array
        values for all vertices (-1 for vertices not along label borders)

    Examples
    --------
    >>> # Extract depth values along label boundaries in sulci (mask):
    >>> import os
    >>> from mindboggle.labels.labels import extract_borders_2nd_surface
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> mask_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> values_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> #
    >>> border_file, border_values = extract_borders_2nd_surface(labels_file, mask_file, values_file)
    >>> #
    >>> plot_vtk(border_file)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.labels.labels import extract_borders

    # Load labeled surface file
    faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(labels_file,
                                                                    return_first=True, return_array=True)

    # Detect boundaries
    neighbor_lists = find_neighbors(faces, npoints)
    indices_boundaries, label_pairs, foo = extract_borders(range(npoints),
                                                           labels, neighbor_lists)

    # Filter values with label boundaries
    border_values = -1 * np.ones(npoints)
    if values_file:
        values, name = read_scalars(values_file, return_first=True, return_array=True)
        border_values[indices_boundaries] = values[indices_boundaries]
    else:
        border_values[indices_boundaries] = 1

    # Mask values (for mask >-1)
    if mask_file:
        mask_values, name = read_scalars(mask_file)
    else:
        mask_values = []

    # Write out label boundary vtk file
    border_file = os.path.join(os.getcwd(), 'borders_' + os.path.basename(labels_file))
    rewrite_scalars(labels_file, border_file, border_values, \
                    'label_borders_in_mask', mask_values)

    return border_file, border_values

#-----------------------------------------------------------------------------
# Majority vote on multiple labels
#-----------------------------------------------------------------------------
def vote_labels(label_lists):
    """
    For each vertex, vote on the majority label.

    Parameters
    ----------
    label_lists : list of lists of integers
        vertex labels assigned by each atlas

    Returns
    -------
    labels_max : list of integers
        majority labels for vertices
    label_counts : list of integers
        number of different labels for vertices
    label_votes : list of integers
        number of votes for the majority labels

    Examples
    --------
    >>> from collections import Counter
    >>> X = [1,1,2,3,4,2,1,2,1,2,1,2]
    >>> Votes = Counter(X)
    >>> Votes
    Counter({1: 5, 2: 5, 3: 1, 4: 1})
    >>> Votes.most_common(1)
    [(1, 5)]
    >>> Votes.most_common(2)
    [(1, 5), (2, 5)]
    >>> len(Votes)
    4

    """
    from collections import Counter

    print("Begin voting...")
    n_atlases = len(label_lists)  # number of atlases used to label subject
    npoints = len(label_lists[0])
    labels_max = [-1 for i in xrange(npoints)]
    label_counts = [1 for i in xrange(npoints)]
    label_votes = [n_atlases for i in xrange(npoints)]

    consensus_vertices = []
    for vertex in xrange(npoints):
        votes = Counter([label_lists[i][vertex] for i in xrange(n_atlases)])

        labels_max[vertex] = votes.most_common(1)[0][0]
        label_votes[vertex] = votes.most_common(1)[0][1]
        label_counts[vertex] = len(votes)
        if len(votes) == n_atlases:
            consensus_vertices.append(vertex)

    print("Voting done.")

    return labels_max, label_votes, label_counts, consensus_vertices

def majority_vote_label(surface_file, annot_files):
    """
    Load a VTK surface and corresponding FreeSurfer annot files.
    Write majority vote labels, and label counts and votes as VTK files.

    Parameters
    ----------
    surface_file : string
        name of VTK surface file
    annot_files : list of strings
        names of FreeSurfer annot files

    Returns
    -------
    labels_max : list of integers
        majority labels for vertices
    label_counts : list of integers
        number of different labels for vertices
    label_votes : list of integers
        number of votes for the majority labels
    consensus_vertices : list of integers
        indicating which are consensus labels
    maxlabel_file : string
        name of VTK file containing majority vote labels
    labelcounts_file : string
        name of VTK file containing number of different label counts
    labelvotes_file : string
        name of VTK file containing number of votes per majority label

    """
    from os import path, getcwd
    import nibabel as nb
    import pyvtk
    from mindboggle.labels.labels import vote_labels
    from mindboggle.utils.io_file import string_vs_list_check

    # Load multiple label sets
    print("Load annotation files...")
    label_lists = []
    for annot_file in annot_files:
        labels, colortable, names = nb.freesurfer.read_annot(annot_file)
        label_lists.append(labels)
    print("Annotations loaded.")

    # Vote on labels for each vertex
    labels_max, label_votes, label_counts, \
    consensus_vertices = vote_labels(label_lists)

    # Check type to make sure the filename is a string
    # (if a list, return the first element)
    surface_file = string_vs_list_check(surface_file)

    # Save files
    VTKReader = pyvtk.VtkData(surface_file)
    Vertices =  VTKReader.structure.points
    Faces =     VTKReader.structure.polygons

    output_stem = path.join(getcwd(), path.basename(surface_file.strip('.vtk')))
    maxlabel_file = output_stem + '.labels.max.vtk'
    labelcounts_file = output_stem + '.labelcounts.vtk'
    labelvotes_file = output_stem + '.labelvotes.vtk'

    pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(labels_max,\
                                  name='Max (majority labels)'))).\
          tofile(maxlabel_file, 'ascii')

    pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(label_counts,\
                                  name='Counts (number of different labels)'))).\
          tofile(labelcounts_file, 'ascii')

    pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(label_votes,\
                                  name='Votes (number of votes for majority labels)'))).\
          tofile(labelvotes_file, 'ascii')

    return labels_max, label_counts, label_votes, consensus_vertices, \
           maxlabel_file, labelcounts_file, labelvotes_file
