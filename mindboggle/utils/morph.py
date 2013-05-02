#!/usr/bin/env python
"""
Morphological image operations on surface mesh vertices.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#-----------------------------------------------------------------------------
# Dilate
#-----------------------------------------------------------------------------
def dilate(indices, nedges, neighbor_lists):
    """
    Dilate region on a surface mesh.

    Parameters
    ----------
    indices : list of integers
        indices of vertices to dilate
    nedges : integer
        number of edges to dilate across
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    dilated_indices : list of integers
        indices of dilated vertices

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.morph import dilate
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> nedges = 3
    >>> # Select a single fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11 #11
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> #
    >>> dilated_indices = dilate(indices, nedges, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> IDs = -1 * np.ones(len(folds))
    >>> IDs[dilated_indices] = 2
    >>> IDs[indices] = 1
    >>> rewrite_scalars(depth_file, 'dilate.vtk', IDs, 'dilated_fold', IDs)
    >>> plot_vtk('dilate.vtk')

    """
    from mindboggle.utils.mesh import find_neighborhood

    N = find_neighborhood(neighbor_lists, indices, nedges)

    dilated_indices = indices[:]
    dilated_indices.extend(N)

    return dilated_indices

#-----------------------------------------------------------------------------
# Erode
#-----------------------------------------------------------------------------
def erode(indices, nedges, neighbor_lists):
    """
    Erode region on a surface mesh.

    Parameters
    ----------
    indices : list of integers
        indices of vertices to erode
    nedges : integer
        number of edges to erode across
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    eroded_indices : list of integers
        indices of eroded vertices

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.morph import erode
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> nedges = 3
    >>> # Select a single fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11 #11
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> #
    >>> eroded_indices = erode(indices, nedges, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> IDs = -1 * np.ones(len(folds))
    >>> IDs[indices] = 1
    >>> IDs[eroded_indices] = 2
    >>> rewrite_scalars(depth_file, 'erode.vtk', IDs, 'eroded_fold', IDs)
    >>> plot_vtk('erode.vtk')

    """
    from mindboggle.utils.mesh import find_neighborhood

    N1 = find_neighborhood(neighbor_lists, indices, nedges=1)
    N2 = find_neighborhood(neighbor_lists, N1, nedges)

    eroded_indices = list(frozenset(indices).difference(N2))

    return eroded_indices

#-----------------------------------------------------------------------------
# Fill holes
#-----------------------------------------------------------------------------
def fill_holes(regions, neighbor_lists, values=[], exclude_range=[]):
    """
    Fill holes in regions on a surface mesh by using region boundaries.

    NOTE: assumes one set of connected vertices per region

    Steps ::

        1. Segment region vertex neighbors into connected vertices (region boundaries).
        2. Remove the largest region boundary, presumably the
           outer contour of the region, leaving smaller boundaries,
           presumably the contours of holes within the region.
        3. Call label_holes() to fill holes with surrounding region numbers.

    Parameters
    ----------
    regions : numpy array of integers
        region numbers for all vertices (default -1)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    values : list of integers
        values for vertices, for use in determining which holes to remove
    exclude_range : list of two floats
        hole is not filled if it contains values within this range
        (prevents cases where surface connected by folds mistaken for holes)

    Returns
    -------
    regions : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors, remove_faces
    >>> from mindboggle.utils.morph import fill_holes
    >>> from mindboggle.labels.segment import segment
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> #
    >>> # Select one fold
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> n_fold = np.unique(folds)[1]
    >>> folds[folds != n_fold] = -1
    >>> #
    >>> # Make two holes in fold (values of -1 and excluded values)
    >>> # Hole 1:
    >>> # Find a vertex whose removal (with its neighbors) would create a hole
    >>> I = np.where(folds==n_fold)[0]
    >>> for index1 in I:
    >>>     N1 = neighbor_lists[index1]
    >>>     stop = True
    >>>     for n in N1:
    >>>         if any(folds[neighbor_lists[n]] == -1):
    >>>             stop = False
    >>>             break
    >>>         else:
    >>>             for f in neighbor_lists[n]:
    >>>                 if any(folds[neighbor_lists[f]] == -1):
    >>>                     stop = False
    >>>                     break
    >>>     if stop:
    >>>         break
    >>> folds[index1] = -1
    >>> folds[N1] = -1
    >>> # Hole 2:
    >>> I = np.where(folds==n_fold)[0]
    >>> for index2 in I:
    >>>     N2 = neighbor_lists[index2]
    >>>     stop = True
    >>>     for n in N2:
    >>>         if any(folds[neighbor_lists[n]] == -1):
    >>>             stop = False
    >>>             break
    >>>         else:
    >>>             for f in neighbor_lists[n]:
    >>>                 if any(folds[neighbor_lists[f]] == -1):
    >>>                     stop = False
    >>>                     break
    >>>     if stop:
    >>>         break
    >>> folds[index2] = -1
    >>> folds[N2] = -1
    >>> values = np.zeros(len(folds))
    >>> values[index2] = 100
    >>> values[N2] = 200
    >>> #
    >>> # Write holes to vtk file and view:
    >>> holes = folds.copy()
    >>> holes[index1] = 10
    >>> holes[N1] = 20
    >>> holes[index2] = 30
    >>> holes[N2] = 40
    >>> indices = [i for i,x in enumerate(holes) if x > -1]
    >>> write_vtk('test_holes.vtk', points, indices, lines,
    >>>           remove_faces(faces, indices), [holes.tolist()], ['holes'])
    >>> from mindboggle.utils.plots import plot_vtk
    >>> plot_vtk('test_holes.vtk')
    >>> #
    >>> # Fill Hole 1 but not Hole 2:
    >>> # (because values has an excluded value in the hole)
    >>> regions = np.copy(folds)
    >>> regions = fill_holes(regions, neighbor_lists, values, [99,101])
    >>> #
    >>> # Write results to vtk file and view:
    >>> indices = [i for i,x in enumerate(regions) if x > -1]
    >>> write_vtk('test_fill_holes.vtk', points, indices, lines,
    >>>           remove_faces(faces, indices), regions.tolist(), 'regions')
    >>> from mindboggle.utils.plots import plot_vtk
    >>> plot_vtk('test_fill_holes.vtk')

    """
    import numpy as np
    from mindboggle.utils.mesh import label_holes
    from mindboggle.labels.segment import segment

    # Make sure argument is a numpy array
    if not isinstance(regions, np.ndarray):
        regions = np.array(regions)

    def label_holes(holes, regions, neighbor_lists):
        """
        Fill holes in regions on a surface mesh.

        Parameters
        ----------
        holes : list or array of integers
            hole numbers for all vertices (default -1)
        regions : numpy array of integers
            region numbers for all vertices (default -1)
        neighbor_lists : list of lists of integers
            each list contains indices to neighboring vertices for each vertex

        Returns
        -------
        regions : numpy array of integers
            region numbers for all vertices (default -1)

        """
        import numpy as np

        # Make sure argument is a numpy array
        if not isinstance(regions, np.ndarray):
            regions = np.array(regions)

        # Identify the vertices for each hole
        hole_numbers = [x for x in np.unique(holes) if x > -1]
        for n_hole in hole_numbers:
            I = [i for i,x in enumerate(holes) if x == n_hole]

            # Identify neighbors to these vertices
            N=[]; [N.extend(neighbor_lists[i]) for i in I]
            if N:

                # Assign the hole the maximum region ID number of its neighbors
                regions[I] = max([regions[x] for x in N])

        return regions

    #-------------------------------------------------------------------------
    # Find boundaries to holes
    #-------------------------------------------------------------------------
    hole_boundaries = -1 * np.ones(len(regions))

    # Identify vertices for each region
    region_numbers = [x for x in np.unique(regions) if x > -1]
    count = 0
    for n_region in region_numbers:
        region_indices = np.where(regions == n_region)[0]

        # Identify neighbors to these vertices and their neighbors
        N = []
        [N.extend(neighbor_lists[x]) for x in region_indices]
        N = list(frozenset(N).difference(region_indices))
        N2 = []
        [N2.extend(neighbor_lists[x]) for x in N]
        N.extend(N2)
        N = list(frozenset(N).difference(region_indices))
        if N:

            # Segment the neighbors into connected vertices (region boundaries)
            boundaries = segment(N, neighbor_lists)

            # Remove the largest region boundary, presumably the
            # outer contour of the region, leaving smaller boundaries,
            # presumably the contours of holes within the region
            boundary_numbers = [x for x in np.unique(boundaries) if x > -1]
            max_size = 0
            max_number = 0
            for n_boundary in boundary_numbers:
                boundary_indices = np.where(boundaries == n_boundary)[0]
                if len(boundary_indices) > max_size:
                    max_size = len(boundary_indices)
                    max_number = n_boundary
            boundaries[boundaries == max_number] = -1
            boundary_numbers = [x for x in boundary_numbers if x != max_number]

            # Add remaining boundaries to holes array
            for n_boundary in boundary_numbers:
                indices = [i for i,x in enumerate(boundaries) if x == n_boundary]
                hole_boundaries[indices] = count
                count += 1

    #-------------------------------------------------------------------------
    # Fill holes
    #-------------------------------------------------------------------------
    # If there are any holes
    if count > 0:
        hole_numbers = [x for x in np.unique(hole_boundaries) if x > -1]
        background = [i for i,x in enumerate(regions) if x == -1]

        # Grow seeds from hole boundaries to fill holes (background value -1)
        for n_hole in hole_numbers:
            seed_list = np.where(hole_boundaries == n_hole)[0].tolist()
            seed_lists = [list(frozenset(background).intersection(seed_list))]
            hole = segment(background, neighbor_lists, 1, seed_lists)

            # Label the vertices for each hole by surrounding region number
            # if hole does not include values within exclude_range:
            if len(exclude_range) == 2:
                Ihole = np.where(hole > -1)[0]
                #if not len(frozenset(values[Ihole]).intersection(exclude_range)):
                if not [x for x in values[Ihole]
                        if x > exclude_range[0] and x < exclude_range[1]]:
                    regions = label_holes(hole, regions, neighbor_lists)
            else:
                regions = label_holes(hole, regions, neighbor_lists)

    return regions

#-----------------------------------------------------------------------------
# Test for simple points
#-----------------------------------------------------------------------------
def topo_test(index, values, neighbor_lists):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    Parameters
    ----------
    index : integer
        index of vertex
    values : numpy array of integers or floats
        values for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    sp : Boolean
        simple point or not?
    n_inside : integer
        number of neighboring vertices greater than threshold

    """
    import numpy as np

    # Make sure argument is a numpy array
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    # Find neighbors to the input vertex, and binarize them
    # into those greater or less than the class boundary threshold for HMMF (0.5)
    # ("inside" and "outside"); count the number of inside and outside neighbors
    I_neighbors = neighbor_lists[index]
    neighbor_values = values[I_neighbors]
    inside = [I_neighbors[i] for i,x in enumerate(neighbor_values) if x > 0.5]
    n_inside = len(inside)
    n_outside = len(I_neighbors) - n_inside

    # If the number of inside or outside neighbors is zero,
    # than the vertex IS NOT a simple point
    if n_outside * n_inside == 0:
        sp = False
    # Or if either the number of inside or outside neighbors is one,
    # than the vertex IS a simple point
    elif n_outside == 1 or n_inside == 1:
        sp = True
    # Otherwise, test to see if all of the inside neighbors share neighbors
    # with each other, in which case the vertex IS a simple point
    else:
        # For each neighbor exceeding the threshold,
        # find its neighbors that also exceed the threshold,
        # then store these neighbors' indices in a sublist of "N"
        labels = range(1, n_inside + 1)
        N = []
        for i_in in range(n_inside):
            new_neighbors = neighbor_lists[inside[i_in]]
            new_neighbors = [x for x in new_neighbors
                             if values[x] > 0.5 if x != index]
            new_neighbors.extend([inside[i_in]])
            N.append(new_neighbors)

        # Consolidate labels of connected vertices:
        # Loop through neighbors (lists within "N"),
        # reassigning the labels for the lists until each label's
        # list(s) has a unique set of vertices
        change = True
        while change:
            change = False

            # Loop through pairs of inside neighbors
            # and continue if their two labels are different
            for i in range(n_inside - 1):
                for j in range(i + 1, n_inside):
                    if labels[i] != labels[j]:
                        # Assign the two subsets the same label
                        # if they share at least one vertex,
                        # and continue looping
                        if frozenset(N[i]).intersection(N[j]):
                            labels[i] = max([labels[i], labels[j]])
                            labels[j] = labels[i]
                            change = True

        # The vertex is a simple point if all of its neighbors
        # (if any) share neighbors with each other (one unique label)
        D = []
        if len([D.append(x) for x in labels if x not in D]) == 1:
            sp = True
        else:
            sp = False

    return sp, n_inside
