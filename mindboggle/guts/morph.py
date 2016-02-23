#!/usr/bin/env python
"""
Morphological image operations on surface mesh vertices.

Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
        indices of original vertices with dilated vertices

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.morph import dilate
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_travel_depth'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> nedges = 3
    >>> # Select a single fold:
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> dilated_indices = dilate(indices, nedges, neighbor_lists)
    >>> (len(indices), len(dilated_indices))
    (1065, 1540)

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars
    >>> IDs = -1 * np.ones(len(folds)) # doctest: +SKIP
    >>> IDs[dilated_indices] = 2 # doctest: +SKIP
    >>> IDs[indices] = 1 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'dilate.vtk', IDs, 'dilated_fold', IDs) # doctest: +SKIP
    >>> plot_surfaces('dilate.vtk') # doctest: +SKIP

    """
    from mindboggle.guts.mesh import find_neighborhood

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
        indices of original vertices without eroded vertices

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.morph import erode
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> nedges = 3
    >>> # Select a single fold:
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> eroded_indices = erode(indices, nedges, neighbor_lists)
    >>> (len(indices), len(eroded_indices))
    (1065, 680)

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> IDs = -1 * np.ones(len(folds)) # doctest: +SKIP
    >>> IDs[indices] = 1 # doctest: +SKIP
    >>> IDs[eroded_indices] = 2 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'erode.vtk', IDs, 'eroded_fold', IDs) # doctest: +SKIP
    >>> plot_surfaces('erode.vtk') # doctest: +SKIP

    """
    from mindboggle.guts.mesh import find_neighborhood

    N1 = find_neighborhood(neighbor_lists, indices, nedges=1)
    N2 = find_neighborhood(neighbor_lists, N1, nedges)

    eroded_indices = list(frozenset(indices).difference(N2))

    return eroded_indices

#-----------------------------------------------------------------------------
# Erode
#-----------------------------------------------------------------------------
def extract_edge(indices, neighbor_lists):
    """
    Erode region on a surface mesh to extract the region's edge.

    Parameters
    ----------
    indices : list of integers
        indices of vertices to erode
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    edge_indices : list of integers
        indices of eroded vertices

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.morph import extract_edge
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> # Select a single fold:
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> edge_indices = extract_edge(indices, neighbor_lists)
    >>> (len(indices), len(edge_indices))
    (1065, 131)

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> IDs = -1 * np.ones(len(folds)) # doctest: +SKIP
    >>> IDs[indices] = 1 # doctest: +SKIP
    >>> IDs[edge_indices] = 2 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'edge.vtk', IDs, 'edges_of_fold', IDs) # doctest: +SKIP
    >>> plot_surfaces('edge.vtk') # doctest: +SKIP

    """
    from mindboggle.guts.mesh import find_neighborhood

    N1 = find_neighborhood(neighbor_lists, indices, nedges=1)
    N2 = find_neighborhood(neighbor_lists, N1, nedges=1)
    edge_indices = list(set(N2).intersection(indices))

    return edge_indices

#-----------------------------------------------------------------------------
# Fill holes
#-----------------------------------------------------------------------------
def fill_holes(regions, neighbor_lists, values=[], exclude_range=[],
               background_value=-1):
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
        region numbers for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    values : list of integers
        values for vertices, for use in determining which holes to remove
    exclude_range : list of two floats
        hole is not filled if it contains values within this range
        (prevents cases where surface connected by folds mistaken for holes)
    background_value : integer
        background value

    Returns
    -------
    regions : numpy array of integers
        region numbers for all vertices

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.morph import fill_holes
    >>> from mindboggle.guts.mesh import find_neighbors, keep_faces
    >>> from mindboggle.mio.vtks import read_scalars, read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> background_value = -1
    >>> # Select one fold
    >>> folds, name = read_scalars(folds_file, return_first=True,
    ...                            return_array=True)
    >>> points, indices, lines, faces, f1,f2, npoints, f3 = read_vtk(vtk_file,
    ...     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> n_fold = np.unique(folds)[1]
    >>> folds[folds != n_fold] = background_value
    >>> # Make two holes in fold (background value and excluded values)
    >>> # Hole 1:
    >>> # Find a vertex whose removal (with its neighbors) would create a hole
    >>> I = np.where(folds==n_fold)[0]
    >>> for index1 in I:
    ...     N1 = neighbor_lists[index1]
    ...     stop = True
    ...     for n in N1:
    ...         if any(folds[neighbor_lists[n]] == background_value):
    ...             stop = False
    ...             break
    ...         else:
    ...             for f in neighbor_lists[n]:
    ...                 if any(folds[neighbor_lists[f]] == background_value):
    ...                     stop = False
    ...                     break
    ...     if stop:
    ...         break
    >>> folds[index1] = background_value
    >>> folds[N1] = background_value
    >>> # Hole 2:
    >>> I = np.where(folds==n_fold)[0]
    >>> for index2 in I:
    ...     N2 = neighbor_lists[index2]
    ...     stop = True
    ...     for n in N2:
    ...         if any(folds[neighbor_lists[n]] == background_value):
    ...             stop = False
    ...             break
    ...         else:
    ...             for f in neighbor_lists[n]:
    ...                 if any(folds[neighbor_lists[f]] == background_value):
    ...                     stop = False
    ...                     break
    ...     if stop:
    ...         break
    >>> folds[index2] = background_value
    >>> folds[N2] = background_value
    >>> values = np.zeros(len(folds))
    >>> values[index2] = 100
    >>> values[N2] = 200
    >>> # Write to vtk files:
    >>> holes = folds.copy()
    >>> holes[index1] = 10
    >>> holes[N1] = 20
    >>> holes[index2] = 30
    >>> holes[N2] = 40
    >>> indices = [i for i,x in enumerate(holes) if x != background_value]
    >>> indices[0:10]
    [259, 535, 536, 539, 540, 541, 545, 546, 547, 548]

    View holes (skip test):

    >>> from mindboggle.mio.vtks import write_vtk # doctest: +SKIP
    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> write_vtk('holes.vtk', points, indices, lines,
    ...           keep_faces(faces, indices),
    ...           [holes.tolist()], ['holes'], 'int') # doctest: +SKIP
    >>> plot_surfaces('holes.vtk') # doctest: +SKIP

    >>> # Fill Hole 1 but not Hole 2:
    >>> # (because values has an excluded value in the hole)
    >>> regions = np.copy(folds)
    >>> exclude_range = [99,101]
    >>> regions = fill_holes(regions, neighbor_lists, values, exclude_range,
    ...                      background_value)
    >>> indices2 = [i for i,x in enumerate(regions) if x != background_value]
    >>> filter(lambda x:x not in indices2, indices)
    [925, 926, 1292, 1293, 1294, 1305, 1306]

    View filled hole (one of two holes filled) (skip test):

    >>> write_vtk('fill_holes.vtk', points, indices2, lines,
    ...     keep_faces(faces, indices2), regions.tolist(),
    ...     'regions', 'int') # doctest: +SKIP
    >>> plot_surfaces('fill_holes.vtk') # doctest: +SKIP

    """
    import numpy as np
    from mindboggle.guts.segment import segment

    # Make sure argument is a numpy array
    if not isinstance(regions, np.ndarray):
        regions = np.array(regions)

    def label_holes(holes, regions, neighbor_lists):
        """
        Fill holes in regions on a surface mesh.

        Parameters
        ----------
        holes : list or array of integers
            hole numbers for all vertices
        regions : numpy array of integers
            region numbers for all vertices
        neighbor_lists : list of lists of integers
            each list contains indices to neighboring vertices for each vertex

        Returns
        -------
        regions : numpy array of integers
            region numbers for all vertices

        """
        import numpy as np

        # Make sure argument is a numpy array
        if not isinstance(regions, np.ndarray):
            regions = np.array(regions)

        # Identify the vertices for each hole
        hole_numbers = [x for x in np.unique(holes) if x != background_value]
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
    hole_boundaries = background_value * np.ones(len(regions))

    # Identify vertices for each region
    region_numbers = [x for x in np.unique(regions) if x != background_value]
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
            boundary_numbers = [x for x in np.unique(boundaries)
                                if x != background_value]
            max_size = 0
            max_number = 0
            for n_boundary in boundary_numbers:
                border_indices = np.where(boundaries == n_boundary)[0]
                if len(border_indices) > max_size:
                    max_size = len(border_indices)
                    max_number = n_boundary
            boundaries[boundaries == max_number] = background_value
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
        hole_numbers = [x for x in np.unique(hole_boundaries)
                        if x != background_value]
        background = [i for i,x in enumerate(regions)
                      if x == background_value]

        # Grow seeds from hole boundaries to fill holes
        for n_hole in hole_numbers:
            seed_list = np.where(hole_boundaries == n_hole)[0].tolist()
            seed_lists = [list(frozenset(background).intersection(seed_list))]
            hole = segment(background, neighbor_lists, 1, seed_lists)

            # Label the vertices for each hole by surrounding region number
            # if hole does not include values within exclude_range:
            if len(exclude_range) == 2:
                Ihole = np.where(hole != background_value)[0]
                #if not len(frozenset(values[Ihole]).intersection(exclude_range)):
                if not [x for x in values[Ihole]
                        if x > exclude_range[0] if x < exclude_range[1]]:
                    regions = label_holes(hole, regions, neighbor_lists)
            else:
                regions = label_holes(hole, regions, neighbor_lists)

    return regions


def close_surface_pair(faces, points1, points2, scalars, background_value=-1):
    """
    Close a surface patch by connecting its border vertices with
    corresponding vertices in a second surface file.

    Assumes no lines or indices when reading VTK files in.

    Note ::

        Scalar values different than background define the surface patch.
        The two sets of points have a 1-to-1 mapping; they are from
        two surfaces whose corresponding vertices are shifted in position.
        For pial vs. gray-white matter, the two surfaces are not parallel,
        so connecting the vertices leads to intersecting faces.

    Parameters
    ----------
    faces : list of lists of integers
        each sublist contains 3 indices of vertices that form a face
        on a surface mesh
    points1 : list of lists of floats
        each sublist contains 3-D coordinates of a vertex on a surface mesh
    points2 : list of lists of floats
        points from second surface with 1-to-1 correspondence with points1
    scalars : numpy array of integers
        labels used to find foreground vertices
    background_value : integer
        scalar value for background vertices

    Returns
    -------
    closed_faces : list of lists of integers
        indices of vertices that form a face on the closed surface mesh
    closed_points : list of lists of floats
        3-D coordinates from points1 and points2
    closed_scalars : list of integers
        scalar values for points1 and points2

    Examples
    --------
    >>> # Build a cube by closing two parallel planes:
    >>> from mindboggle.guts.morph import close_surface_pair
    >>> # Build plane:
    >>> background_value = -1
    >>> n = 10  # plane edge length
    >>> points1 = []
    >>> for x in range(n):
    ...     for y in range(n):
    ...          points1.append([x,y,0])
    >>> points2 = [[x[0],x[1],1] for x in points1]
    >>> scalars = [background_value for x in range(len(points1))]
    >>> p = n*(n-1)/2 - 1
    >>> for i in [p, p+1, p+n, p+n+1]:
    ...     scalars[i] = 1
    >>> faces = []
    >>> for x in range(n-1):
    ...     for y in range(n-1):
    ...         faces.append([x+y*n,x+n+y*n,x+n+1+y*n])
    ...         faces.append([x+y*n,x+1+y*n,x+n+1+y*n])
    >>> #write_vtk('plane.vtk', points1, [], [], faces, scalars)
    >>> #plot_surfaces('plane.vtk')
    >>> closed_faces, closed_points, closed_scalars = close_surface_pair(faces,
    ...     points1, points2, scalars, background_value)
    >>> closed_faces[0:4]
    [[44, 54, 55], [44, 45, 55], [144, 154, 155], [144, 145, 155]]

    View cube (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import write_vtk # doctest: +SKIP
    >>> write_vtk('cube.vtk', closed_points, [],[], closed_faces,
    ...     closed_scalars, 'int') # doctest: +SKIP
    >>> plot_surfaces('cube.vtk') # doctest: +SKIP

    """
    import sys
    import numpy as np

    from mindboggle.guts.mesh import find_neighbors, keep_faces
    from mindboggle.guts.segment import extract_borders

    if isinstance(scalars, list):
        scalars = np.array(scalars)

    N = len(points1)
    closed_points = points1 + points2

    # Find all vertex neighbors and surface patch border vertices:
    neighbor_lists = find_neighbors(faces, N)
    I = np.where(scalars != background_value)[0]
    scalars[scalars == background_value] = background_value + 1
    scalars[I] = background_value + 2
    scalars = scalars.tolist()
    borders, u1, u2 = extract_borders(range(N), scalars, neighbor_lists)
    if not len(borders):
        sys.exit('There are no border vertices!')
    borders = [x for x in borders if x in I]

    # Reindex copy of faces and combine with original (both zero-index):
    indices = range(N)
    indices2 = range(N, 2 * N)
    reindex = dict([(index, indices2[i]) for i, index in enumerate(indices)])
    faces = keep_faces(faces, I)
    faces2 = [[reindex[i] for i in face] for face in faces]
    closed_faces = faces + faces2

    # Connect border vertices between surface patches and add new faces:
    add_faces = []
    taken_already = []
    for index in borders:
        if index not in taken_already:
            neighbors = list(set(neighbor_lists[index]).intersection(borders))
            taken_already.append(index)
            #taken_already.extend([index] + neighbors)
            for neighbor in neighbors:
                add_faces.append([index, index + N, neighbor])
                add_faces.append([index + N, neighbor, neighbor + N])
    closed_faces = closed_faces + add_faces

    closed_scalars = scalars * 2

    return closed_faces, closed_points, closed_scalars


def topo_test(index, values, neighbor_lists):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    "Simple" is not to be mistaken with the following usage:
    "A vertex is usually assigned one of five possible classifications:
    simple, complex, boundary, interior edge, or corner vertex.
    A simple vertex is surrounded by a closed fan of triangles".

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
        number of neighboring vertices with a value greater than threshold

    Examples
    --------
    >>> # Square with a center vertex:
    >>> # indices [[0,1,2],[3,4,6],[7,8,9]] = 0 and indices [2,4,6] = 1:
    >>> import numpy as np
    >>> from mindboggle.guts.morph import topo_test
    >>> values = np.array([0,0,1,0,1,0,1,0,0])
    >>> neighbor_lists = [[1,3],[0,2,3,4],[1,4,5],
    ...                   [0,1,4,6],[1,2,3,5,6,7],[2,4,7,8],
    ...                   [3,4,7],[4,5,6,8],[5,7]]
    >>> sps = []
    >>> for index in range(9):
    ...     sp, n_inside = topo_test(index, values, neighbor_lists)
    ...     sps.append(sp)
    >>> sps
    [False, True, True, True, False, True, True, True, False]

    """
    import numpy as np

    # Make sure argument is a numpy array:
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    # Find neighbors to the input vertex, and binarize them
    # into those greater or less than a class boundary threshold equal to 0.5
    # ("inside" and "outside"); count inside and outside neighbors:
    I_neighbors = neighbor_lists[index]
    neighbor_values = values[I_neighbors]
    inside = [I_neighbors[i] for i,x in enumerate(neighbor_values) if x > 0.5]
    n_inside = len(inside)
    n_outside = len(I_neighbors) - n_inside

    # If the number of inside or outside neighbors is zero,
    # than the vertex IS NOT a simple point:
    if n_outside * n_inside == 0:
        sp = False
    # Or if either the number of inside or outside neighbors is one,
    # than the vertex IS a simple point:
    elif n_outside == 1 or n_inside == 1:
        sp = True
    # Otherwise, test to see if all of the inside neighbors share neighbors
    # with each other, in which case the vertex IS a simple point:
    else:
        # For each neighbor exceeding the threshold,
        # find its neighbors that also exceed the threshold,
        # then store these neighbors' indices in a sublist of "N":
        labels = range(1, n_inside + 1)
        N = []
        for i_in in range(n_inside):
            new_neighbors = neighbor_lists[inside[i_in]]
            new_neighbors = [x for x in new_neighbors
                             if values[x] > 0.5 if x != index]
            new_neighbors.extend([inside[i_in]])
            N.append(new_neighbors)

        # Consolidate labels of connected vertices --
        # Loop through neighbors (lists within "N"),
        # reassigning the labels for the lists until each label's
        # list(s) has a unique set of vertices:
        change = True
        while change:
            change = False

            # Loop through pairs of inside neighbors
            # and continue if their two labels are different:
            for i in range(n_inside - 1):
                for j in range(i + 1, n_inside):
                    if labels[i] != labels[j]:
                        # Assign the two subsets the same label
                        # if they share at least one vertex,
                        # and continue looping:
                        if set(N[i]).intersection(N[j]):
                            labels[i] = max([labels[i], labels[j]])
                            labels[j] = labels[i]
                            change = True

        # The vertex is a simple point if all of its neighbors
        # (if any) share neighbors with each other (one unique label):
        D = []
        if len([D.append(x) for x in labels if x not in D]) == 1:
            sp = True
        else:
            sp = False

    return sp, n_inside


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()