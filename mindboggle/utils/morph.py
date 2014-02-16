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
        indices of original vertices with dilated vertices

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.morph import dilate
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
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
    >>> rewrite_scalars(vtk_file, 'dilate.vtk', IDs, 'dilated_fold', IDs)
    >>> plot_surfaces('dilate.vtk')

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
        indices of original vertices without eroded vertices

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.morph import erode
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
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
    >>> rewrite_scalars(vtk_file, 'erode.vtk', IDs, 'eroded_fold', IDs)
    >>> plot_surfaces('erode.vtk')

    """
    from mindboggle.utils.mesh import find_neighborhood

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
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.morph import extract_edge
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> # Select a single fold:
    >>> fold_number = 11 #11
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> #
    >>> edge_indices = extract_edge(indices, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> IDs = -1 * np.ones(len(folds))
    >>> IDs[indices] = 1
    >>> IDs[edge_indices] = 2
    >>> rewrite_scalars(vtk_file, 'extract_edge.vtk', IDs, 'edge', IDs)
    >>> plot_surfaces('extract_edge.vtk')

    """
    from mindboggle.utils.mesh import find_neighborhood

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
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors, remove_faces
    >>> from mindboggle.utils.morph import fill_holes
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> #
    >>> background_value = -1
    >>> # Select one fold
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> faces, lines, indices, points, npoints, scalars, name, input_vtk = read_vtk(vtk_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> n_fold = np.unique(folds)[1]
    >>> folds[folds != n_fold] = background_value
    >>> #
    >>> # Make two holes in fold (background value and excluded values)
    >>> # Hole 1:
    >>> # Find a vertex whose removal (with its neighbors) would create a hole
    >>> I = np.where(folds==n_fold)[0]
    >>> for index1 in I:
    >>>     N1 = neighbor_lists[index1]
    >>>     stop = True
    >>>     for n in N1:
    >>>         if any(folds[neighbor_lists[n]] == background_value):
    >>>             stop = False
    >>>             break
    >>>         else:
    >>>             for f in neighbor_lists[n]:
    >>>                 if any(folds[neighbor_lists[f]] == background_value):
    >>>                     stop = False
    >>>                     break
    >>>     if stop:
    >>>         break
    >>> folds[index1] = background_value
    >>> folds[N1] = background_value
    >>> # Hole 2:
    >>> I = np.where(folds==n_fold)[0]
    >>> for index2 in I:
    >>>     N2 = neighbor_lists[index2]
    >>>     stop = True
    >>>     for n in N2:
    >>>         if any(folds[neighbor_lists[n]] == background_value):
    >>>             stop = False
    >>>             break
    >>>         else:
    >>>             for f in neighbor_lists[n]:
    >>>                 if any(folds[neighbor_lists[f]] == background_value):
    >>>                     stop = False
    >>>                     break
    >>>     if stop:
    >>>         break
    >>> folds[index2] = background_value
    >>> folds[N2] = background_value
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
    >>> indices = [i for i,x in enumerate(holes) if x != background_value]
    >>> write_vtk('holes.vtk', points, indices, lines,
    >>>           remove_faces(faces, indices), [holes.tolist()], ['holes'], 'int')
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> plot_surfaces('holes.vtk')
    >>> #
    >>> # Fill Hole 1 but not Hole 2:
    >>> # (because values has an excluded value in the hole)
    >>> regions = np.copy(folds)
    >>> exclude_range = [99,101],
    >>> regions = fill_holes(regions, neighbor_lists, values, exclude_range, background_value)
    >>> #
    >>> # Write results to vtk file and view:
    >>> indices = [i for i,x in enumerate(regions) if x != background_value]
    >>> write_vtk('fill_holes.vtk', points, indices, lines,
    >>>           remove_faces(faces, indices), regions.tolist(), 'regions', 'int')
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> plot_surfaces('fill_holes.vtk')

    """
    import numpy as np
    from mindboggle.utils.segment import segment

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
    >>> # Example 1: build a cube by closing two parallel planes:
    >>> import os
    >>> from mindboggle.utils.morph import close_surface_pair
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> from mindboggle.utils.io_vtk import write_vtk
    >>> # Build plane:
    >>> background_value = -1
    >>> n = 10  # plane edge length
    >>> points1 = []
    >>> for x in range(n):
    >>>     for y in range(n):
    >>>         points1.append([x,y,0])
    >>> points2 = [[x[0],x[1],1] for x in points1]
    >>> scalars = [background_value for x in range(len(points1))]
    >>> p = n*(n-1)/2 - 1
    >>> for i in [p, p+1, p+n, p+n+1]:
    >>>     scalars[i] = 1
    >>> faces = []
    >>> for x in range(n-1):
    >>>     for y in range(n-1):
    >>>         faces.append([x+y*n,x+n+y*n,x+n+1+y*n])
    >>>         faces.append([x+y*n,x+1+y*n,x+n+1+y*n])
    >>> #write_vtk('plane.vtk', points1, [], [], faces, scalars)
    >>> #plot_surfaces('plane.vtk') # doctest: +SKIP
    >>> closed_faces, closed_points, closed_scalars = close_surface_pair(faces, points1, points2, scalars, background_value)
    >>> # View:
    >>> write_vtk('cube.vtk', closed_points, [], [], closed_faces, closed_scalars, 'int')
    >>> plot_surfaces('cube.vtk') # doctest: +SKIP
    >>> #
    >>> # Example 2: Gray and white cortical brain surfaces:
    >>> import os
    >>> from mindboggle.utils.morph import close_surface_pair
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, read_points, write_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> patch_surface1 = 'fold.pial.vtk'
    >>> whole_surface2 = 'fold.white.vtk'
    >>> # Select a single fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> points1 = read_points(folds_file)
    >>> scalars, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> scalars[scalars != fold_number] = -1
    >>> white_surface = os.path.join(path, 'arno', 'freesurfer', 'lh.white.vtk')
    >>> faces, u1, u2, points2, N, u3, u4, u5 = read_vtk(white_surface)
    >>> background_value = -1
    >>> closed_faces, closed_points, closed_scalars = close_surface_pair(faces, points1, points2, scalars, background_value)
    >>> # View:
    >>> write_vtk('closed.vtk', closed_points, [], [], closed_faces, closed_scalars, name, 'int')
    >>> plot_surfaces('closed.vtk') # doctest: +SKIP

    """
    import sys
    import numpy as np

    from mindboggle.utils.mesh import find_neighbors, remove_faces
    from mindboggle.utils.segment import extract_borders

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
    faces = remove_faces(faces, I)
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


def close_surface_pair_from_files(patch_surface1, whole_surface2,
                                  background_value=-1, output_vtk=''):
    """
    Close a surface patch by connecting its border vertices with
    corresponding vertices in a second surface file.

    Assumes no lines or indices when reading VTK files in.

    Note ::

        The first VTK file contains scalar values different than background
        for a surface patch.  The second VTK file contains the (entire)
        surface whose corresponding vertices are shifted in position.
        For pial vs. gray-white matter, the two surfaces are not parallel,
        so connecting the vertices leads to intersecting faces.

    Parameters
    ----------
    patch_surface1 : string
        vtk file with surface patch of non-background scalar values
    whole_surface2 : string
        second vtk file with 1-to-1 vertex correspondence with patch_surface1
        (whole surface so as to derive vertex neighbor lists)
    background_value : integer
        scalar value for background vertices
    output_vtk : string
        output vtk file name with closed surface patch

    Returns
    -------
    output_vtk : string
        output vtk file name with closed surface patch

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.morph import close_surface_pair_from_files
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, read_points, write_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> patch_surface1 = 'fold.pial.vtk'
    >>> whole_surface2 = 'fold.white.vtk'
    >>> # Select a single fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> points = read_points(folds_file)
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> folds[folds != fold_number] = -1
    >>> white_surface = os.path.join(path, 'arno', 'freesurfer', 'lh.white.vtk')
    >>> faces, u1, u2, points2, N, u3, u4, u5 = read_vtk(white_surface)
    >>> write_vtk(patch_surface1, points, [], [], faces, folds, name)
    >>> write_vtk(whole_surface2, points2, [], [], faces, folds, name)
    >>> background_value = -1
    >>> output_vtk = ''
    >>> close_surface_pair_from_files(patch_surface1, whole_surface2, background_value, output_vtk)
    >>> # View:
    >>> plot_surfaces('closed.vtk') # doctest: +SKIP

    """
    import os
    import numpy as np

    from mindboggle.utils.io_vtk import read_vtk, write_vtk
    from mindboggle.utils.morph import close_surface_pair

    # Read VTK surface mesh files:
    u1, u2, u3, points1, N, scalars, name, u4 = read_vtk(patch_surface1,
                                                         True, True)
    faces, u1, u2, points2, N, u3, u4, u5 = read_vtk(whole_surface2,
                                                     True, True)

    # Close surface:
    closed_faces, closed_points, closed_scalars = close_surface_pair(faces,
        points1, points2, scalars, background_value)

    # Write output file:
    if not output_vtk:
        output_vtk = os.path.join(os.getcwd(), 'closed.vtk')
    # closed_scalars is a list
    if np.ndim(closed_scalars) == 1:
        scalar_type = type(closed_scalars[0]).__name__
    elif np.ndim(closed_scalars) == 2:
        scalar_type = type(closed_scalars[0][0]).__name__
    else:
        print("Undefined scalar type!")
    write_vtk(output_vtk, closed_points, [], [], closed_faces, closed_scalars,
              name, scalar_type=scalar_type)

    return output_vtk


#-----------------------------------------------------------------------------
# Test for simple points
#-----------------------------------------------------------------------------
def topo_test(index, values, neighbor_lists):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    "Simple" is not to be mistaken with the following usage:
    "A vertex is usually assigned one of five possible classifications:
    simple, complex, boundary, interior edge, or corner vertex.
     A simple vertex is surrounded by a closed fan of triangles".

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
