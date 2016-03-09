#!/usr/bin/env python
"""
Functions for realigning label boundaries based on fundus lines.

Authors:
    - Oliver Hinds, 2013  (ohinds@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def realign_boundaries_to_fundus_lines(
    surf_file, init_label_file, fundus_lines_file, thickness_file,
    out_label_file=None):
    """
    Fix label boundaries to fundus lines.

    Parameters
    ----------
    surf_file : file containing the surface geometry in vtk format
    init_label_file : file containing scalars that represent the
                      initial guess at labels
    fundus_lines_file : file containing scalars representing fundus lines.
    thickness_file: file containing cortical thickness scalar data
    (for masking out the medial wall only)
    out_label_file : if specified, the realigned labels will be writen to
                     this file

    Returns
    -------
    numpy array representing the realigned label for each surface vertex.
    """

    import numpy as np
    from mindboggle.guts.segment import extract_borders
    import mindboggle.guts.graph as go
    from mindboggle.mio.vtks import read_vtk, read_scalars, write_vtk
    from mindboggle.guts.mesh import find_neighbors
    import propagate_fundus_lines

    ## read files
    points, indices, lines, faces, scalars, scalar_names, num_points, \
        input_vtk = read_vtk(surf_file, return_first=True, return_array=True)
    indices = range(num_points)

    init_labels, _ = read_scalars(init_label_file,
                                  return_first=True, return_array=True)

    fundus_lines, _ = read_scalars(fundus_lines_file,
                                   return_first=True, return_array=True)

    thickness, _ = read_scalars(thickness_file,
                             return_first=True, return_array=True)

    # remove labels from vertices with zero thickness (get around
    # DKT40 annotations having the label '3' for all the Corpus
    # Callosum vertices).
    cc_inds = [x for x in indices if thickness[x] < 0.001]
    init_labels[cc_inds] = 0

    ## setup seeds from initial label boundaries
    neighbor_lists = find_neighbors(faces, num_points)

    # extract all vertices that are on a boundary between labels
    boundary_indices, label_pairs, _ = extract_borders(
        indices, init_labels, neighbor_lists,
        return_label_pairs=True)

    # split boundary vertices into segments with common boundary pairs.
    boundary_segments = {}
    for boundary_index, label_pair in zip(boundary_indices, label_pairs):
        key = ((label_pair[0], label_pair[1]) if label_pair[0] < label_pair[1]
               else (label_pair[1], label_pair[0]))
        if key not in boundary_segments:
            boundary_segments[key] = []

        boundary_segments[key].append(boundary_index)

    boundary_matrix, boundary_matrix_keys = _build_boundary_matrix(
        boundary_segments, num_points)

    # build the affinity matrix
    affinity_matrix = go.weight_graph(np.array(points), indices,
        np.array(faces), sigma=10, add_to_graph=False, verbose=False)

    ## propagate boundaries to fundus line vertices
    learned_matrix = _propagate_labels(
       affinity_matrix, boundary_matrix, boundary_indices, 100, 1)

    # assign labels to fundus line vertices based on highest probability
    new_boundaries = -1 * np.ones(init_labels.shape)
    fundus_line_indices = [i for i, x in enumerate(fundus_lines) if x > 0.5]

    # tile the surface into connected components delimited by fundus lines
    closed_fundus_lines, _, _ = propagate_fundus_lines.propagate_fundus_lines(
        points, faces, fundus_line_indices, thickness)

    closed_fundus_line_indices = np.where(closed_fundus_lines > 0)[0]

    # split surface into connected components
    connected_component_faces = _remove_boundary_faces(
        points, faces, closed_fundus_line_indices)

    # label components based on most probable label assignment
    new_labels = _label_components(
        connected_component_faces, num_points, boundary_indices, learned_matrix,
        boundary_matrix_keys)

    # propagate new labels to fill holes
    label_matrix, label_map = _build_label_matrix(new_labels)
    new_learned_matrix = _propagate_labels(
        affinity_matrix, label_matrix,
        [i for i in range(num_points) if new_labels[i] >= 0], 100, 1)

    # assign most probable labels
    for idx in [i for i in range(num_points) if new_labels[i] == -1]:
        max_idx = np.argmax(new_learned_matrix[idx])
        new_labels[idx] = label_map[max_idx]

    # save
    if out_label_file is not None:
        write_vtk(out_label_file, points, faces=faces,
                  scalars=[int(x) for x in new_labels], scalar_type='int')

    return new_labels

def _build_boundary_matrix(boundary_segments, num_points):
    import numpy as np

    boundary_segment_matrix = np.zeros((num_points, len(boundary_segments)))

    realignment_mapping = {}
    for label_index, (key, val) in enumerate(list(boundary_segments.items())):
        realignment_mapping[label_index] = key
        boundary_segment_matrix[val, :] = -1
        boundary_segment_matrix[val, label_index] = 1

    return boundary_segment_matrix, realignment_mapping

def _build_label_matrix(vertex_labels):
    import numpy as np

    label_map = sorted(set(vertex_labels))

    label_matrix = np.zeros((len(vertex_labels), len(label_map)))

    for idx, label in enumerate(vertex_labels):
        if label >= 0:
            label_matrix[idx, :] = -1
            label_index = np.where(label_map == label)[0]
            label_matrix[idx, label_index] = 1

    return label_matrix, label_map

def _propagate_labels(affinity_matrix, label_matrix, seed_indices,
                      max_iterations, tolerance):
    """
    This function was lifted from labels/rebound.py.

    Run iterative weighted average algorithm to propagate
    labels to unlabeled vertices.

    Parameters
    ----------
    max_iterations:  int (number of iterations)
    tolerance:       float (threshold for terminating algorithm)

    Returns
    -------
    learned_matrix: np array
        (n x C matrix of probabilities that vertex belongs to a given label)

    Features: Hard label clamps, probabilistic solution.
    See: Zhu and Ghahramani, 2002.

    """

    """ The first approach to be considered in the semi-supervised
    learning case is to propagate labels on the graph.  A simple
    algorithm of this sort has been propoosed by Zhu and
    Ghahramani (2002), and starts (like all such algorithms) with
    a set of n vertices, l of which are labeled, and u unlabeled.
    The algorithm takes as its input the affinity matrix W
    (self.affinity_matrix).  From the affinity matrix, one may
    construct the diagonal degree matrix, which is a measure of
    the total weight (or number of edges) attached to a vertex.
    """

    import mindboggle.guts.graph as go
    import numpy as np
    from scipy.sparse import csr_matrix
    import time

    DDM = go.diagonal_degree_matrix(affinity_matrix, inverse=True)


    """ Now, we can actually proceed to perform the iterative
    algorithm.  At each timestep, the labels will be updated to
    reflect the weighted average of adjacent vertices. An important
    caveat of this algorithm is that the labeled vertices remain
    fixed, or clamped.  They should not be changed, and will need to
    be reset.  We accomplish the reset by recalling that
    self.seed_labels stores the indexes of seed labels, and
    self.Labels contains the actual labels.  The algorithm repeats
    itself until either convergence or max_iterations (which will
    prevent excessive computation time).  We must also take care to
    solve the multi-label problem.  To do so, we employ a one-vs-all
    framework, where each label is considered independently, and set
    against the rest of the labels.  More specifically,
    self.label_matrix is an n x C matrix, where each row represents a
    vertex and each column represents label membership.  We can go
    column by column and process the algorithm iteratively. So, we'd
    start at the first column and see which vertices get labeled. Then
    we'd move to the next column and label more vertices.  Because it
    is possible (likely) that some vertices will not receive any
    label, and also to account for probabilistic labeling, we will
    assign a probability of a vertex receiving a label. Then we can
    report these probabilities.  So, to begin, let us first construct
    this probabilistic label assignment: This matrix will store a 1
    for 100% probability, 0 for 0%, and fractional values for the
    rest.  We will rename self.label_matrix for this purpose."""

    learned_matrix = label_matrix

    """ We will later change the -1s to 0s.  As vertices get labeled,
    we assign a confidence measure to the labeling and store the value
    in this matrix.  Now, let us go column by column, and run the
    weighted averaging algorithm.  For each column, you're studying
    one label. Therefore, when updating self.learned_matrix, you'll be
    working with one column at a time too.  If a label gets vertex,
    keep the fractional value, do not simply round to 1 to assign
    membership."""

    i = 0 # record of label number
    for column in learned_matrix.T:

        t0 = time.time()
        print('Number of initial members for label {0}: {1}'.format(
              i, np.nonzero(column==1)[0].size))

        Y_hat_now = csr_matrix(column).transpose()
        converged = False
        counter = 0
        while not converged and counter < max_iterations:
            # column matrix
            Y_hat_next = (DDM * affinity_matrix * Y_hat_now).todense()
            # reset
            Y_hat_next[seed_indices, 0] = column[seed_indices]
            # check convergence
            converged = (np.sum(np.abs(Y_hat_now.todense() - Y_hat_next))
                         < tolerance)
            Y_hat_now = csr_matrix(Y_hat_next)
            counter += 1

        # Print out the number of iterations, so that we get a sense
        # for future runs.  It is also an indication of whether the
        # algorithm converged.

        if counter == max_iterations:
            print('Done in {0:.2f} seconds (the algorithm did not converge)'.
                  format(time.time()-t0))
        else:
            print('Done in {0:.2f} seconds ({1} iterations)'.
                  format(time.time()-t0, counter))

        learned_matrix[:,i] = Y_hat_now.todense().flatten()

        #print('There were {0} initial seed vertices for this label'.
        # format(self.count_assigned_members(i))
        #print('The file actually had {0} vertices for this label'.
        # format(self.count_real_members(self.label_mapping[i]))
        #print('Using only those vertices which crossed the threshold,
        # there are now: '.format(self.count_real_members(i))
        #pylab.plot(self.learned_matrix[:,i])
        #pylab.show()

        i += 1

    """ Before reporting the probabilistic assignment, we change all -1's,
    which indicates 0 probability that the vertex has that label.

    Note: because the labels were initially numbered -1 and 1, there
    will be 'probabilities' below 0.  So, to obtain a sort of
    probability distribution which preserves order, we will add 1 to
    each number, and then divide by 2. Thus -1 --> 0, 1 --> 1 and
    everything else keeps its order."""

    learned_matrix += 1
    learned_matrix = learned_matrix / 2

    """ self.learned_matrix is now complete."""
    return learned_matrix

def _remove_boundary_faces(points, faces, boundary_indices):
    """Given a surface represented by points and faces, remove all
    faces that contain a boundary vertex. Return a new faces list
    without these faces.
    """
    result_faces = []

    for face in faces:
        num_boundary_indices = 0
        for vertex in face:
            if vertex in boundary_indices:
                num_boundary_indices += 1

        if num_boundary_indices == 0:
            result_faces.append(face)

    return result_faces

def _label_components(component_faces, num_points, boundary_indices,
                      boundary_probability_matrix, boundary_matrix_keys):
    """Label the connected components of a surface with the most
    probable label based on the boundary_probability_matrix.
    """

    import numpy as np
    from mindboggle.guts.mesh import find_neighbors

    neighbor_lists = find_neighbors(component_faces, num_points)

    result_labels = -1 * np.ones((num_points))

    # find all the connected components
    point_visited = num_points * [False]
    components = {}
    component_boundaries = {}
    print("Finding connected components")
    while True:
        first_vertex = None
        try:
            first_vertex = next(i for i, v in enumerate(point_visited) if not v)
        except:
            break

        open_vertices = [first_vertex]
        point_visited[first_vertex] = True
        component_vertices = []
        component_boundary_vertices = []
        while len(open_vertices) > 0:
            this_vertex = open_vertices.pop()
            component_vertices.append(this_vertex)
            if this_vertex in boundary_indices:
                component_boundary_vertices.append(this_vertex)
            for neighbor in neighbor_lists[this_vertex]:
                if not point_visited[neighbor]:
                    open_vertices.append(neighbor)
                    point_visited[neighbor] = True
        components[len(component_vertices)] = component_vertices
        component_boundaries[len(component_vertices)] = \
            component_boundary_vertices

    # compute the most probable label for each connected
    # component. Only boundary indices are considered when computing
    # label probability.
    # Note: Here we assume that components and component_boundaries
    # have the same keys.
    used_labels = []
    print("Computing most probable labels")
    #for component in sorted(components.keys(), None, None, True):
    for component in sorted(list(components), None, None, True):
        label_likelihoods = {}
        for vertex in component_boundaries[component]:
            for index, boundary_prob in \
                enumerate(boundary_probability_matrix[vertex]):
                labels = boundary_matrix_keys[index]

                for label in labels:
                    # if label in used_labels:
                    #     continue

                    if label not in label_likelihoods:
                        label_likelihoods[label] = 0

                    label_likelihoods[label] += boundary_prob

        # assign the most likely label
        max_label = None
        max_label_likelihood = None
        for key, val in list(label_likelihoods.items()):
            if max_label is None or val > max_label_likelihood:
                max_label = key
                max_label_likelihood = val

        if max_label is not None:
            result_labels[components[component]] = max_label
            used_labels.append(max_label)

    return result_labels

def main(argv):
    realign_boundaries_to_fundus_lines(argv[1],
                                       argv[2],
                                       argv[3],
                                       argv[4],
                                       argv[5])

if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv))
