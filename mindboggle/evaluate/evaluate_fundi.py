#!/usr/bin/python

"""
Evaluate fundi by computing closest distances from label boundary vertices
corresponding to fundi in the Desikan-Killiany-Tourville cortical labeling
protocol.

Example
-------
$ python evaluate_fundi.py fundi.vtk folds.vtk lh.labels.DKT25.manual.vtk


Authors:
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def compute_fundus_distances(label_boundary_fundi, fundi, folds, points, n_fundi):
    """
    Create a fundus distance matrix.

    Compute the minimum distance from each label boundary vertex corresponding
    to a fundus in the Desikan-Killiany-Tourville cortical labeling protocol
    to all of the fundus vertices in the same fold.

    Returns
    -------
    distances : numpy array
        distance value for each vertex (zero where there is no vertex)
    distance_matrix : numpy array
        rows are for vertices and columns for fundi (-1 default value)

    """
    import numpy as np
    from mindboggle.measure.measure_functions import compute_point_distance

    n_points = len(points)

    distances = np.zeros(n_points)
    distance_matrix = -1 * np.ones((n_points, n_fundi))

    sum_distances = 0
    num_distances = 0

    # For each label boundary fundus point
    for i_label_point, fundus_ID in enumerate(label_boundary_fundi):
        if fundus_ID > 0:

            # Find (indices of) fundus points in the same fold
            I_fundus_points = [i for i,x in enumerate(fundi)
                               if x > 0
                               if folds[i] == folds[i_label_point]]

            # Find the closest fundus point to the label boundary fundus point
            d, i = compute_point_distance(points[i_label_point],
                                          points[I_fundus_points])
            distances[i_label_point] = d
            distance_matrix[i_label_point, fundus_ID - 1] = d

            sum_distances += d
            num_distances += 1

            if i_label_point % 1000 == 0 and num_distances > 0:
                percent_done = 100.0 * float(i_label_point) / float(n_points)
                mean_distance = sum_distances / num_distances
                print('Done: {0} pct, mean dist {1}, {2}'.format(
                      percent_done, mean_distance, i_label_point))

    mean_distance = sum_distances / num_distances
    print('Done: 100 pct, mean dist {0}'.format(mean_distance))

    return distances, distance_matrix


if __name__ == "__main__":

    import sys
    import numpy as np
    from mindboggle.utils.io_vtk import load_scalars, write_scalar_lists
    from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    from mindboggle.utils.mesh_operations import detect_boundaries, find_neighbors

    fundi_file = sys.argv[1]
    folds_file = sys.argv[2]
    labels_file = sys.argv[3]

    print('***')
    print('Input fundi:' + fundi_file)
    print('Input folds:' + folds_file)
    print('Input labels:' + labels_file)
    print('***')

    # Load fundi, folds, labels
    points, faces, fundi, n_vertices = load_scalars(fundi_file, return_arrays=True)
    points, faces, folds, n_vertices = load_scalars(folds_file, return_arrays=True)
    points, faces, labels, n_vertices = load_scalars(labels_file, return_arrays=True)
    n_points = len(points)

    # List of indices to fold vertices
    fold_indices = [i for i,x in enumerate(folds) if x > 0]

    # Calculate neighbor lists for all points
    print('Find neighbors to all vertices...')
    neighbor_lists = find_neighbors(faces, n_points)

    # Prepare list of all unique sorted label pairs in the labeling protocol
    print('Prepare a list of unique, sorted label pairs in the protocol...')
    label_pair_lists = sulcus_boundaries()
    n_fundi = len(label_pair_lists)

    # Find label boundary points in any of the folds
    print('Find label boundary points in any of the folds...')
    boundary_indices, boundary_label_pairs, unique_boundary_label_pairs = \
        detect_boundaries(fold_indices, labels, neighbor_lists)
    if not len(boundary_indices):
        sys.exit('There are no label boundary points!')

    # Initialize an array of label boundaries fundus IDs
    # (label boundary vertices that define sulci in the labeling protocol)
    print('Build an array of label boundary fundus IDs...')
    label_boundary_fundi = np.zeros(n_points)

    # For each list of sorted label pairs (corresponding to a sulcus)
    for isulcus, label_pairs in enumerate(label_pair_lists):
        print('  Sulcus ' + str(isulcus + 1))

        # Keep the boundary points with label pair labels
        fundus_indices = [x for i,x in enumerate(boundary_indices)
                          if np.unique(boundary_label_pairs[i]).tolist()
                              in label_pairs]

        # Store the points as sulcus IDs in the fundus IDs array
        label_boundary_fundi[fundus_indices] = isulcus + 1

    # Construct a distance matrix
    print('Construct a boundary-to-fundus distance matrix...')
    distances, distance_matrix = compute_fundus_distances(label_boundary_fundi,
        fundi, folds, points, n_fundi)

    # Compute mean distances for each fundus
    fundus_mean_distances = np.zeros(n_fundi)
    fundus_std_distances = np.zeros(n_fundi)
    for ifundus in range(n_fundi):
        fundus_distances = distance_matrix[:, ifundus]
        fundus_distances = [x for x in fundus_distances if x != -1]
        fundus_mean_distances[ifundus] = np.mean(fundus_distances)
        fundus_std_distances[ifundus] = np.std(fundus_distances)

    print('Fundus mean distances')
    print(fundus_mean_distances)
    print('Fundus standard deviations of distances')
    print(fundus_std_distances)

    # Write resulting fundus-label boundary distances to VTK file
    print('Write distances to a VTK file for visualization...')
    write_scalar_lists('evaluate_fundi.vtk', points, range(n_points), faces,
                       [distances], ['fundus-label boundary distances'])
