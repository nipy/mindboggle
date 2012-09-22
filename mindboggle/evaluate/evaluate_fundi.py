#!/usr/bin/python

"""
Evaluate fundi by computing closest distances from label boundary vertices
corresponding to fundi in the Desikan-Killiany-Tourville cortical labeling
protocol.

Example
-------

$ python evaluate_fundi.py <fundi_file> <folds_file> <labels_file>


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

    In the resulting distance matrix, the leftmost column contains distances
    between label boundaries and fundi for each point, and subsequent columns
    represent fundi, with the distance values repeated only for fundus points
    (values default to -1).

    """
    import numpy as np
    from utils.mesh_operations import closest_distance

    n_points = len(points)

    distances = np.zeros(n_points)
    distance_matrix = -1 * np.ones((n_points, n_fundi))

    sum_distances = 0
    num_distances = 0

    # For each label boundary fundus point
    for i_label_point, fundus_ID in enumerate(label_boundary_fundi):
        if fundus_ID > 0:

            # Find indices of fundus points in the same fold
            I_fundus_points = [i for i,x in enumerate(fundi)
                               if x > 0
                               if folds[i] == folds[i_label_point]]

            # Find the closest fundus point to label boundary fundus point
            d = closest_distance(points[i_label_point],
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
    from utils.io_vtk import load_scalar
    from info.sulcus_boundaries import sulcus_boundaries
    from utils.mesh_operations import detect_boundaries, find_all_neighbors

    fundi_file = sys.argv[1]
    folds_file = sys.argv[2]
    labels_file = sys.argv[3]

    print('***')
    print('Input fundi:' + fundi_file)
    print('Input folds:' + folds_file)
    print('Input labels:' + labels_file)
    print('***')

    # Load fundi, folds, labels
    points, faces, fundi = load_scalar(fundi_file, return_arrays=1)
    points, faces, folds = load_scalar(folds_file, return_arrays=1)
    points, faces, labels = load_scalar(labels_file, return_arrays=1)
    n_points = len(points)

    # List of lists of folds
    unique_fold_IDs = set(folds)
    fold_lists = []
    for id in unique_fold_IDs:
        if id > 0:
            fold = [i for i,x in enumerate(folds) if x == id]
            fold_lists.append(fold)

    # Calculate neighbor lists for all points
    neighbor_lists = find_all_neighbors(faces)

    # Prepare list of all unique sorted label pairs in the labeling protocol
    label_pair_lists = sulcus_boundaries()
    n_fundi = len(label_pair_lists)

    # Initialize an array of label boundaries fundus IDs
    # (label boundary vertices that define sulci in the labeling protocol)
    label_boundary_fundi = np.zeros(n_points)

    # For each list of label pairs (corresponding to a sulcus)
    for isulcus, label_pairs in enumerate(label_pair_lists):

        # For each fold
        for fold in fold_lists:

            # Find label boundary points within the fold
            boundary_points, boundary_labels = detect_boundaries(labels,
                 fold, neighbor_lists)

            # Keep the boundary points with label pair labels
            fundus_points = [x for i,x in enumerate(boundary_points)
                             if boundary_labels[i] in label_pairs]

            # Store the points as sulcus IDs in the fundus IDs array
            label_boundary_fundi[fundus_points] = isulcus + 1

    # Create a distance matrix:
    # The leftmost column contains distances between label boundaries
    # and fundi for each point, and subsequent columns represent fundi,
    # with the distance values repeated only for fundus points (else -1).
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
    write_scalars('evaluate_fundi.vtk', points, range(n_points), faces,
                  [distances], ['fundus-label boundary distances'])
