#!/usr/bin/python

"""
Aggregate surface labels.


Author:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def aggregate_surface_labels(vtk_file, aggregate_labels_list, label_string):
    """
    Aggregate surface labels.

    Input:
    vtk_file:  input labeled VTK file
    aggregate_labels_file:  text file with two columns of label numbers --
                            all regions receive the 2nd label per row.
    old_string:  old ending of vtk_file name (e.g., 'DKT31.vtk')
    new_string:  new ending of vtk_file name (e.g., 'DKT25.vtk')

    """

    from os import path, getcwd
    import numpy as np
    from utils.io_vtk import load_scalar, write_scalars
    from utils.io_file import read_list_2strings

    # Load labeled vtk surfaces
    Points, Faces, Scalars = load_scalar(vtk_file)
    Scalars = np.array(Scalars)
    Vertices = range(1, len(Points) + 1)

    # Load label lists
    label_lists = read_list_2strings(aggregate_labels_list)
    labels_to_replace = label_lists[0]
    aggregate_labels = label_lists[1]
    for i, aggregate_label in enumerate(aggregate_labels):
            
        # Find which vertices have the label
        indices = np.where(Scalars == int(labels_to_replace[i]))[0]
        Scalars[indices] = int(aggregate_label)

    aggregate_label_file = path.join(getcwd(),
                                     vtk_file.strip('vtk') + label_string)

    write_scalars(aggregate_label_file, Points, Vertices, Faces,
                  [Scalars.tolist()], [label_string])

    return aggregate_labels_vtk
