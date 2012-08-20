#!/usr/bin/python

"""
Combine surface labels in a VTK file.

For doing the same in a .annot file:
- https://mail.nmr.mgh.harvard.edu/pipermail//freesurfer/2010-June/014620.html
- mris_translate_annotation <subject> <hemi> <in annot> <translation file> <out annot>
- <translation file>: text file that lists the labels (one per line)
  you want to group, and the new label you want to create.  You have to use
  the RGB codes; each line will provide the input and output RGB values:
  221     220     60      223     220     60
  221     220     160     223     220     60
  221     220     100     223     220     60


Author:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def relabel_surface(vtk_file, relabel_list, new_string):
    """
    Combine surface labels.

    Input:
    vtk_file:  input labeled VTK file
    relabel_list:  text file with two columns of label numbers --
                   all regions receive the 2nd label per row.
    new_string:  new ending of vtk_file name (e.g., 'labels.DKT25.vtk')

    """

    import os
    import numpy as np
    import utils.io_vtk as io_vtk
    import utils.io_file as io_file

    # Load labeled vtk surfaces
    Points, Faces, Scalars = io_vtk.load_scalar(vtk_file)
    Scalars = np.array(Scalars)
    Vertices = range(1, len(Points) + 1)

    # Load label lists
    label_lists = io_file.read_list_2strings(relabel_list)
    labels_to_replace = label_lists[0]
    new_labels = label_lists[1]
    for i, new_label in enumerate(new_labels):

        # Find which vertices have the label
        indices = np.where(Scalars == int(labels_to_replace[i]))[0]
        Scalars[indices] = int(new_label)

    relabeled_vtk = os.path.join(os.getcwd(),
                                 os.path.basename(vtk_file).split('.')[0] + \
                                 '.' + new_string)
    print(relabeled_vtk)
    io_vtk.write_scalars(relabeled_vtk, Points, Vertices, Faces,
                         [Scalars.tolist()], ['Labels'])

    return relabeled_vtk
