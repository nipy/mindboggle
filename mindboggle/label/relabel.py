#!/usr/bin/python

"""
Relabel surface or volume or annot files.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def relabel_volume(input_file, old_labels, new_labels):
    """
    Relabel volume labels in a nifti file.

    Parameters
    ----------
    input_file : labeled nibabel-readable (e.g., nifti) file
    old_labels : lists of old labels
    new_labels : lists of new labels

    """

    import os
    import numpy as np
    import nibabel as nb

    # Load labeled image volume and extract data as 1-D array
    vol = nb.load(input_file)
    xfm = vol.get_affine()
    data = vol.get_data().ravel()

    # Initialize output
    new_data = data.copy()

    # Loop through labels
    for ilabel, label in enumerate(old_labels):
        label = int(label)
        new_label = int(new_labels[ilabel])

        # Relabel
        if new_label != label:
            new_data[np.where(data==label)[0]] = new_label

    # Reshape to original dimensions
    new_data = np.reshape(new_data, vol.shape)

    # Save relabeled file
    output_file = os.path.join(os.getcwd(), os.path.basename(input_file))
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    return output_file

def relabel_surface(vtk_file, relabel_list, new_string):
    """
    Relabel surface in a VTK file.

    Parameters
    ----------
    vtk_file : input labeled VTK file
    relabel_list : text file with two columns of label numbers --
                   all regions receive the 2nd label per row.
    new_string : new ending of vtk_file name (e.g., 'labels.DKT25.vtk')

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
    labels_to_replace, new_labels = io_file.read_columns(relabel_list, 2)
    for i, new_label in enumerate(new_labels):

        # Find which vertices have the label
        indices = np.where(Scalars == int(labels_to_replace[i]))[0]
        Scalars[indices] = int(new_label)

    relabeled_vtk = os.path.join(os.getcwd(),
                                 os.path.basename(vtk_file).split('.')[0] + \
                                 '.' + new_string)
    io_vtk.write_scalars(relabeled_vtk, Points, Vertices, Faces,
                         [Scalars.tolist()], ['Labels'])

    return relabeled_vtk

def relabel_annot_file(hemi, subject, annot_name, new_annot_name, relabel_file):
    """
    Combine surface labels in a .annot file.

    https://mail.nmr.mgh.harvard.edu/pipermail//freesurfer/2010-June/014620.html
      
     `mris_translate_annotation <subject> <hemi> <in annot> <translation file> <out annot>`
      
      ``translation file``: text file that lists the labels (one per line)
      you want to group, and the new label you want to create.  You have to use
      the RGB codes; each line will provide the input and output RGB values::

            221     220     60      223     220     60
            221     220     160     223     220     60
            221     220     100     223     220     60

    Parameters
    ----------
    hemi : hemisphere [string]
    subject : subject name
    annot_name : name of .annot file (without pre- or post-pends)
    relabel_file : text file with old and new RGB values
    new_annot_name : new .annot name

    """

    import os
    from nipype.interfaces.base import CommandLine

    cli = CommandLine(command='mris_translate_annotation')
    cli.inputs.args = ' '.join([subject, hemi, annot_name, relabel_file, new_annot_name])
    cli.cmdline
    cli.run()

    return new_annot_name
