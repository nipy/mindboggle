#!/usr/bin/env python
"""
Relabel surface or volume or annot files.


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def relabel_volume(input_file, old_labels, new_labels):
    """
    Relabel volume labels.

    Parameters
    ----------
    input_file : string
        labeled nibabel-readable (e.g., nifti) file
    old_labels : list of integers
        old labels
    new_labels : list of integers
        new labels

    Examples
    --------
    >>> # Convert DKT31 to DKT25 labels
    >>> import os
    >>> from mindboggle.utils.io_table import read_columns
    >>> from mindboggle.labels.relabel import relabel_volume
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> relabel_file = os.path.join(data_path, 'info', 'labels.volume.DKT31to25.txt')
    >>> old_labels, new_labels = read_columns(relabel_file, 2)
    >>> relabel_volume(input_file, old_labels, new_labels)

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

def remove_volume_labels(input_file, labels_to_remove):
    """
    Remove labels from an image volume.

    Parameters
    ----------
    input_file : string
        labeled nibabel-readable (e.g., nifti) file
    labels_to_remove : list of integers
        labels to remove

    Examples
    --------
    >>> # Remove subcortical labels
    >>> import os
    >>> from mindboggle.labels.relabel import remove_volume_labels
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> labels_to_remove = range(1,300) # Remove noncortical (+aseg) labels
    >>> labels_to_remove.extend([1000,1001,2000,2001])
    >>> remove_volume_labels(input_file, labels_to_remove)

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
    for ilabel, label in enumerate(labels_to_remove):
        label = int(label)
        # Relabel
        new_data[np.where(data==label)[0]] = 0

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
    vtk_file : string
        input labeled VTK file
    relabel_list : string
        text file with two columns of label numbers --
        all regions receive the 2nd label per row.
    new_string : string
        new ending of vtk_file name (e.g., 'labels.DKT25.vtk')

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, write_vtk
    from mindboggle.utils.io_table import read_columns

    # Load labeled vtk surfaces
    faces, lines, indices, points, npoints, scalars, \
        name, input_vtk = read_vtk(vtk_file, return_first=True, return_array=True)
    indices = range(1, npoints + 1)

    # Load label lists
    labels_to_replace, new_labels = read_columns(relabel_list, 2)
    for i, new_label in enumerate(new_labels):

        # Find which vertices have the label
        indices = np.where(scalars == int(labels_to_replace[i]))[0]
        scalars[indices] = int(new_label)

    relabeled_vtk = os.path.join(os.getcwd(),
                                 os.path.basename(vtk_file).split('.')[0] + \
                                 '.' + new_string)
    write_vtk(relabeled_vtk, points, indices, lines, faces,
              [scalars.tolist()], ['Labels'])

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
    hemi : string
        hemisphere ['lh' or 'rh']
    subject : string
        subject name
    annot_name : string
        name of .annot file (without pre- or post-pends)
    relabel_file : string
        text file with old and new RGB values
    new_annot_name : string
        new .annot name

    """
    from nipype.interfaces.base import CommandLine

    cli = CommandLine(command='mris_translate_annotation')
    cli.inputs.args = ' '.join([subject, hemi, annot_name, relabel_file, new_annot_name])
    cli.cmdline
    cli.run()

    return new_annot_name

def overwrite_volume_labels(source, target, output_file='', ignore_labels=[0]):
    """
    Overwrite target labels with source labels (same volume dimensions).

    Parameters
    ----------
    source : string
        labeled nibabel-readable (e.g., nifti) file
    target : string
        labeled nibabel-readable (e.g., nifti) file
    output_file : string
        labeled nibabel-readable (e.g., nifti) file
    ignore_labels : list
        list of source labels to ignore

    Returns
    -------
    output_file : string
        labeled nibabel-readable (e.g., nifti) file

    Examples
    --------
    >>> # Overwrite DKT25 with DKT31 labels
    >>> import os
    >>> from mindboggle.labels.relabel import overwrite_volume_labels
    >>> from mindboggle.utils.plots import plot_volumes
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> source = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> target = os.path.join(data_path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> output_file = ''
    >>> ignore_labels = [0]
    >>> output_file = overwrite_volume_labels(source, target, output_file, ignore_labels)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   os.path.basename(source).split('.')[0] +
                                   '_and_' +
                                   os.path.basename(target).split('.')[0] +
                                   '.nii.gz')

    # Load labeled image volumes:
    vol_source = nb.load(source)
    vol_target = nb.load(target)
    xfm = vol_target.get_affine()
    data_source = vol_source.get_data().ravel()
    data_target = vol_target.get_data().ravel()

    # Initialize output:
    new_data = data_target.copy()

    # Overwrite target labels with source labels:
    IX = [(i,x) for i,x in enumerate(data_source) if x not in ignore_labels]
    I = [x[0] for x in IX]
    X = [x[1] for x in IX]
    new_data[I] = X

    # Reshape to original dimensions:
    new_data = np.reshape(new_data, vol_target.shape)

    # Save relabeled file:
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    return output_file

