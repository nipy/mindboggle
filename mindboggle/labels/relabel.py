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
    >>> from mindboggle.labels.relabel import relabel_volume
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> old_labels = [1010,1023,1026,1027,1019,1020,2010,2023,2026,2027,2019,2020]
    >>> new_labels = [1002,1002,1002,1003,1018,1018,2002,2002,2002,2003,2018,2018]
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

def relabel_surface(vtk_file, hemi='', old_labels=[], new_labels=[],
                    output_file=''):
    """
    Relabel surface in a VTK file.

    Parameters
    ----------
    vtk_file : string
        input labeled VTK file
    hemi : string
        hemisphere ('lh' or 'rh' or '')
    old_labels : list of integers
        old labels
    new_labels : list of integers
        new labels
    output_file : string
        new vtk file name

    Returns
    -------
    output_file : string
        new vtk file name

    Examples
    --------
    >>> import os
    >>> from mindboggle.labels.relabel import relabel_surface
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> hemi = 'lh'
    >>> old_labels = []
    >>> new_labels = []
    >>> output_file = ''
    >>> #
    >>> relabel_surface(vtk_file, hemi, old_labels, new_labels, output_file)
    >>> # View
    >>> plot_vtk('relabeled_lh.labels.DKT25.manual.vtk')

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, write_vtk

    # Load labeled vtk surfaces:
    faces, lines, indices, points, npoints, scalars, \
        name, input_vtk = read_vtk(vtk_file, return_first=True, return_array=True)

    # Add a hemisphere value to each label:
    if hemi:
        ulabels = np.unique(scalars)
        for label in ulabels:
            I = np.where(scalars == int(label))[0]
            if hemi == 'lh':
                scalars[I] = 1000 + int(label)
            elif hemi == 'rh':
                scalars[I] = 2000 + int(label)
    # OR replace each old label with a corresponding new label:
    else:
        for ilabel, new_label in enumerate(new_labels):
            I = np.where(scalars == int(old_labels[ilabel]))[0]
            scalars[I] = int(new_label)

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'relabeled_' + os.path.basename(vtk_file))
    write_vtk(output_file, points, indices, lines, faces,
              [scalars.tolist()], ['Labels'])

    return output_file

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
