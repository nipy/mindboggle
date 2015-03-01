#!/usr/bin/env python
"""
Relabel surface or volume or annot files.


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

<<<<<<< HEAD

def relabel_volume(input_file, old_labels, new_labels, output_file=''):
=======
def relabel_volume(input_file, old_labels, new_labels):
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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
<<<<<<< HEAD
    output_file : string
        output file name

    Returns
    -------
    output_file : string
        output file name
=======
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    Examples
    --------
    >>> # Convert DKT31 to DKT25 labels
<<<<<<< HEAD
    >>> from mindboggle.labels.relabel import relabel_volume
    >>> from mindboggle.LABELS import DKTprotocol
    >>> from mindboggle.utils.plots import plot_volumes
    >>> # Convert DKT31 to DKT25 protocol:
    >>> old_labels = [1010,1023,1026,1027,1019,1020,2010,2023,2026,2027,2019,2020]
    >>> new_labels = [1002,1002,1002,1003,1018,1018,2002,2002,2002,2003,2018,2018]
    >>> # Convert labels to non/cortex segmentation:
    >>> input_file = 'labels.DKT31.manual.MNI152.nii.gz'
    >>> #dkt = DKTprotocol()
    >>> #old_labels = dkt.cerebrum_cortex_numbers + dkt.cerebrum_noncortex_numbers
    >>> #new_labels = dkt.cerebrum_cortex_numbers + [0 for x in dkt.cerebrum_noncortex_numbers]
    >>> #new_labels = [2 for x in dkt.cerebrum_cortex_numbers] + [3 for x in dkt.cerebrum_noncortex_numbers]
    >>> output_file = 'labels.DKT25.manual.MNI152.nii.gz'
    >>> output_file = relabel_volume(input_file, old_labels, new_labels, output_file)
    >>> # View
    >>> plot_volumes(output_file)
=======
    >>> import os
    >>> from mindboggle.utils.io_file import read_columns
    >>> from mindboggle.labels.relabel import relabel_volume
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> relabel_file = os.path.join(data_path, 'info', 'labels.volume.DKT31to25.txt')
    >>> old_labels, new_labels = read_columns(relabel_file, 2)
    >>> relabel_volume(input_file, old_labels, new_labels)
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

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
<<<<<<< HEAD
    if not output_file:
        output_file = os.path.join(os.getcwd(), os.path.basename(input_file))
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        s = "relabel_volume() did not create " + output_file + "."
        raise(IOError(s))

    return output_file


def remove_volume_labels(input_file, labels_to_remove, output_file='',
                         second_file=''):
    """
    Remove labels from an image volume
    (or corresponding voxels in a 2nd volume).
=======
    output_file = os.path.join(os.getcwd(), os.path.basename(input_file))
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    return output_file

def remove_volume_labels(input_file, labels_to_remove):
    """
    Remove labels from an image volume.
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    Parameters
    ----------
    input_file : string
        labeled nibabel-readable (e.g., nifti) file
    labels_to_remove : list of integers
        labels to remove
<<<<<<< HEAD
    output_file : string
        output file name
    second_file : string
        second nibabel-readable file (erase voxels in this file instead)

    Returns
    -------
    output_file : string
        output file name
=======
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    Examples
    --------
    >>> # Remove subcortical labels
    >>> import os
    >>> from mindboggle.labels.relabel import remove_volume_labels
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
<<<<<<< HEAD
    >>> second_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> labels_to_remove = range(1,300) # Remove noncortical (+aseg) labels
    >>> labels_to_remove.extend([1000,1001,2000,2001])
    >>> labels_to_remove.extend(range(2000,2036)) # Remove right cortical labels
    >>> output_file = ''
    >>> remove_volume_labels(input_file, labels_to_remove, output_file, second_file)
=======
    >>> labels_to_remove = range(1,300) # Remove noncortical (+aseg) labels
    >>> labels_to_remove.extend([1000,1001,2000,2001])
    >>> remove_volume_labels(input_file, labels_to_remove)
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    """
    import os
    import numpy as np
    import nibabel as nb

<<<<<<< HEAD
    #-------------------------------------------------------------------------
    # Load labeled image volume and extract data as 1-D array:
    #-------------------------------------------------------------------------
=======
    # Load labeled image volume and extract data as 1-D array
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    vol = nb.load(input_file)
    xfm = vol.get_affine()
    data = vol.get_data().ravel()

<<<<<<< HEAD
    #-------------------------------------------------------------------------
    # If second file specified, erase voxels whose corresponding
    # voxels in the input_file have labels in labels_to_remove:
    #-------------------------------------------------------------------------
    if second_file:
        # Load second image volume and extract data as 1-D array:
        vol = nb.load(second_file)
        xfm = vol.get_affine()
        new_data = vol.get_data().ravel()
        if not output_file:
            output_file = os.path.join(os.getcwd(),
                                       os.path.basename(second_file))
    #-------------------------------------------------------------------------
    # If second file not specified, remove labels in labels_to_remove:
    #-------------------------------------------------------------------------
    else:
        new_data = data.copy()
        if not output_file:
            output_file = os.path.join(os.getcwd(),
                                       os.path.basename(input_file))

    #-------------------------------------------------------------------------
    # Erase voxels as specified above:
    #-------------------------------------------------------------------------
    ulabels = np.unique(data)
    for label in ulabels:
        label = int(label)
        if label in labels_to_remove:
            new_data[np.where(data == label)[0]] = 0

    #-------------------------------------------------------------------------
    # Reshape to original dimensions:
    #-------------------------------------------------------------------------
    new_data = np.reshape(new_data, vol.shape)

    #-------------------------------------------------------------------------
    # Save relabeled file:
    #-------------------------------------------------------------------------
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        s = "remove_volume_labels() did not create " + output_file + "."
        raise(IOError(s))

    return output_file


def keep_volume_labels(input_file, labels_to_keep, output_file='',
                       second_file=''):
    """
    Keep only given labels in an image volume (or use to mask second volume).

    Parameters
    ----------
    input_file : string
        labeled nibabel-readable (e.g., nifti) file
    labels_to_keep : list of integers
        labels to keep
    output_file : string
        output file name
    second_file : string
        second nibabel-readable file (keep/erase voxels in this file instead)

    Returns
    -------
    output_file : string
        output file name

    Examples
    --------
    >>> # Remove right hemisphere labels
    >>> import os
    >>> from mindboggle.labels.relabel import keep_volume_labels
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> second_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> labels_to_keep = range(1000, 1036)
    >>> output_file = ''
    >>> keep_volume_labels(input_file, labels_to_keep, output_file, second_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    #-------------------------------------------------------------------------
    # Load labeled image volume and extract data as 1-D array:
    #-------------------------------------------------------------------------
    vol = nb.load(input_file)
    xfm = vol.get_affine()
    data = vol.get_data().ravel()

    #-------------------------------------------------------------------------
    # If second file specified, erase voxels whose corresponding
    # voxels in the input_file have labels not in labels_to_keep:
    #-------------------------------------------------------------------------
    if second_file:
        # Load second image volume and extract data as 1-D array:
        vol = nb.load(second_file)
        xfm = vol.get_affine()
        new_data = vol.get_data().ravel()
        if not output_file:
            output_file = os.path.join(os.getcwd(),
                                       os.path.basename(second_file))
    #-------------------------------------------------------------------------
    # If second file not specified, remove labels not in labels_to_keep:
    #-------------------------------------------------------------------------
    else:
        new_data = data.copy()
        if not output_file:
            output_file = os.path.join(os.getcwd(),
                                       os.path.basename(input_file))

    #-------------------------------------------------------------------------
    # Erase voxels as specified above:
    #-------------------------------------------------------------------------
    ulabels = np.unique(data)
    for label in ulabels:
        label = int(label)
        if label not in labels_to_keep:
            new_data[np.where(data == label)[0]] = 0

    #-------------------------------------------------------------------------
    # Reshape to original dimensions:
    #-------------------------------------------------------------------------
    new_data = np.reshape(new_data, vol.shape)

    #-------------------------------------------------------------------------
    # Save relabeled file:
    #-------------------------------------------------------------------------
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        s = "keep_volume_labels() did not create " + output_file + "."
        raise(IOError(s))

    return output_file


def relabel_surface(vtk_file, old_labels=[], new_labels=[],
                    erase_remaining=True, erase_labels=[], erase_value=-1,
                    output_file=''):
=======
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    """
    Relabel surface in a VTK file.

    Parameters
    ----------
    vtk_file : string
        input labeled VTK file
<<<<<<< HEAD
    old_labels : list of integers
        old labels (empty list if labels drawn from vtk scalars)
    new_labels : list of integers
        new labels (empty list if labels drawn from vtk scalars)
    erase_remaining : Boolean
        set all values not in old_labels to erase_value?
    erase_labels : list of integers
        values to erase (set to erase_value)
    erase_value : integer
        set vertices with labels in erase_labels to this value
    output_file : string
        new vtk file name

    Returns
    -------
    output_file : string
        new vtk file name

    Examples
    --------
    >>> from mindboggle.labels.relabel import relabel_surface
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> vtk_file = 'rh.labels.DKT31.manual.vtk'
    >>> old_labels = [10, 23, 26, 19, 20, 27]
    >>> new_labels = [2, 2, 2, 18, 18, 3]
    >>> erase_remaining = False
    >>> erase_labels = [0]
    >>> erase_value = -1
    >>> output_file = ''
    >>> #
    >>> relabel_surface(vtk_file, old_labels, new_labels, erase_remaining, erase_labels, erase_value, output_file)
    >>> # View
    >>> plot_surfaces('relabeled_rh.labels.DKT31.manual.vtk')
=======
    relabel_list : string
        text file with two columns of label numbers --
        all regions receive the 2nd label per row.
    new_string : string
        new ending of vtk_file name (e.g., 'labels.DKT25.vtk')
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, write_vtk
<<<<<<< HEAD

    # Load labeled vtk surfaces:
    faces, lines, indices, points, npoints, scalars, \
        name, input_vtk = read_vtk(vtk_file, return_first=True, return_array=True)
    new_scalars = scalars[:]

    # Raise an error if inputs set incorrectly:
    if (new_labels and not old_labels) or \
            (len(old_labels) != len(new_labels)) or \
            (erase_remaining and not old_labels):
        raise IOError("Please check inputs for relabel_surface().")

    # Loop through unique labels in scalars:
    ulabels = np.unique(scalars)
    for label in ulabels:
        I = np.where(scalars == label)[0]

        # If label in erase_labels list, replace with erase_value:
        if label in erase_labels:
            new_scalars[I] = erase_value

        # If label in old_labels list, replace with corresponding new label:
        elif label in old_labels:
            new_scalars[I] = new_labels[old_labels.index(label)]

        # If label unaccounted for and erase_remaining, set to erase_value:
        elif erase_remaining:
            new_scalars[I] = erase_value

    # Ensure that the new scalars are integer values:
    new_scalars = [int(x) for x in new_scalars]

    # Write output VTK file:
    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'relabeled_' + os.path.basename(vtk_file))
    write_vtk(output_file, points, indices, lines, faces,
              [new_scalars], ['Labels'], scalar_type='int')
    if not os.path.exists(output_file):
        s = "relabel_surface() did not create " + output_file + "."
        raise(IOError(s))

    return output_file


def overwrite_volume_labels(source, target, output_file='', ignore_labels=[0],
                            erase_labels=True, background_value=-1):
    """
    For every label in a source image, optionally erase all voxels in the
    target image with this label (if erase_labels is True), and
    for every voxel in the source image with this label,
    assign the label to the corresponding voxel in the target image.
    The source and target images must have the same volume dimensions.

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
    erase_labels : Boolean
        erase target labels (that are in source) before overwriting?
    background_value : integer
        background value (if erase_labels==True)

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
    >>> erase_labels = True
    >>> background_value = -1
    >>> output_file = overwrite_volume_labels(source, target, output_file, ignore_labels, erase_labels, background_value)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    if not output_file:
        output_file = os.path.join(os.getcwd(), os.path.basename(source) +
                                   '_to_' + os.path.basename(target))
    # Load labeled image volumes:
    vol_source = nb.load(source)
    vol_target = nb.load(target)
    if vol_source.shape != vol_target.shape:
        raise(IOError('{0} and {1} need to be the same shape.'.
                      format(source, target)))
    xfm = vol_target.get_affine()
    data_source = vol_source.get_data().ravel()
    data_target = vol_target.get_data().ravel()

    # Initialize output:
    new_data = data_target.copy()

    # Find indices with labels in source:
    IX = [(i,x) for i,x in enumerate(data_source) if x not in ignore_labels]
    I = [x[0] for x in IX]
    X = [x[1] for x in IX]

    # Erase target labels (that are in source) before overwriting:
    if erase_labels:
        rm_labels = np.unique(X)
        Irm = [i for i,x in enumerate(data_target) if x in rm_labels]
        new_data[Irm] = background_value

    # Overwrite target labels with source labels:
    new_data[I] = X

    # Reshape to original dimensions:
    new_data = np.reshape(new_data, vol_target.shape)

    # Save relabeled file:
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        s = "overwrite_volume_labels() did not create " + output_file + "."
        raise(IOError(s))

    return output_file

=======
    from mindboggle.utils.io_file import read_columns

    # Load labeled vtk surfaces
    faces, lines, indices, points, npoints, scalars, \
        name = read_vtk(vtk_file, return_first=True, return_array=True)
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
