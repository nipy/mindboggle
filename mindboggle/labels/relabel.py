#!/usr/bin/env python
"""
Relabel surface or volume or annot files.


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def relabel_volume(input_file, old_labels, new_labels, output_file=''):
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
    output_file : string
        output file name

    Returns
    -------
    output_file : string
        output file name

    Examples
    --------
    >>> # Convert DKT31 to DKT25 labels
    >>> import os
    >>> from mindboggle.labels.relabel import relabel_volume
    >>> from mindboggle.LABELS import DKTprotocol
    >>> from mindboggle.utils.plots import plot_volumes
    >>> # Convert DKT31 to DKT25 protocol:
    >>> #data_path = os.environ['MINDBOGGLE_DATA']
    >>> #input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> #old_labels = [1010,1023,1026,1027,1019,1020,2010,2023,2026,2027,2019,2020]
    >>> #new_labels = [1002,1002,1002,1003,1018,1018,2002,2002,2002,2003,2018,2018]
    >>> # Convert labels to non/cortex segmentation:
    >>> input_file = os.path.join(os.environ['HOME'], 'mindboggled/OASIS-TRT-20-1/labels/FreeSurfer_labels.nii.gz')
    >>> dkt = DKTprotocol()
    >>> old_labels = dkt.cortex_numbers + dkt.noncortex_numbers
    >>> new_labels = [2 for x in dkt.cortex_numbers] + [3 for x in dkt.noncortex_numbers]
    >>> output_file = ''
    >>> output_file = relabel_volume(input_file, old_labels, new_labels, output_file)
    >>> # View
    >>> plot_volumes(output_file)

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
    if not output_file:
        output_file = os.path.join(os.getcwd(), os.path.basename(input_file))
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file


def remove_volume_labels(input_file, labels_to_remove, output_file=''):
    """
    Remove labels from an image volume.

    Parameters
    ----------
    input_file : string
        labeled nibabel-readable (e.g., nifti) file
    labels_to_remove : list of integers
        labels to remove
    output_file : string
        output file name

    Returns
    -------
    output_file : string
        output file name

    Examples
    --------
    >>> # Remove subcortical labels
    >>> import os
    >>> from mindboggle.labels.relabel import remove_volume_labels
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> input_file = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> labels_to_remove = range(1,300) # Remove noncortical (+aseg) labels
    >>> labels_to_remove.extend([1000,1001,2000,2001])
    >>> output_file = ''
    >>> remove_volume_labels(input_file, labels_to_remove, output_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load labeled image volume and extract data as 1-D array:
    vol = nb.load(input_file)
    xfm = vol.get_affine()
    data = vol.get_data().ravel()

    # Initialize output:
    new_data = data.copy()

    # Loop through labels:
    for label in labels_to_remove:
        label = int(label)
        # Relabel:
        new_data[np.where(data==label)[0]] = 0

    # Reshape to original dimensions:
    new_data = np.reshape(new_data, vol.shape)

    # Save relabeled file:
    if not output_file:
        output_file = os.path.join(os.getcwd(), os.path.basename(input_file))
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file


def keep_volume_labels(input_file, labels_to_keep, output_file=''):
    """
    Keep only given labels in an image volume.

    Parameters
    ----------
    input_file : string
        labeled nibabel-readable (e.g., nifti) file
    labels_to_keep : list of integers
        labels to keep
    output_file : string
        output file name

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
    >>> labels_to_keep = range(1000, 1036)
    >>> output_file = ''
    >>> keep_volume_labels(input_file, labels_to_keep, output_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load labeled image volume and extract data as 1-D array:
    vol = nb.load(input_file)
    xfm = vol.get_affine()
    data = vol.get_data().ravel()

    # Initialize output:
    new_data = data.copy()

    # Loop through unique labels:
    ulabels = np.unique(data)
    for label in ulabels:
        label = int(label)
        if label not in labels_to_keep:
            # Relabel if not a label to keep:
            new_data[np.where(data == label)[0]] = 0

    # Reshape to original dimensions:
    new_data = np.reshape(new_data, vol.shape)

    # Save relabeled file:
    if not output_file:
        output_file = os.path.join(os.getcwd(), os.path.basename(input_file))
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

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
        if set, add 1000 to left and 2000 to right hemisphere labels;
    old_labels : list of integers
        old labels (empty list if labels drawn from vtk scalars);
        may be used in conjunction with hemi
    new_labels : list of integers
        new labels (empty list if labels drawn from vtk scalars);
        may be used in conjunction with hemi
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
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> hemi = 'lh'
    >>> old_labels = []
    >>> new_labels = []
    >>> output_file = ''
    >>> #
    >>> relabel_surface(vtk_file, hemi, old_labels, new_labels, output_file)
    >>> # View
    >>> plot_surfaces('relabeled_lh.labels.DKT25.manual.vtk')

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, write_vtk

    # Load labeled vtk surfaces:
    faces, lines, indices, points, npoints, scalars, \
        name, input_vtk = read_vtk(vtk_file, return_first=True, return_array=True)

    # Add a hemisphere value to each unique label drawn from scalars:
    if hemi and not old_labels and not new_labels:
        ulabels = np.unique(scalars)
        for label in ulabels:
            I = np.where(scalars == int(label))[0]
            if hemi == 'lh':
                scalars[I] = 1000 + int(label)
            elif hemi == 'rh':
                scalars[I] = 2000 + int(label)
    # OR replace each old label with a corresponding new label
    # (hemisphere setting optionally adds 1000 or 2000 to the new label):
    else:
        for ilabel, new_label in enumerate(new_labels):
            I = np.where(scalars == int(old_labels[ilabel]))[0]
            if hemi == 'lh':
                scalars[I] = 1000 + int(new_label)
            elif hemi == 'rh':
                scalars[I] = 2000 + int(new_label)
            else:
                scalars[I] = int(new_label)

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'relabeled_' + os.path.basename(vtk_file))
    write_vtk(output_file, points, indices, lines, faces,
              [scalars.tolist()], ['Labels'])

    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file


def overwrite_volume_labels(source, target, output_file='', ignore_labels=[0],
                            erase_labels=True, background_value=-1):
    """
    For every label in a source image, optionally erase all voxels in the
    target image with this label (if erase_labels is True), and
    for every voxel in the source image with this label,
    assign the label to the corresponding voxel in the target image.

    Note::

        Assumes same volume dimensions.

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
        raise(IOError(output_file + " not found"))

    return output_file

