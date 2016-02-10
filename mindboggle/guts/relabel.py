#!/usr/bin/env python
"""
Relabel surface or volume or annot files.


Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
    >>> from mindboggle.guts.relabel import relabel_volume
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_file = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(input_file, input_file + '.nii.gz')
    >>> input_file = input_file + '.nii.gz'
    >>> dkt = DKTprotocol()
    >>> old_labels = dkt.cerebrum_cortex_numbers + dkt.cerebrum_noncortex_numbers
    >>> ctx = [5000 for x in dkt.cerebrum_cortex_numbers]
    >>> nonctx = [6000 for x in dkt.cerebrum_noncortex_numbers]
    >>> new_labels = ctx + nonctx
    >>> output_file = ''
    >>> output_file = relabel_volume(input_file, old_labels, new_labels,
    ...                              output_file)

    View nifti file (skip test):

    >>> from mindboggle.mio.plots import plot_volumes # doctest: +SKIP
    >>> plot_volumes(output_file) # doctest: +SKIP

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
        raise IOError("relabel_volume() did not create " + output_file + ".")

    return output_file


def remove_volume_labels(input_file, labels_to_remove, output_file='',
                         second_file=''):
    """
    Remove labels from an image volume
    (or corresponding voxels in a 2nd volume).

    Parameters
    ----------
    input_file : string
        labeled nibabel-readable (e.g., nifti) file
    labels_to_remove : list of integers
        labels to remove
    output_file : string
        output file name
    second_file : string
        second nibabel-readable file (erase voxels in this file instead)

    Returns
    -------
    output_file : string
        output file name

    Examples
    --------
    >>> # Remove subcortical labels
    >>> import os
    >>> from mindboggle.guts.relabel import remove_volume_labels
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_file = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(input_file, input_file + '.nii.gz')
    >>> input_file = input_file + '.nii.gz'
    >>> second_file = ''
    >>> labels_to_remove = range(1,300) # Remove noncortex (+aseg) labels
    >>> labels_to_remove.extend([1000,1001,2000,2001])
    >>> labels_to_remove.extend(range(2000,2036)) # Remove right cortex labels
    >>> output_file = ''
    >>> output_file = remove_volume_labels(input_file, labels_to_remove,
    ...                                    output_file, second_file)

    View nifti file (skip test):

    >>> from mindboggle.mio.plots import plot_volumes # doctest: +SKIP
    >>> plot_volumes(output_file) # doctest: +SKIP

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
        raise IOError("remove_volume_labels() did not create " + output_file
                      + ".")

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
    >>> from mindboggle.guts.relabel import keep_volume_labels
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_file = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(input_file, input_file + '.nii.gz')
    >>> input_file = input_file + '.nii.gz'
    >>> second_file = ''
    >>> labels_to_keep = range(1000, 1036)
    >>> output_file = ''
    >>> output_file = keep_volume_labels(input_file, labels_to_keep,
    ...                                  output_file, second_file)

    View nifti file (skip test):

    >>> from mindboggle.mio.plots import plot_volumes # doctest: +SKIP
    >>> plot_volumes(output_file) # doctest: +SKIP

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
        raise IOError("keep_volume_labels() did not create " + output_file + ".")

    return output_file


def relabel_surface(vtk_file, hemi='', old_labels=[], new_labels=[],
                    erase_remaining=True, erase_labels=[], erase_value=-1,
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
    >>> import numpy as np
    >>> from mindboggle.guts.relabel import relabel_surface
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> hemi = 'lh'
    >>> old_labels = [1003,1009,1030]
    >>> new_labels = [0,500,1000]
    >>> erase_remaining = True
    >>> erase_labels = [0]
    >>> erase_value = -1
    >>> output_file = ''
    >>> output_file = relabel_surface(vtk_file, hemi, old_labels, new_labels,
    ...     erase_remaining, erase_labels, erase_value, output_file)
    >>> labels, name = read_scalars(output_file, True, True)
    >>> np.unique(labels)
    array([  -1, 1000, 1500, 2000])

    View relabeled surface file (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces(output_file) # doctest: +SKIP

    """
    import os
    import numpy as np
    from mindboggle.mio.vtks import read_vtk, write_vtk

    # Load labeled vtk surfaces:
    points, indices, lines, faces, scalars, scalar_names, npoints, \
        input_vtk = read_vtk(vtk_file, return_first=True, return_array=True)
    new_scalars = scalars[:]

    # Raise an error if inputs set incorrectly:
    if (new_labels and not old_labels) or \
       (hemi and hemi not in ['lh','rh']) or \
       (new_labels and len(old_labels) != len(new_labels)) or \
       (erase_remaining and not old_labels):
        raise IOError("Please check inputs for relabel_surface().")

    # Loop through unique labels in scalars:
    ulabels = np.unique(scalars)
    for label in ulabels:
        I = np.where(scalars == label)[0]

        # If label in erase_labels list, replace with erase_value:
        if label in erase_labels:
            new_scalars[I] = erase_value

        # If label in old_labels list, replace with corresponding new label,
        # and if hemi set, add 1000 or 2000 to the new label:
        elif label in old_labels and (len(old_labels) == len(new_labels)):
            new_label = new_labels[old_labels.index(label)]
            if hemi == 'lh':
                new_scalars[I] = 1000 + new_label
            elif hemi == 'rh':
                new_scalars[I] = 2000 + new_label
            else:
                new_scalars[I] = new_label

        # If labels not set then optionally add hemi value:
        elif hemi and not new_labels:
            if hemi == 'lh':
                new_scalars[I] = 1000 + label
            elif hemi == 'rh':
                new_scalars[I] = 2000 + label

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
        raise IOError("relabel_surface() did not create " + output_file + ".")

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
    >>> from mindboggle.guts.relabel import overwrite_volume_labels
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> source = fetch_data(urls['freesurfer_labels'])
    >>> target = fetch_data(urls['ants_labels'])
    >>> os.rename(source, source + '.nii.gz')
    >>> source = source + '.nii.gz'
    >>> os.rename(target, target + '.nii.gz')
    >>> target = target + '.nii.gz'
    >>> output_file = ''
    >>> ignore_labels = [0]
    >>> erase_labels = False
    >>> background_value = -1
    >>> output_file = overwrite_volume_labels(source, target, output_file,
    ...     ignore_labels, erase_labels, background_value)

    View nifti file (skip test):

    >>> from mindboggle.mio.plots import plot_volumes # doctest: +SKIP
    >>> plot_volumes(output_file) # doctest: +SKIP

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
        raise IOError('{0} and {1} need to be the same shape.'.
                      format(source, target))
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
        raise IOError("overwrite_volume_labels() did not create {0}."
                      .format(output_file))

    return output_file


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()