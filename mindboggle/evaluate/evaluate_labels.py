#!/usr/bin/python

"""
Compute surface and volume label overlaps.

Compute the Dice and Jaccard overlap measures for each labeled region
of two labeled surfaces or image volumes, for example one that has been
manually labeled and one that has been automatically labeled.

For surface overlap, this program simply calls Joachim Giard's code.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def measure_surface_overlap(command, labels_file1, labels_file2):
    """
    Measure surface overlap using Joachim Giard's code.

    Parameters
    ----------
    command : surface overlap C++ executable command
    labels_file1 : ``vtk file`` with index labels for scalar values
    labels_file2 : ``vtk file`` with index labels for scalar values

    Returns
    -------
    overlap_file : string
        name of output text file with overlap results

    Examples
    --------
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import measure_surface_overlap
    >>> ccode_path = os.environ['MINDBOGGLE_TOOLS']
    >>> command = os.path.join(ccode_path, 'surface_overlap', 'SurfaceOverlapMain')
    >>> data_path = os.path.join(os.environ['MINDBOGGLE_DATA'], 'rescan_labels')
    >>> appnd = '.lh.pial.labels.DKT31.manual.vtk'
    >>> # Two misaligned label files:
    >>> file1 = os.path.join(data_path, 'rescan_labels', 'MMRR-21-2' + appnd)
    >>> file2 = os.path.join(data_path, 'rescan_labels', 'MMRR-21-2_rescan' + appnd)
    >>> measure_surface_overlap(command, file1, file2)

    """
    import os
    from nipype.interfaces.base import CommandLine

    overlap_filename = os.path.basename(labels_file1) + '_and_' + \
                       os.path.basename(labels_file2) + '.txt'
    overlap_file = os.path.join(os.getcwd(), overlap_filename)
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([labels_file1, labels_file2, overlap_file])
    cli.cmdline
    cli.run()

    return overlap_file


def measure_volume_overlap(labels, atlas_file, input_file):
    """
    Measure overlap between individual label regions
    in source and target nifti (nii.gz) images.

    Parameters
    ----------
    labels : list of label indices
    atlas_file : source image, consisting of index-labeled pixels/voxels
    input_file : target image, consisting of index-labeled pixels/voxels

    Returns
    -------
    overlaps : numpy array
        overlap values
    out_file : string
        output text file name with overlap values

    """
    import os
    import numpy as np
    import nibabel as nb

    save_output = True

    # Load labeled image volumes
    input_data = nb.load(input_file).get_data().ravel()
    atlas_data = nb.load(atlas_file).get_data().ravel()
    #print(input_file + ' ' + atlas_file)
    #print(np.unique(input_data))
    #print(np.unique(atlas_data))

    # Initialize output
    overlaps = np.zeros((len(labels), 3))

    # Loop through labels
    for ilabel, label in enumerate(labels):
        label = int(label)
        overlaps[ilabel, 0] = label

        # Find which voxels contain the label in each volume
        input_indices = np.where(input_data==label)[0]
        atlas_indices = np.where(atlas_data==label)[0]
        input_label_sum = len(input_indices)
        atlas_label_sum = len(atlas_indices)

        # Determine their intersection and union
        intersect_label_sum = len(np.intersect1d(input_indices, atlas_indices))
        union_label_sum = len(np.union1d(input_indices, atlas_indices))
        #print('{0} {1} {2} {3} {4}'.format(label, input_label_sum, atlas_label_sum,
        #                              intersect_label_sum, union_label_sum))

        # There must be at least one voxel with the label in each volume
        if input_label_sum * atlas_label_sum > 0:

            # Compute Dice and Jaccard coefficients
            dice = np.float(2.0 * intersect_label_sum) / +\
                   (input_label_sum + atlas_label_sum)
            jacc = np.float(intersect_label_sum) / union_label_sum
            overlaps[ilabel, 1:3] = [dice, jacc]
            print('label: {0}, dice: {1:.2f}, jacc: {2:.2f}'.format(
                  label, dice, jacc))

        # NOTE:  untested:
        if save_output:
            input_name = os.path.splitext(os.path.basename(input_file))[0]
            atlas_name = os.path.splitext(os.path.basename(atlas_file))[0]
            out_file = os.path.join(os.getcwd(), 'labelvolume_dice_jacc_' +
                                    input_name + '_vs_' + atlas_name + '.txt')
            np.savetxt(out_file, overlaps, fmt='%d %.4f %.4f',
                       delimiter='\t', newline='\n')

    return overlaps, out_file

