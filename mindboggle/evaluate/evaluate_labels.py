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

    """
#    import os
    from nipype.interfaces.base import CommandLine

    #import numpy as np
    #from utils import io_vtk
    #Points, Faces, Scalars1, n_vertices = io_vtk.load_scalar(labels_file1, return_arrays=1)
    #Points, Faces, Scalars2, n_vertices = io_vtk.load_scalar(labels_file2, return_arrays=1)
    #print(np.unique(Scalars1))
    #print(np.unique(Scalars2))

    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([labels_file1, labels_file2])
    cli.cmdline

#    name1 = os.path.splitext(os.path.basename(labels_file1))[0]
#    name2 = os.path.splitext(os.path.basename(labels_file2))[0]
#    out_file = os.path.join(os.getcwd(), 'labelvolume_dice_jacc_for_' +
#                            name1 + '_vs_' + name2 + '.txt')

    return cli.run()


def measure_volume_overlap(labels, atlas_file, input_file):
    """
    Measure overlap between individual label regions
    in source and target nifti (nii.gz) images.

    Parameters
    ----------
    labels : list of label indices
    atlas_file : source image, consisting of index-labeled pixels/voxels
    input_file : target image, consisting of index-labeled pixels/voxels

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

