#!/usr/bin/python

"""
Compute volume label overlaps.

Compute the Dice and Jaccard overlap measures for each labeled region
of two labeled image volumes, for example one that has been manually labeled
(atlas_file) and one that has been automatically labeled (input_file).
Perform this computation for each label in a list.


Author:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""


def measure_volume_overlap(labels, atlas_file, input_file):
    """
    Measure overlap between individual label regions
    in source and target nifti (nii.gz) images.

    Input:
    arg1:  list of label indices
    arg2, arg3:  source, target images, consisting of index-labeled pixels/voxels

    """

    import os
    import numpy as np
    import nibabel as nb

    save_output = True

    # Load labeled image volumes
    input_data = nb.load(input_file).get_data().ravel()
    atlas_data = nb.load(atlas_file).get_data().ravel()
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
        print('{} {} {} {} {}'.format(label, input_label_sum, atlas_label_sum,
                                      intersect_label_sum, union_label_sum))

        # There must be at least one voxel with the label in each volume
        if input_label_sum * atlas_label_sum > 0:

            # Compute Dice and Jaccard coefficients
            dice = np.float(2.0 * intersect_label_sum) / (input_label_sum + atlas_label_sum)
            jacc = np.float(intersect_label_sum) / union_label_sum
            overlaps[ilabel, 1:3] = [dice, jacc]
            print('label: {}, dice: {:.2f}, jacc: {:.2f}'.format(label, dice, jacc))

        # NOTE:  untested:
        if save_output:
            input_name = os.path.splitext(os.path.basename(input_file))[0]
            atlas_name = os.path.splitext(os.path.basename(atlas_file))[0]
            out_file = os.path.join(os.getcwd(), input_name + '_vs_' + atlas_name + '.txt')
            np.savetxt(out_file, overlaps, fmt='%d %.4f %.4f', delimiter='\t', newline='\n')

    return overlaps, out_file
