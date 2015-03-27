#!/usr/bin/env python

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

def measure_surface_overlap(command, labels_file1, labels_file2, overlap_file):
    """
    Measure surface overlap using Joachim Giard's code.

    Parameters
    ----------
    command : string
        surface overlap C++ executable command
    labels_file1 : string
        ``vtk file`` with index labels for scalar values
    labels_file2 : string
        ``vtk file`` with index labels for scalar values
    overlap_file : string
        (optional) output file name

    Returns
    -------
    overlap_file : string
        name of output text file with overlap results

    Examples
    --------
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import measure_surface_overlap
    >>> from mindboggle.mindboggle import hashes_url
    >>> from mindboggle.mio.fetch_data import fetch_check_data
    >>> hashes, url, cache_env, cache = hashes_url()
    >>> ccode_path = os.environ['MINDBOGGLE_TOOLS']
    >>> command = os.path.join(ccode_path, 'surface_overlap', 'SurfaceOverlapMain')
    >>> label_file1 = 'lh.labels.DKT25.manual.vtk'
    >>> label_file2 = 'lh.labels.DKT31.manual.vtk'
    >>> file1 = fetch_check_data(label_file1, url, hashes, cache_env, cache)
    >>> file2 = fetch_check_data(label_file2, url, hashes, cache_env, cache)
    >>> overlap_file = ''
    >>> measure_surface_overlap(command, file1, file2, overlap_file)

    """
    import os
    from nipype.interfaces.base import CommandLine

    if not overlap_file:
        overlap_file = os.path.basename(labels_file1) + '_and_' + \
                           os.path.basename(labels_file2) + '.txt'
    overlap_file = os.path.join(os.getcwd(), overlap_file)
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([labels_file1, labels_file2, overlap_file])
    cli.cmdline
    cli.run()

    return overlap_file


def measure_volume_overlap(labels, file1, file2, overlap_file):
    """
    Measure overlap between individual label regions
    in source and target nifti (nii.gz) images.

    Parameters
    ----------
    labels : list
        label indices
    file1 : string
        source image, consisting of index-labeled pixels/voxels
    file2 : string
        target image, consisting of index-labeled pixels/voxels
    overlap_file : string
        (optional) output file name

    Returns
    -------
    overlaps : numpy array
        overlap values
    overlap_file : string
        output text file name with overlap values

    Examples
    --------
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import measure_volume_overlap
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> file1 = os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> file2 = os.path.join(path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> dkt = DKTprotocol()
    >>> overlap_file = ''
    >>> measure_volume_overlap(dkt.label_numbers, file1, file2, overlap_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    save_output = True

    # Load labeled image volumes
    file1_data = nb.load(file1).get_data().ravel()
    file2_data = nb.load(file2).get_data().ravel()
    #print('Files: ' + file1 + ' ' + file2)
    #print(np.unique(file1_data))
    #print(np.unique(file2_data))

    # Initialize output
    overlaps = np.zeros((len(labels), 3))
    if save_output:
        file1_name = os.path.splitext(os.path.basename(file1))[0]
        file2_name = os.path.splitext(os.path.basename(file2))[0]
        if not overlap_file:
            overlap_file = os.path.join(os.getcwd(), 'labelvolume_dice_jacc_' +
                                        file2_name + '_vs_' + file1_name + '.txt')

    # Loop through labels
    for ilabel, label in enumerate(labels):
        label = int(label)
        overlaps[ilabel, 0] = label

        # Find which voxels contain the label in each volume
        file1_indices = np.where(file1_data==label)[0]
        file2_indices = np.where(file2_data==label)[0]
        file1_label_sum = len(file1_indices)
        file2_label_sum = len(file2_indices)

        # Determine their intersection and union
        intersect_label_sum = len(np.intersect1d(file2_indices, file1_indices))
        union_label_sum = len(np.union1d(file2_indices, file1_indices))
        #print('{0} {1} {2} {3} {4}'.format(label, file2_label_sum, file1_label_sum,
        #                              intersect_label_sum, union_label_sum))

        # There must be at least one voxel with the label in each volume
        if file2_label_sum * file1_label_sum > 0:

            # Compute Dice and Jaccard coefficients
            dice = np.float(2.0 * intersect_label_sum) / +\
                   (file2_label_sum + file1_label_sum)
            jacc = np.float(intersect_label_sum) / union_label_sum
            overlaps[ilabel, 1:3] = [dice, jacc]
            print('label: {0}, dice: {1:.2f}, jacc: {2:.2f}'.format(
                  label, dice, jacc))

    # Save output:
    if save_output:
        np.savetxt(overlap_file, overlaps, fmt='%d %.4f %.4f',
                   delimiter='\t', newline='\n')

    return overlaps, overlap_file


#-----------------------------------------------------------------------------
# Run evaluate_labels() on Mindboggle-101 data
# to compare manual and automated labels.
#-----------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    import numpy as np

    from mindboggle.mio.labels import DKTprotocol
    from mindboggle.evaluate.evaluate_labels import measure_volume_overlap
    from mindboggle.evaluate.evaluate_labels import measure_surface_overlap

    #-------------------------------------------------------------------------
    # Settings:
    #-------------------------------------------------------------------------
    label_method = 'freesurfer'  # 'ants'
    use_ants_segmentation = True  # False

    #-------------------------------------------------------------------------
    # File names, paths:
    #-------------------------------------------------------------------------
    if use_ants_segmentation:
        ants_str = ''
    else:
        ants_str = '_no_ants'
    if label_method == 'ants':
        volstem = 'ants_filled_labels'
    else:
        volstem = 'freesurfer_wmparc_filled_labels'
        if not use_ants_segmentation:
            volstem += '_in_freesurfer_segmentation'
    volfile = volstem + '.nii.gz'

    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20, 1,1,2,2,12]
    mindboggled = '/mnt/nfs-share/Mindboggle101/mindboggled/auto' + ants_str
    labels_dir = '/mnt/nfs-share/Mindboggle101/mindboggled/manual' + ants_str
    autosurf = 'freesurfer_cortex_labels.vtk'
    mansurf = 'relabeled_labels.DKT31.manual.vtk'

    #-------------------------------------------------------------------------
    # Evaluate surface labels:
    #-------------------------------------------------------------------------
    ccode_path = os.environ['MINDBOGGLE_TOOLS']
    command = os.path.join(ccode_path, 'surface_overlap',
                           'SurfaceOverlapMain')
    surfs = ['left_cortical_surface', 'right_cortical_surface']
    for iname, name in enumerate(names):
        number = numbers[iname]
        for n in range(1, number+1):
            subject = name+'-'+str(n)
            for isurf, surf in enumerate(surfs):
                mdir = os.path.join(mindboggled, subject, 'labels', surf)
                ldir = os.path.join(labels_dir, subject, 'labels', surf)
                file1 = os.path.join(mdir, autosurf)
                file2 = os.path.join(ldir, mansurf)
                overlap_file = "{0}_{1}_{2}_overlap.csv".\
                    format(subject, autosurf, mansurf)
#                measure_surface_overlap(command, file1, file2, overlap_file)

    #-------------------------------------------------------------------------
    # Evaluate volume labels:
    #-------------------------------------------------------------------------
    dkt = DKTprotocol()
    for iname, name in enumerate(names):
        number = numbers[iname]
        for n in range(1, number+1):
            subject = name+'-'+str(n)
            file1 = os.path.join(mindboggled, subject, 'labels', volfile)
            file2 = os.path.join(labels_dir, subject, 'labels', volfile)
            overlap_file = "{0}_{1}_overlap.csv".format(subject, volstem)
            measure_volume_overlap(dkt.label_numbers, file1, file2, overlap_file)
