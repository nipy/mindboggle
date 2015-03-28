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


def evaluate_volume_overlaps(labels, file1, file2, 
                             output_file='', save_output=True):
    """
    Compute overlap between individual label regions
    in source and target nifti (nii.gz) images.

    Parameters
    ----------
    labels : list
        label indices
    file1 : string
        source image, consisting of index-labeled pixels/voxels
    file2 : string
        target image, consisting of index-labeled pixels/voxels
    output_file : string
        (optional) output file name
    save_output : Boolean
        save output file?

    Returns
    -------
    dice_overlaps : numpy array
        Dice overlap values
    jacc_overlaps : numpy array
        Jaccard overlap values
    output_file : string
        output text file name with overlap values

    Examples
    --------
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import evaluate_volume_overlaps
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> path = '/homedir/mindboggled'
    >>> file1 = os.path.join(path, 'Twins-2-1', 'labels', 'freesurfer_wmparc_filled_labels.nii.gz')
    >>> file2 = os.path.join(path, 'Twins-2-1', 'labels', 'freesurfer_wmparc_filled_labels.nii.gz')
    >>> dkt = DKTprotocol()
    >>> labels = dkt.cerebrum_cortex_DKT31_numbers
    >>> output_file = ''
    >>> save_output = True
    >>> evaluate_volume_overlaps(labels, file1, file2, output_file=output_file, save_output=save_output)

    """
    import nibabel as nb

    from mindboggle.guts.compute import compute_overlaps

    # Load labeled image volumes:
    list1 = nb.load(file1).get_data().ravel()
    list2 = nb.load(file2).get_data().ravel()

    dice_overlaps, jacc_overlaps, output_file = compute_overlaps(labels,
        list1, list2, output_file=output_file, save_output=save_output)

    return dice_overlaps, jacc_overlaps, output_file


def evaluate_surface_overlaps(labels, index, table1, table2,
                              output_file='', save_output=True):
    """
    Measure surface overlap per label by comparing Mindboggle vertices tables.

    Parameters
    ----------
    labels : list
        cortical label indices to measure surface overlap
    index : integer
        index (starting from zero) to column of table containing label indices
    table1 : string
        table with index labels for scalar values
    table2 : string
        table with index labels for scalar values
    output_file : string
        (optional) output file name

    Returns
    -------
    dice_overlaps : numpy array
        Dice overlap values
    jacc_overlaps : numpy array
        Jaccard overlap values
    output_file : string
        (optional) output file name
    save_output : Boolean
        save output file?

    Examples
    --------
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import evaluate_surface_overlaps
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> dkt = DKTprotocol()
    >>> labels = dkt.cerebrum_cortex_DKT31_numbers
    >>> index = 0
    >>> path = '/homedir/mindboggled'
    >>> table1 = os.path.join(path, 'Twins-2-1', 'tables', 'left_cortical_surface', 'vertices.csv')
    >>> table2 = os.path.join(path, 'Twins-2-1', 'tables', 'left_cortical_surface', 'vertices.csv')
    >>> output_file = ''
    >>> save_output = True
    >>> evaluate_surface_overlaps(labels, index, table1, table2, output_file=output_file, save_output=save_output)

    """
    import pandas as pd

    from mindboggle.guts.compute import compute_overlaps

    # Load surface label tables:
    df1 = pd.read_csv(table1)
    df2 = pd.read_csv(table2)
    list1 = df1.iloc[:, index]
    list2 = df2.iloc[:, index]

    dice_overlaps, jacc_overlaps, output_file = compute_overlaps(labels,
        list1, list2, output_file=output_file, save_output=save_output)

    return dice_overlaps, jacc_overlaps, output_file


# def evaluate_surface_overlaps_cpp(command, labels_file1, labels_file2,
#                                   output_file):
#     """
#     Measure surface overlap using Joachim Giard's code.
#
#     Note: Currently failing on VTK versions of manually labeled files,
#     so evaluate_surface_overlaps() is included.
#
#     Parameters
#     ----------
#     command : string
#         surface overlap C++ executable command
#     labels_file1 : string
#         ``vtk file`` with index labels for scalar values
#     labels_file2 : string
#         ``vtk file`` with index labels for scalar values
#     output_file : string
#         (optional) output file name
#
#     Returns
#     -------
#     output_file : string
#         name of output text file with overlap results
#
#     Examples
#     --------
#     >>> import os
#     >>> from mindboggle.evaluate.evaluate_labels import evaluate_surface_overlaps_cpp
#     >>> from mindboggle.mindboggle import hashes_url
#     >>> from mindboggle.mio.fetch_data import fetch_check_data
#     >>> hashes, url, cache_env, cache = hashes_url()
#     >>> ccode_path = os.environ['MINDBOGGLE_TOOLS']
#     >>> command = os.path.join(ccode_path, 'surface_overlap', 'SurfaceOverlapMain')
#     >>> label_file1 = 'lh.labels.DKT25.manual.vtk'
#     >>> label_file2 = 'lh.labels.DKT31.manual.vtk'
#     >>> file1 = fetch_check_data(label_file1, url, hashes, cache_env, cache)
#     >>> file2 = fetch_check_data(label_file2, url, hashes, cache_env, cache)
#     >>> output_file = ''
#     >>> evaluate_surface_overlaps_cpp(command, file1, file2, output_file)
#
#     """
#     import os
#     from nipype.interfaces.base import CommandLine
#
#     if not output_file:
#         output_file = os.path.basename(labels_file1) + '_and_' + \
#                            os.path.basename(labels_file2) + '.txt'
#     output_file = os.path.join(os.getcwd(), output_file)
#     cli = CommandLine(command = command)
#     cli.inputs.args = ' '.join([labels_file1, labels_file2, output_file])
#     cli.cmdline
#     cli.run()
#
#     return output_file


#-----------------------------------------------------------------------------
# Run evaluate_labels.py on Mindboggle-101 data
# to compare manual and automated volume labels and surface labels.
#-----------------------------------------------------------------------------
if __name__ == "__main__":

    import os

    from mindboggle.mio.labels import DKTprotocol
    from mindboggle.evaluate.evaluate_labels import evaluate_volume_overlaps
    from mindboggle.evaluate.evaluate_labels import evaluate_surface_overlaps

    dkt = DKTprotocol()

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
                index = 0
                output_file = "{0}_{1}_{2}_volume_label_overlaps.csv".\
                    format(subject, autosurf, mansurf)
                evaluate_surface_overlaps(dkt.cerebrum_cortex_DKT31_numbers,
                                          index, file1, file2, output_file)

    #-------------------------------------------------------------------------
    # Evaluate volume labels:
    #-------------------------------------------------------------------------
    for iname, name in enumerate(names):
        number = numbers[iname]
        for n in range(1, number+1):
            subject = name+'-'+str(n)
            file1 = os.path.join(mindboggled, subject, 'labels', volfile)
            file2 = os.path.join(labels_dir, subject, 'labels', volfile)
            output_file = "{0}_{1}_surface_label_overlaps.csv".format(
                subject, volstem)
            evaluate_volume_overlaps(dkt.label_numbers,
                                     file1, file2, output_file)
