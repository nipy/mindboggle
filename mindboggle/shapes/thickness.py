#!/usr/bin/python
"""
Compute a simple thickness measure for each labeled cortex region.


Authors:
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def thickinthehead(segmented_file, labeled_file, cortex_value=2,
                   noncortex_value=3, labels=[], out_dir='', resize=True,
                   propagate=True, output_table=False, use_c3d=False):
    """
    Compute a simple thickness measure for each labeled cortex region.

    Note: cortex, noncortex, labeled files are all the same coregistered brain.

    Steps ::

        1. Extract noncortex and cortex.
        2. Propagate labels through cortex (or simply multiply).
        3. Resample cortex and noncortex files.
        4. Extract outer and inner borders of cortex voxels.
        5. Optionally call ImageMath to propagate labels to fill cortex mask.
        6. Compute thickness as a ratio of label volume and edge volume.

    Parameters
    ----------
    segmented_file : str
        image volume with cortex and noncortex labels
    labeled_file : str
        corresponding image volume with index labels
    cortex_value : integer
        cortex label value in segmented_file
    noncortex_value : integer
        noncortex label value in segmented_file
    labels : list of integers
        label indices
    out_dir : str
        output directory
    resize : Boolean
        resize (2x) segmented_file for more accurate thickness estimates?
    propagate : Boolean
        propagate labels through cortex?
    output_table : False
        output a table with labels and thickness values?
    use_c3d : Boolean
        use convert3d? (otherwise ANTs ImageMath)

    Returns
    -------
    thickness_table : numpy array containing integers and floats
        label indices, volumes, areas, and thickness values
    table_file : string
        name of output table file

    Examples
    --------
    >>> from mindboggle.shapes.thickness import thickinthehead
    >>> segmented_file = '/Users/arno/Data/antsCorticalThickness/OASIS-TRT-20-1/tmp23314/tmpBrainSegmentation.nii.gz'
    >>> labeled_file = 'labels.nii.gz'
    >>> cortex_value = 2
    >>> noncortex_value = 3
    >>> #labels = [2]
    >>> labels = range(1002,1036) + range(2002,2036)
    >>> labels.remove(1004)
    >>> labels.remove(2004)
    >>> labels.remove(1032)
    >>> labels.remove(2032)
    >>> labels.remove(1033)
    >>> labels.remove(2033)
    >>> out_dir = '.'
    >>> resize = True
    >>> propagate = False
    >>> output_table = True
    >>> use_c3d = False
    >>> thickness_table, table_file = thickinthehead(segmented_file, labeled_file, cortex_value, noncortex_value, labels, out_dir, resize, propagate, output_table, use_c3d)

    """
    import os
    import numpy as np
    import nibabel as nb

    from mindboggle.utils.utils import execute

    #-------------------------------------------------------------------------
    # Output files:
    #-------------------------------------------------------------------------
    if out_dir:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
    else:
        out_dir = os.getcwd()
    cortex = os.path.join(out_dir, 'cortex.nii.gz')
    noncortex = os.path.join(out_dir, 'noncortex.nii.gz')
    temp = os.path.join(out_dir, 'temp.nii.gz')
    inner_edge = os.path.join(out_dir, 'cortex_inner_edge.nii.gz')
    use_outer_edge = True
    if use_outer_edge:
        outer_edge = os.path.join(out_dir, 'cortex_outer_edge.nii.gz')

    #-------------------------------------------------------------------------
    # Extract noncortex and cortex:
    #-------------------------------------------------------------------------
    if use_c3d:
        cmd = ['c3d', segmented_file,
               '-threshold', str(noncortex_value), str(noncortex_value), '1 0',
               '-o', noncortex]
        execute(cmd)
        cmd = ['c3d', segmented_file,
               '-threshold', str(cortex_value), str(cortex_value), '1 0',
               '-o', cortex]
        execute(cmd)
    else:
        cmd = ['ThresholdImage 3', segmented_file,
               noncortex, str(noncortex_value), str(noncortex_value), '1 0']
        execute(cmd)
        cmd = ['ThresholdImage 3', segmented_file,
               cortex, str(cortex_value), str(cortex_value), '1 0']
        execute(cmd)

    #-------------------------------------------------------------------------
    # Propagate labels through cortex (or simply multiply):
    #-------------------------------------------------------------------------
    if propagate and not use_c3d:
        cmd = ['ImageMath', '3', cortex, 'PropagateLabelsThroughMask',
               cortex, labeled_file]
        execute(cmd)
    else:
        if use_c3d:
            cmd = ['c3d', cortex, labeled_file, '-multiply', '-o', cortex]
            execute(cmd)
        else:
            cmd = ['ImageMath 3', cortex, 'm', cortex, labeled_file]
            execute(cmd)

    #-------------------------------------------------------------------------
    # Resample cortex and noncortex files:
    #-------------------------------------------------------------------------
    if resize:
        rescale = 2.0
        rescale_percent = str(rescale * 100)

        if use_c3d:
            cmd = ['c3d', cortex, '-interpolation nearestneighbor',
                   '-resample '+rescale_percent+'%', '-o', cortex]
            execute(cmd)
            cmd = ['c3d', noncortex, '-interpolation nearestneighbor',
                   '-resample '+rescale_percent+'%', '-o', noncortex]
            execute(cmd)
        else:
            dims = ' '.join([str(1/rescale), str(1/rescale), str(1/rescale)])
            cmd = ['ResampleImageBySpacing 3', cortex, cortex, dims, '0 0 1']
            execute(cmd)
            cmd = ['ResampleImageBySpacing 3', noncortex, noncortex, dims, '0 0 1']
            execute(cmd)

    #-------------------------------------------------------------------------
    # Extract cortex inner and outer border voxels:
    #-------------------------------------------------------------------------
    if use_c3d:
        cmd = ['c3d', noncortex, '-dilate 1 1x1x1vox', cortex, '-multiply',
               '-o', inner_edge]
        execute(cmd)
        if use_outer_edge:
            cmd = ['c3d', cortex, '-binarize -erode 1 1x1x1vox',
                   '-threshold 1 1 0 1', cortex, '-multiply', '-o', outer_edge]
            execute(cmd)
            cmd = ['c3d', inner_edge, '-binarize -threshold 1 1 0 1',
                   outer_edge, '-multiply', '-o', outer_edge]
            execute(cmd)
    else:
        cmd = ['ImageMath 3', inner_edge, 'MD', noncortex, '1']
        execute(cmd)
        cmd = ['ImageMath 3', inner_edge, 'm', cortex, inner_edge]
        execute(cmd)
        if use_outer_edge:
            cmd = ['ThresholdImage 3', cortex, outer_edge, '1 10000 1 0']
            execute(cmd)
            cmd = ['ImageMath 3', outer_edge, 'ME', outer_edge, '1']
            execute(cmd)
            cmd = ['ThresholdImage 3', outer_edge, outer_edge, '1 1 0 1']
            execute(cmd)
            cmd = ['ImageMath 3', outer_edge, 'm', cortex, outer_edge]
            execute(cmd)
            cmd = ['ThresholdImage 3', inner_edge, temp, '1 10000 1 0']
            execute(cmd)
            cmd = ['ThresholdImage 3', temp, temp, '1 1 0 1']
            execute(cmd)
            cmd = ['ImageMath 3', outer_edge, 'm', temp, outer_edge]
            execute(cmd)

    #-------------------------------------------------------------------------
    # Load data:
    #-------------------------------------------------------------------------
    compute_real_volume = True
    if compute_real_volume:
        img = nb.load(cortex)
        hdr = img.get_header()
        vv = np.prod(hdr.get_zooms())
        cortex_data = img.get_data().ravel()
    else:
        vv = 1
        cortex_data = nb.load(cortex).get_data().ravel()
    inner_edge_data = nb.load(inner_edge).get_data().ravel()
    if use_outer_edge:
        outer_edge_data = nb.load(outer_edge).get_data().ravel()

    #-------------------------------------------------------------------------
    # Loop through labels:
    #-------------------------------------------------------------------------
    if not labels:
        labeled_data = nb.load(labeled_file).get_data().ravel()
        labels = np.unique(labeled_data)
    labels = [int(x) for x in labels]
    if output_table:
        label_volume_area_thickness = np.zeros((len(labels), 4))
        label_volume_area_thickness[:,0] = labels
    for ilabel, label in enumerate(labels):

        #---------------------------------------------------------------------
        # Compute thickness as a ratio of label volume and edge volume:
        #---------------------------------------------------------------------
        label_cortex_volume = vv * len(np.where(cortex_data==label)[0])
        label_inner_edge_volume = vv * len(np.where(inner_edge_data==label)[0])
        if label_inner_edge_volume:
            if use_outer_edge:
                label_outer_edge_volume = \
                    vv * len(np.where(outer_edge_data==label)[0])
                label_area = (label_inner_edge_volume +
                              label_outer_edge_volume) / 2.0
            else:
                label_area = label_inner_edge_volume
            thickness = label_cortex_volume / label_area
            if output_table:
                label_volume_area_thickness[ilabel,1] = label_cortex_volume
                label_volume_area_thickness[ilabel,2] = label_area
                label_volume_area_thickness[ilabel,3] = thickness
            print('label {0} volume: cortex={1:2.2f}, inner={2:2.2f}, '
                  'outer={3:2.2f}, area={4:2.2f}, thickness={5:2.2f}mm'.
                  format(label, label_cortex_volume, label_inner_edge_volume,
                  label_outer_edge_volume, label_area, thickness))

    if output_table:
        table_file = os.path.join(out_dir, 'label_volume_area_thickness.csv')
        np.savetxt(table_file, label_volume_area_thickness,
                   fmt='%d %2.4f %2.4f %2.4f', delimiter='\t', newline='\n')
    else:
        table_file = None

    return label_volume_area_thickness, table_file


def run_thickinthehead(subjects, labels, out_dir='', atropos_dir='',
                       atropos_stem='', label_dir='', label_filename=''):
    """
    Combine FreeSurfer volume outputs (no surface-based outputs) to obtain
    a table of simple thickness measures for each labeled region.

    Steps ::

        1. Convert FreeSurfer volumes to nifti format in their original space.
        2. Combine FreeSurfer (and optionally ANTs) cortex & noncortex segmentations.
        3. Compute simple thickness per (optionally non-FreeSurfer) label.
        4. Store thickness values in a table.

    Requires ::

        FreeSurfer's mri_vol2vol

    Parameters
    ----------
    subjects : list of strings
        names of FreeSurfer subject directories
    labels : list of integers
        label indices
    out_dir : str (optional)
        output directory
    atropos_dir : str (optional)
        directory containing subject subdirectories with cortex files
    atropos_stem : str (optional)
        stem prepending name of antsCorticalThickness.sh output files
    label_dir : str (optional)
        directory containing subject subdirectories with label files
    label_filename : str (optional)
        label file name within <label_dir>/<subject>/

    Returns
    -------
    thickness_table : numpy array
        thickness values
    table_file : string
        name of file containing thickness values

    Examples
    --------
    >>> from mindboggle.shapes.thickness import run_thickinthehead
    >>> subjects = ['OASIS-TRT-20-1']
    >>> labels = range(1002,1036) + range(2002,2036)
    >>> labels.remove(1004)
    >>> labels.remove(2004)
    >>> labels.remove(1032)
    >>> labels.remove(2032)
    >>> labels.remove(1033)
    >>> labels.remove(2033)
    >>> out_dir = 'thickness_outputs'
    >>> atropos_dir = ''
    >>> atropos_stem = ''
    >>> label_dir = ''
    >>> label_filename = ''
    >>> thickness_table, table_file = run_thickinthehead(subjects, labels, out_dir, atropos_dir, atropos_stem, label_dir, label_filename)

    """
    import os, sys
    import numpy as np
    import nibabel as nb

    from mindboggle.utils.freesurfer import combine_segmentations
    from mindboggle.shapes.thickness import thickinthehead
    from mindboggle.utils.utils import execute

    subjects_dir = os.environ['SUBJECTS_DIR']

    use_c3d = False
    resize = True
    thickness_table = np.zeros((len(labels), len(subjects)+1))
    thickness_table[:,0] = labels
    for isubject, subject in enumerate(subjects):

        #---------------------------------------------------------------------
        # Outputs:
        #---------------------------------------------------------------------
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        out_subdir = os.path.join(out_dir, subject)
        if not os.path.exists(out_subdir):
            os.mkdir(out_subdir)
        labeled_file = os.path.join(out_subdir, 'labeled.nii.gz')

        #---------------------------------------------------------------------
        # Convert FreeSurfer label volume to nifti format (if no label file):
        #---------------------------------------------------------------------
        if label_filename:
            labeled_file = os.path.join(label_dir, subject, label_filename)
        else:
            aparc = os.path.join(subjects_dir, subject,
                                 'mri', 'aparc+aseg.mgz')
            orig = os.path.join(subjects_dir, subject,
                                'mri', 'orig', '001.mgz')
            if not os.path.exists(orig):
                orig = os.path.join(subjects_dir, subject,
                                    'mri', 'rawavg.mgz')
                if not os.path.exists(orig):
                    sys.exit('Could not find ' + orig +
                             ' for subject ' + subject)

            cmd = ['mri_vol2vol --mov', aparc, '--targ', orig,
                   '--interp nearest', '--regheader --o', labeled_file]
            execute(cmd)
            cmd = ['c3d', labeled_file, '-replace 2 0 41 0',
                   '-o', labeled_file]
            execute(cmd)

        #---------------------------------------------------------------------
        # Combine FreeSurfer and Atropos segmentations:
        #---------------------------------------------------------------------
        second_segmentation_file = ''
        aseg = ''
        filled = ''
        cortex_value = 2
        noncortex_value = 3
        segmented_file = combine_segmentations(subject, aseg, filled,
            out_subdir, second_segmentation_file, cortex_value,
            noncortex_value, use_c3d)
        #---------------------------------------------------------------------
        # Tabulate thickness values:
        #---------------------------------------------------------------------
        propagate = True
        output_table = False
        label_volume_area_thickness, u1 = thickinthehead(segmented_file,
                                    labeled_file, cortex_value, 
                                    noncortex_value, labels, out_subdir, 
                                    resize, propagate, output_table, use_c3d)

        thickness_table[:, isubject+1] = label_volume_area_thickness[3]

    #-------------------------------------------------------------------------
    # Save results:
    #-------------------------------------------------------------------------
    table_file = os.path.join(out_dir, 'thicknesses.csv')
    formatting = ' '.join(['%2.4f' for x in subjects])
    np.savetxt(table_file, thickness_table, fmt='%d ' + formatting,
               delimiter='\t', newline='\n')

    return thickness_table, table_file


if __name__ == "__main__":

    from mindboggle.shapes.thickness import run_thickinthehead

    subjects = []
    for i in range(1,21):
        subjects.append('OASIS-TRT-20-'+str(i))
        subjects.append('OASIS-TRT2-20-'+str(i))

    #-------------------------------------------------------------------------
    # Labels:
    #-------------------------------------------------------------------------
    labels = range(1002,1036) + range(2002,2036)
    labels.remove(1004)
    labels.remove(2004)
    labels.remove(1032)
    labels.remove(2032)
    labels.remove(1033)
    labels.remove(2033)

    out_dir = 'thickness_outputs_fs-ants'
    atropos_dir = '/homedir/Data/Brains/OASIS-TRT-20/antsCorticalThickness'
    atropos_stem = 'tmp'
    label_dir = ''
    label_filename = ''
    thickness_table, table_file = run_thickinthehead(subjects,
        labels, out_dir, atropos_dir, atropos_stem, label_dir, label_filename)
