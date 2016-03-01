#!/usr/bin/env python
"""
Compute shape measures from volume images.


Authors:
    - Arno Klein, 2013-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def volume_per_brain_region(input_file, include_labels=[], exclude_labels=[],
                            label_names=[], save_table=False,
                            output_table='', verbose=False):
    """
    Compute volume per labeled region in a nibabel-readable image.

    Note: Results are truncated at three decimal places because we found that
    volume label propagation led to differences in the third decimal place.

    Parameters
    ----------
    input_file : string
        name of image file, consisting of index-labeled pixels/voxels
    include_labels : list of integers
        labels to include
        (if empty, use unique numbers in image volume file)
    exclude_labels : list of integers
        labels to be excluded
    label_names : list of strings
        label names corresponding to labels to include
    save_table : Boolean
        save output table file with labels and volume values?
    output_table : string
        name of output table file with labels and volume values
    verbose : Boolean
        print statements?

    Returns
    -------
    unique_labels : list of integers
        unique label numbers (default -1)
    volumes : list of floats
        volume for each label (default -1)
    output_table : string
        name of output volume table file (if output_table not empty)

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.shapes.volume_shapes import volume_per_brain_region
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_file = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(input_file, input_file + '.nii.gz')
    >>> input_file += '.nii.gz'
    >>> dkt = DKTprotocol()
    >>> include_labels = dkt.label_numbers
    >>> exclude_labels = []
    >>> label_names = dkt.label_names
    >>> save_table = True
    >>> output_table = 'volumes.csv'
    >>> verbose = False
    >>> unique_labels, volumes, table = volume_per_brain_region(input_file,
    ...     include_labels, exclude_labels, label_names, save_table,
    ...     output_table, verbose)
    >>> print(np.array_str(np.array(volumes[0:5]),
    ...       precision=5, suppress_small=True))
    [  971.99797  2413.99496  2192.99543  8328.98262  2940.99386]
    >>> print(np.array_str(np.array(volumes[5:10]),
    ...       precision=5, suppress_small=True))
    [  1997.99583  10905.97725  11318.97639  10789.97749   2700.99437]

    """
    import os
    import numpy as np
    import nibabel as nb

    from mindboggle.guts.compute import count_per_label

    # Load labeled image volumes:
    img = nb.load(input_file)
    hdr = img.get_header()
    volume_per_voxel = np.product(hdr.get_zooms())
    labels = img.get_data().ravel()

    unique_labels, counts = count_per_label(labels, include_labels,
                                            exclude_labels)
    volumes = [volume_per_voxel * x for x in counts]

    # Output table:
    if save_table:
        if output_table:
            output_table = os.path.join(os.getcwd(), output_table)
        else:
            output_table = os.path.join(os.getcwd(),
                                        'volume_for_each_label.csv')
        fid = open(output_table, 'w')
        if len(label_names) == len(unique_labels):
            fid.write("name, ID, volume\n")
        else:
            fid.write("ID, volume\n")

        # Loop through labels:
        for ilabel, label in enumerate(unique_labels):
            if volumes[ilabel]:

                if len(label_names) == len(unique_labels):
                    if verbose:
                        print('{0} ({1}) volume = {2:2.3f}mm^3\n'.format(
                              label_names[ilabel], label, volumes[ilabel]))
                    fid.write('{0}, {1}, {2:2.3f}\n'.format(
                              label_names[ilabel], label, volumes[ilabel]))
                else:
                    if verbose:
                        print('{0} volume = {1:2.3f}mm^3\n'.format(
                              label, volumes[ilabel]))
                    fid.write('{0}, {1:2.3f}\n'.format(label, volumes[ilabel]))
    else:
        output_table = ''

    return unique_labels, volumes, output_table


def thickinthehead(segmented_file, labeled_file, cortex_value=2,
                   noncortex_value=3, labels=[], names=[], resize=True,
                   propagate=True, output_dir='', save_table=False,
                   output_table='', ants_path='', verbose=False):
    """
    Compute a simple thickness measure for each labeled cortex region volume.

    Note::

      - Cortex, noncortex, & label files are from the same coregistered brain.
      - Calls ANTs functions: ImageMath, Threshold, ResampleImageBySpacing
      - There may be slight discrepancies between volumes computed by
        thickinthehead() and volumes computed by volume_for_each_label();
        in 31 of 600+ ADNI 1.5T images, some volume_for_each_label() volumes
        were slightly larger (in the third decimal place), presumably due to
        label propagation through the cortex in thickinthehead().
        This is more pronounced in ANTs vs. FreeSurfer-labeled volumes.

    Example preprocessing steps ::

      1. Run Freesurfer and antsCorticalThickness.sh on T1-weighted image.
      2. Convert FreeSurfer volume labels (e.g., wmparc.mgz or aparc+aseg.mgz)
         to cortex (2) and noncortex (3) segments using relabel_volume()
         function [refer to LABELS.rst or FreeSurferColorLUT labels file].
      3. Convert ANTs Atropos-segmented volume (tmpBrainSegmentation.nii.gz)
         to cortex and noncortex segments, by converting 1-labels to 0 and
         4-labels to 3 with the relabel_volume() function
         (the latter is to include deep-gray matter with noncortical tissues).
      4. Combine FreeSurfer and ANTs segmentation volumes to obtain a single
         cortex (2) and noncortex (3) segmentation file using the function
         combine_2labels_in_2volumes(). This function takes the union of
         cortex voxels from the segmentations, the union of the noncortex
         voxels from the segmentations, and overwrites intersecting cortex
         and noncortex voxels with noncortex (3) labels.
         ANTs tends to include more cortical gray matter at the periphery of
         the brain than Freesurfer, and FreeSurfer tends to include more white
         matter that extends deep into gyral folds than ANTs, so the above
         attempts to remedy their differences by overlaying ANTs cortical gray
         with FreeSurfer white matter.
      5. Optional, see Step 2 below:
         Fill segmented cortex with cortex labels and noncortex with
         noncortex labels using the PropagateLabelsThroughMask() function
         (which calls ImageMath ... PropagateLabelsThroughMask in ANTs).
         The labels can be initialized using FreeSurfer (e.g. wmparc.mgz)
         or ANTs (by applying the nonlinear inverse transform generated by
         antsCorticalThickness.sh to labels in the Atropos template space).
         [Note: Any further labeling steps may be applied, such as
         overwriting cerebrum with intersecting cerebellum labels.]

    Steps ::

        1. Extract noncortex and cortex.
        2. Either mask labels with cortex or fill cortex with labels.
        3. Resample cortex and noncortex files from 1x1x1 to 0.5x0.5x0.5
           to better represent the contours of the boundaries of the cortex.
        4. Extract outer and inner boundary voxels of the cortex,
           by eroding 1 (resampled) voxel for cortex voxels (2) bordering
           the outside of the brain (0) and bordering noncortex (3).
        5. Estimate middle cortical surface area by the average volume
           of the outer and inner boundary voxels of the cortex.
        6. Compute the volume of a labeled region of cortex.
        7. Estimate the thickness of the labeled cortical region as the
           volume of the labeled region (#6) divided by the surface area (#5).

    Parameters
    ----------
    segmented_file : string
        image volume with cortex and noncortex (and any other) labels
    labeled_file : string
        corresponding image volume with index labels
    cortex_value : integer
        cortex label value in segmented_file
    noncortex_value : integer
        noncortex label value in segmented_file
    labels : list of integers
        label indices
    names : list of strings
        label names
    resize : Boolean
        resize (2x) segmented_file for more accurate thickness estimates?
    propagate : Boolean
        propagate labels through cortex?
    output_dir : string
        output directory
    save_table : Boolean
        save output table file with label volumes and thickness values?
    output_table : string
        name of output table file with label volumes and thickness values
    verbose : Boolean
        print statements?

    Returns
    -------
    label_volume_thickness : list of lists of integers and floats
        label indices, volumes, and thickness values (default -1)
    output_table : string
        name of output table file with label volumes and thickness values

    Examples
    --------
    >>> # Example simply using ants segmentation and labels:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.shapes.volume_shapes import thickinthehead
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> segmented_file = fetch_data(urls['ants_segmentation'])
    >>> os.rename(segmented_file, segmented_file + '.nii.gz')
    >>> segmented_file += '.nii.gz'
    >>> labeled_file = fetch_data(urls['ants_labels'])
    >>> os.rename(labeled_file, labeled_file + '.nii.gz')
    >>> labeled_file += '.nii.gz'
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
    >>> names = []
    >>> resize = True
    >>> propagate = False
    >>> output_dir = ''
    >>> save_table = True
    >>> output_table = ''
    >>> verbose = False
    >>> label_volume_thickness, output_table = thickinthehead(segmented_file,
    ...     labeled_file, cortex_value, noncortex_value, labels, names,
    ...     resize, propagate, output_dir, save_table, output_table, verbose)
    >>> print(np.array_str(np.array(label_volume_thickness[0][0:10]),
    ...       precision=5, suppress_small=True))
    [ 1002.  1003.  1005.  1006.  1007.  1008.  1009.  1010.  1011.  1012.]
    >>> print(np.array_str(np.array(label_volume_thickness[1][0:5]),
    ...       precision=5, suppress_small=True))
    [  3136.99383   7206.98582   3257.99359   1950.99616  12458.97549]
    >>> print(np.array_str(np.array(label_volume_thickness[2][0:5]),
    ...       precision=5, suppress_small=True))
    [ 3.8639   3.69637  2.56334  4.09336  4.52592]

    """
    import os
    import numpy as np
    import nibabel as nb

    from mindboggle.guts.utilities import execute

    #-------------------------------------------------------------------------
    # Output files:
    #-------------------------------------------------------------------------
    if output_dir:
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
    else:
        output_dir = os.getcwd()
    cortex = os.path.join(output_dir, 'cortex.nii.gz')
    noncortex = os.path.join(output_dir, 'noncortex.nii.gz')
    temp = os.path.join(output_dir, 'temp.nii.gz')
    inner_edge = os.path.join(output_dir, 'cortex_inner_edge.nii.gz')
    use_outer_edge = True
    if use_outer_edge:
        outer_edge = os.path.join(output_dir, 'cortex_outer_edge.nii.gz')

    if save_table:
        if output_table:
            output_table = os.path.join(os.getcwd(), output_table)
        else:
            output_table = os.path.join(os.getcwd(),
                                        'thickinthehead_for_each_label.csv')
        fid = open(output_table, 'w')
        if names:
            fid.write("name, ID, thickness (thickinthehead)\n")
        else:
            fid.write("ID, thickness (thickinthehead)\n")
    else:
        output_table = ''

    #-------------------------------------------------------------------------
    # ants command paths:
    #-------------------------------------------------------------------------
    ants_thresh = 'ThresholdImage'
    ants_math = 'ImageMath'
    ants_resample = 'ResampleImageBySpacing'

    #-------------------------------------------------------------------------
    # Extract noncortex and cortex:
    #-------------------------------------------------------------------------
    cmd = [ants_thresh, '3', segmented_file, noncortex,
           str(noncortex_value), str(noncortex_value), '1 0']
    execute(cmd, 'os')
    cmd = [ants_thresh, '3', segmented_file, cortex,
           str(cortex_value), str(cortex_value), '1 0']
    execute(cmd, 'os')

    #-------------------------------------------------------------------------
    # Either mask labels with cortex or fill cortex with labels:
    #-------------------------------------------------------------------------
    if propagate:
        cmd = [ants_math, '3', cortex, 'PropagateLabelsThroughMask',
               cortex, labeled_file]
        execute(cmd, 'os')
    else:
        cmd = [ants_math, '3', cortex, 'm', cortex, labeled_file]
        execute(cmd, 'os')

    #-------------------------------------------------------------------------
    # Load data and dimensions:
    #-------------------------------------------------------------------------
    if resize:
        rescale = 2.0
    else:
        rescale = 1.0
    compute_real_volume = True
    if compute_real_volume:
        img = nb.load(cortex)
        hdr = img.get_header()
        vv_orig = np.prod(hdr.get_zooms())
        vv = np.prod([x/rescale for x in hdr.get_zooms()])
        cortex_data = img.get_data().ravel()
    else:
        vv = 1/rescale
        cortex_data = nb.load(cortex).get_data().ravel()

    #-------------------------------------------------------------------------
    # Resample cortex and noncortex files from 1x1x1 to 0.5x0.5x0.5
    # to better represent the contours of the boundaries of the cortex:
    #-------------------------------------------------------------------------
    if resize:
        dims = ' '.join([str(1/rescale), str(1/rescale), str(1/rescale)])
        cmd = [ants_resample, '3', cortex, cortex, dims, '0 0 1']
        execute(cmd, 'os')
        cmd = [ants_resample, '3', noncortex, noncortex, dims, '0 0 1']
        execute(cmd, 'os')

    #-------------------------------------------------------------------------
    # Extract outer and inner boundary voxels of the cortex,
    # by eroding 1 (resampled) voxel for cortex voxels (2) bordering
    # the outside of the brain (0) and bordering noncortex (3):
    #-------------------------------------------------------------------------
    cmd = [ants_math, '3', inner_edge, 'MD', noncortex, '1']
    execute(cmd, 'os')
    cmd = [ants_math, '3', inner_edge, 'm', cortex, inner_edge]
    execute(cmd, 'os')
    if use_outer_edge:
        cmd = [ants_thresh, '3', cortex, outer_edge, '1 10000 1 0']
        execute(cmd, 'os')
        cmd = [ants_math, '3', outer_edge, 'ME', outer_edge, '1']
        execute(cmd, 'os')
        cmd = [ants_thresh, '3', outer_edge, outer_edge, '1 1 0 1']
        execute(cmd, 'os')
        cmd = [ants_math, '3', outer_edge, 'm', cortex, outer_edge]
        execute(cmd, 'os')
        cmd = [ants_thresh, '3', inner_edge, temp, '1 10000 1 0']
        execute(cmd, 'os')
        cmd = [ants_thresh, '3', temp, temp, '1 1 0 1']
        execute(cmd, 'os')
        cmd = [ants_math, '3', outer_edge, 'm', temp, outer_edge]
        execute(cmd, 'os')

    #-------------------------------------------------------------------------
    # Load data:
    #-------------------------------------------------------------------------
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
    label_volume_thickness = -1 * np.ones((len(labels), 3))
    label_volume_thickness[:, 0] = labels
    for ilabel, label in enumerate(labels):
        if names:
            name = names[ilabel]

        #---------------------------------------------------------------------
        # Compute thickness as a ratio of label volume and edge volume:
        #   - Estimate middle cortical surface area by the average volume
        #     of the outer and inner boundary voxels of the cortex.
        #   - Compute the volume of a labeled region of cortex.
        #   - Estimate the thickness of the labeled cortical region as the
        #     volume of the labeled region divided by the surface area.
        #---------------------------------------------------------------------
        label_cortex_volume = vv_orig * len(np.where(cortex_data==label)[0])
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
            label_volume_thickness[ilabel, 1] = label_cortex_volume
            label_volume_thickness[ilabel, 2] = thickness

            if save_table:
                if names:
                    if verbose:
                        print('{0} ({1}) thickinthehead thickness = '
                              '{2:2.2f}mm'.format(name, label, thickness))
                    fid.write('{0}, {1}, {2:2.3f}\n'.format(name, label,
                                                            thickness))
                else:
                    if verbose:
                        print('{0} thickinthehead thickness = {1:2.2f}mm'.
                              format(label, thickness))
                    fid.write('{0}, {1:2.3f}\n'.format(label, thickness))

    label_volume_thickness = label_volume_thickness.transpose().tolist()

    return label_volume_thickness, output_table


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()