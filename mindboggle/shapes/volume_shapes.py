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

    Computing the volume per labeled region is very straightforward:
    this function simply multiplies the volume per voxel by the number
    of voxels per region.

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
    save_table : bool
        save output table file with labels and volume values?
    output_table : string
        name of output table file with labels and volume values
    verbose : bool
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
    >>> input_file = fetch_data(urls['freesurfer_labels'], '', '.nii.gz')
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
    >>> [np.float("{0:.{1}f}".format(x, 5))
    ...  for x in [y for y in volumes if y > 0][0:5]]
    [971.99797, 2413.99496, 2192.99543, 8328.98262, 2940.99386]
    >>> [np.float("{0:.{1}f}".format(x, 5))
    ...  for x in [y for y in volumes if y > 0][5:10]]
    [1997.99583, 10905.97725, 11318.97639, 10789.97749, 2700.99437]

    """
    import os
    import numpy as np
    import nibabel as nb
    from io import open

    from mindboggle.guts.compute import count_per_label

    # Load labeled image volumes:
    img = nb.load(input_file)
    volume_per_voxel = np.product(img.header.get_zooms())
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
        fid = open(output_table, 'w', encoding='utf-8')
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


def thickinthehead(segmented_file, labeled_file,
                   cortex_value=2, noncortex_value=3, labels=[], names=[],
                   propagate=False, output_dir='', save_table=False,
                   output_table='', verbose=False):
    """
    Compute a simple thickness measure for each labeled cortex region volume.

    Since Mindboggle accepts FreeSurfer data as input, we include FreeSurfer
    cortical thickness estimates with Mindboggle’s shape measures.
    However, surface mesh reconstruction from MRI data does not always
    produce favorable results. For example, we found that at least a quarter
    of the over one hundred EMBARC brain images we processed through
    FreeSurfer clipped ventral cortical regions, resulting in bad surface
    patches in those regions. For comparison, we built this function called
    thickinthehead which computes a simple thickness measure for each
    cortical region using a segmentation volume rather than surfaces.

    We have revised this algorithm from the original published version.
    We removed upsampling to reduce memory issues for large image volumes,
    and replaced the estimated volume of middle cortical layer
    with an estimate of its surface area. We made these revisions to be less
    susceptible to deviations in voxel size from isometric 1mm^3 voxels
    for which thickinthehead was originally built.

    Steps ::

        1. Extract noncortex and cortex into separate files.
        2. Either mask labels with cortex or fill cortex with labels.
        3. Extract outer and inner boundary voxels of the cortex,
           by morphologically eroding the cortex (=2) by one voxel bordering
           the outside of the brain (=0) and bordering the inside of the brain
           (non-cortex=3).
        4. Estimate middle cortical layer's surface area by the average
           surface area of the outer and inner boundary voxels of the cortex,
           where surface area is roughly estimated as the average face area
           of a voxel times the number of voxels.
        5. Compute the volume of a labeled region of cortex.
        6. Estimate the thickness of the labeled cortical region as the
           volume of the labeled region (#5) divided by the
           estimate of the middle cortical surface area of that region (#4).

    Note::

      - Cortex, noncortex, & label files are from the same coregistered brain.
      - Calls ANTs functions: ImageMath and Threshold
      - There may be slight discrepancies between volumes computed by
        thickinthehead() and volumes computed by volume_per_label();
        in 31 of 600+ ADNI 1.5T images, some volume_per_label() volumes
        were slightly larger (in the third decimal place), presumably due to
        label propagation through the cortex in thickinthehead().
        This is more pronounced in ANTs vs. FreeSurfer-labeled volumes.

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
    propagate : bool
        propagate labels through cortex (or mask labels with cortex)?
    output_dir : string
        output directory
    save_table : bool
        save output table file with label volumes and thickness values?
    output_table : string
        name of output table file with label volumes and thickness values
    verbose : bool
        print statements?

    Returns
    -------
    label_volume_thickness : list of lists of integers and floats
        label indices, volumes, and thickness values (default -1)
    output_table : string
        name of output table file with label volumes and thickness values

    Examples
    --------
    >>> # Example simply using ants segmentation and labels vs. hybrid segmentation:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.shapes.volume_shapes import thickinthehead
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> segmented_file = fetch_data(urls['ants_segmentation'], '', '.nii.gz')
    >>> labeled_file = fetch_data(urls['ants_labels'], '', '.nii.gz')
    >>> cortex_value = 2
    >>> noncortex_value = 3
    >>> #labels = [2]
    >>> labels = list(range(1002,1036)) + list(range(2002,2036))
    >>> labels.remove(1004)
    >>> labels.remove(2004)
    >>> labels.remove(1032)
    >>> labels.remove(2032)
    >>> labels.remove(1033)
    >>> labels.remove(2033)
    >>> names = []
    >>> propagate = False
    >>> output_dir = ''
    >>> save_table = True
    >>> output_table = ''
    >>> verbose = False

    Skip online test because it requires installation of ANTs:

    >>> label_volume_thickness, output_table = thickinthehead(segmented_file,
    ...     labeled_file, cortex_value, noncortex_value, labels, names,
    ...     propagate, output_dir, save_table, output_table, verbose) # doctest: +SKIP
    >>> [np.int("{0:.{1}f}".format(x, 5)) label_volume_thickness[0][0:10]] # doctest: +SKIP
    >>> [np.float("{0:.{1}f}".format(x, 5)) for x in label_volume_thickness[1][0:5]] # doctest: +SKIP
    >>> [np.float("{0:.{1}f}".format(x, 5)) for x in label_volume_thickness[2][0:5]] # doctest: +SKIP

    """
    import os
    import numpy as np
    import nibabel as nb
    from io import open

    from mindboggle.guts.utilities import execute

    # ------------------------------------------------------------------------
    # Output files:
    # ------------------------------------------------------------------------
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
        fid = open(output_table, 'w', encoding='utf-8')
        if names:
            fid.write("name, ID, thickness (thickinthehead)\n")
        else:
            fid.write("ID, thickness (thickinthehead)\n")
    else:
        output_table = ''

    # ------------------------------------------------------------------------
    # Extract noncortex and cortex:
    # ------------------------------------------------------------------------
    cmd = ['ThresholdImage', '3', segmented_file, noncortex,
           str(noncortex_value), str(noncortex_value), '1 0']
    execute(cmd, 'os')
    cmd = ['ThresholdImage', '3', segmented_file, cortex,
           str(cortex_value), str(cortex_value), '1 0']
    execute(cmd, 'os')

    # ------------------------------------------------------------------------
    # Either mask labels with cortex or fill cortex with labels:
    # ------------------------------------------------------------------------
    if propagate:
        cmd = ['ImageMath', '3', cortex, 'PropagateLabelsThroughMask',
               cortex, labeled_file]
        execute(cmd, 'os')
    else:
        cmd = ['ImageMath', '3', cortex, 'm', cortex, labeled_file]
        execute(cmd, 'os')

    # ------------------------------------------------------------------------
    # Load data and dimensions:
    # ------------------------------------------------------------------------
    img = nb.load(cortex)
    cortex_data = img.get_data().ravel()
    voxsize = img.header.get_zooms()
    voxvol = np.prod(voxsize)
    voxarea = (voxsize[0] * voxsize[1] + \
               voxsize[0] * voxsize[2] + \
               voxsize[1] * voxsize[2]) / 3

    # ------------------------------------------------------------------------
    # Extract outer and inner boundary voxels of the cortex,
    # by eroding 1 voxel for cortex voxels (=2) bordering
    # the outside of the brain (=0) and bordering noncortex (=3):
    # ------------------------------------------------------------------------
    cmd = ['ImageMath', '3', inner_edge, 'MD', noncortex, '1']
    execute(cmd, 'os')
    cmd = ['ImageMath', '3', inner_edge, 'm', cortex, inner_edge]
    execute(cmd, 'os')
    if use_outer_edge:
        cmd = ['ThresholdImage', '3', cortex, outer_edge, '1 10000 1 0']
        execute(cmd, 'os')
        cmd = ['ImageMath', '3', outer_edge, 'ME', outer_edge, '1']
        execute(cmd, 'os')
        cmd = ['ThresholdImage', '3', outer_edge, outer_edge, '1 1 0 1']
        execute(cmd, 'os')
        cmd = ['ImageMath', '3', outer_edge, 'm', cortex, outer_edge]
        execute(cmd, 'os')
        cmd = ['ThresholdImage', '3', inner_edge, temp, '1 10000 1 0']
        execute(cmd, 'os')
        cmd = ['ThresholdImage', '3', temp, temp, '1 1 0 1']
        execute(cmd, 'os')
        cmd = ['ImageMath', '3', outer_edge, 'm', temp, outer_edge]
        execute(cmd, 'os')

    # ------------------------------------------------------------------------
    # Load data:
    # ------------------------------------------------------------------------
    inner_edge_data = nb.load(inner_edge).get_data().ravel()
    if use_outer_edge:
        outer_edge_data = nb.load(outer_edge).get_data().ravel()

    # ------------------------------------------------------------------------
    # Loop through labels:
    # ------------------------------------------------------------------------
    if not labels:
        labeled_data = nb.load(labeled_file).get_data().ravel()
        labels = np.unique(labeled_data)
    labels = [int(x) for x in labels]
    label_volume_thickness = -1 * np.ones((len(labels), 3))
    label_volume_thickness[:, 0] = labels
    for ilabel, label in enumerate(labels):
        if names:
            name = names[ilabel]

        # --------------------------------------------------------------------
        # Compute thickness as a ratio of label volume and layer surface area:
        #   - Estimate middle cortical surface area by the average area
        #     of the outer and inner boundary voxels of the cortex.
        #   - Surface area is roughly estimated as the average face area
        #     of a voxel times the number of voxels.
        #   - Compute the volume of a labeled region of cortex.
        #   - Estimate the thickness of the labeled cortical region as the
        #     volume of the labeled region divided by the middle surface area.
        # --------------------------------------------------------------------
        label_cortex_volume = voxvol * len(np.where(cortex_data==label)[0])
        label_inner_edge_area = voxarea * \
                                  len(np.where(inner_edge_data==label)[0])
        if label_inner_edge_area:
            if use_outer_edge:
                label_outer_edge_area = \
                    voxarea * len(np.where(outer_edge_data==label)[0])
                label_area = (label_inner_edge_area +
                              label_outer_edge_area) / 2.0
            else:
                label_area = label_inner_edge_area
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


# ============================================================================
# Doctests
# ============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules