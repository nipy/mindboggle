#!/usr/bin/python
"""
Functions that call ANTs (UPenn's PICSL group) commands.

Authors:
Arno Klein, 2011-2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def fetch_ants_data(segmented_file):
    """
    Fetch antsCorticalThickness.sh output.

    The input argument "segmented_file" is one of the relevant
    antsCorticalThickness.sh output files called by Mindboggle:
        ants_subjects/subject1/antsBrainExtractionMask.nii.gz
        ants_subjects/subject1/antsBrainSegmentation.nii.gz
        ants_subjects/subject1/antsTemplateToSubject0GenericAffine.mat
        ants_subjects/subject1/antsTemplateToSubject1Warp.nii.gz

    Parameters
    ----------
    segmented_file : string
        full path to a subject's antsCorticalThickness.sh segmented file

    Returns
    -------
    mask : string
        antsBrainExtraction.sh brain volume mask for extracting brain volume
    segments : string
        Atropos-segmented brain volume
    affine : string
        subject to template affine transform (antsRegistration)
        Note: transform name contains "TemplateToSubject"
    warp : string
        subject to template nonlinear transform (antsRegistration)
        Note: transform name contains "TemplateToSubject"
    invwarp : string
        inverse of subject to template nonlinear transform

    Examples
    --------
    >>> from mindboggle.utils.ants import fetch_ants_data
    >>> segmented_file = 'ants_subjects/OASIS-TRT-20-1/tmpBrainSegmentation.nii.gz'
    >>> fetch_ants_data(segmented_file)

    """
    import os

    prefix = segmented_file.strip('BrainSegmentation.nii.gz')

    mask = prefix + 'BrainExtractionMask.nii.gz'
    segments = segmented_file
    affine = prefix + 'TemplateToSubject0GenericAffine.mat'
    warp = prefix + 'TemplateToSubject1Warp.nii.gz'
    invwarp = prefix + 'TemplateToSubject1InverseWarp.nii.gz'

    for s in [mask, segments, affine, warp, invwarp]:
        if not os.path.exists(s):
            raise IOError(s + ' does not exist')

    return mask, segments, affine, warp, invwarp


def ComposeMultiTransform(transform_files, inverse_Booleans,
                          output_transform_file='', ext='.txt'):
    """
    Run ANTs ComposeMultiTransform function to create a single transform.

    Parameters
    ----------
    transform_files : list of strings
        transform file names
    inverse_Booleans : list of Booleans
        Booleans to indicate which transforms to take the inverse of
    output_transform_file : string
        transform file name
    ext : string
        '.txt' to save transform file as text, '.mat' for data file

    Returns
    -------
    output_transform_file : string
        single composed transform file name

    Examples
    --------
    >>> from mindboggle.utils.ants import ComposeMultiTransform
    >>> transform_files = ['affine1.mat', 'affine2.mat']
    >>> transform_files = ['/data/Brains/Mindboggle101/antsCorticalThickness/OASIS-TRT-20_volumes/OASIS-TRT-20-1/antsTemplateToSubject0GenericAffine.mat','/data/Brains/Mindboggle101/antsCorticalThickness/OASIS-TRT-20_volumes/OASIS-TRT-20-1/antsTemplateToSubject0GenericAffine.mat']
    >>> inverse_Booleans = [False, False]
    >>> output_transform_file = ''
    >>> ext = '.txt'
    >>> ComposeMultiTransform(transform_files, inverse_Booleans, output_transform_file, ext)

    """
    import os

    from mindboggle.utils.utils import execute

    if not output_transform_file:
        output_transform_file = os.path.join(os.getcwd(), 'affine' + ext)

    xfms = []
    for ixfm, xfm in enumerate(transform_files):
        if inverse_Booleans[ixfm]:
            xfms.append('-i')
        xfms.append(xfm)

    cmd = ['ComposeMultiTransform 3', output_transform_file, ' '.join(xfms)]
    print(cmd)
    execute(cmd, 'os')
    #if not os.path.exists(output_transform_file):
    #    raise(IOError(output_transform_file + " not found"))

    return output_transform_file


def ImageMath(volume1, volume2, operator='m', output_file=''):
    """
    Use the ImageMath function in ANTS to perform operation on two volumes::

        m         : Multiply ---  use vm for vector multiply
        +         : Add ---  use v+ for vector add
        -         : Subtract ---  use v- for vector subtract
        /         : Divide
        ^         : Power
        exp       : Take exponent exp(imagevalue*value)
        addtozero : add image-b to image-a only over points where image-a has zero values
        overadd   : replace image-a pixel with image-b pixel if image-b pixel is non-zero
        abs       : absolute value
        total     : Sums up values in an image or in image1*image2 (img2 is the probability mask)
        mean      :  Average of values in an image or in image1*image2 (img2 is the probability mask)
        vtotal    : Sums up volumetrically weighted values in an image or in image1*image2 (img2 is the probability mask)
        Decision  : Computes result=1./(1.+exp(-1.0*( pix1-0.25)/pix2))
        Neg       : Produce image negative


    Parameters
    ----------
    volume1 : string
        nibabel-readable image volume
    volume2 : string
        nibabel-readable image volume
    operator : string
        ImageMath string corresponding to mathematical operator
    output_file : string
        nibabel-readable image volume

    Returns
    -------
    output_file : string
        name of output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import ImageMath
    >>> from mindboggle.utils.plots import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> volume1 = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> volume2 = os.path.join(path, 'arno', 'mri', 'mask.nii.gz')
    >>> operator = 'm'
    >>> output_file = ''
    >>> output_file = ImageMath(volume1, volume2, operator, output_file)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os
    from mindboggle.utils.utils import execute

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   os.path.basename(volume1) + '_' +
                                   os.path.basename(volume2))

    cmd = ['ImageMath', '3', output_file, operator, volume1, volume2]
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file


def ThresholdImage(volume, output_file='', threshlo=1, threshhi=10000):
    """
    Use the ThresholdImage function in ANTS to threshold image volume::

    Usage: ThresholdImage ImageDimension ImageIn.ext outImage.ext
           threshlo threshhi <insideValue> <outsideValue>

    Parameters
    ----------
    volume : string
        nibabel-readable image volume
    output_file : string
        nibabel-readable image volume
    threshlo : integer
        lower threshold
    threshhi : integer
        upper threshold

    Returns
    -------
    output_file : string
        name of output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import ThresholdImage
    >>> from mindboggle.utils.plots import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> volume = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> output_file = ''
    >>> threshlo = 500
    >>> threshhi = 10000
    >>> output_file = ThresholdImage(volume, output_file, threshlo, threshhi)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os
    from mindboggle.utils.utils import execute

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'threshold_' + os.path.basename(volume))

    cmd = 'ThresholdImage 3 {0} {1} {2} {3}'.format(volume, output_file,
                                                    threshlo, threshhi)
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file


def ANTS(source, target, iterations='30x99x11', output_stem=''):
    """
    Use ANTs to register a source image volume to a target image volume.

    This program uses the ANTs SyN registration method.

    Parameters
    ----------
    source : string
        file name of source image volume
    target : string
        file name of target image volume
    iterations : string
        number of iterations ("0" for affine, "30x99x11" default)
    output_stem : string
        file name stem for output transform matrix

    Returns
    -------
    affine_transform : string
        file name for affine transform matrix
    nonlinear_transform : string
        file name for nonlinear transform nifti file
    nonlinear_inverse_transform : string
        file name for nonlinear inverse transform nifti file
    output_stem : string
        file name stem for output transform matrix

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import ANTS
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> source = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> target = os.path.join(path, 'atlases', 'MNI152_T1_1mm_brain.nii.gz')
    >>> iterations = "0"
    >>> output_stem = ""
    >>> #
    >>> ANTS(source, target, iterations, output_stem)

    """
    import os
    from mindboggle.utils.utils import execute

    if not output_stem:
        src = os.path.basename(source).split('.')[0]
        tgt = os.path.basename(target).split('.')[0]
        output_stem = os.path.join(os.getcwd(), src+'_to_'+tgt)

    cmd = ['ANTS', '3', '-m CC[' + target + ',' + source + ',1,2]',
            '-r Gauss[2,0]', '-t SyN[0.5] -i', iterations,
            '-o', output_stem, '--use-Histogram-Matching',
            '--number-of-affine-iterations 10000x10000x10000x10000x10000']
    execute(cmd, 'os')

    affine_transform = output_stem + 'Affine.txt'
    nonlinear_transform = output_stem + 'Warp.nii.gz'
    nonlinear_inverse_transform = output_stem + 'InverseWarp.nii.gz'

    if not os.path.exists(affine_transform):
        raise(IOError(affine_transform + " not found"))
    if not os.path.exists(nonlinear_transform):
        raise(IOError(nonlinear_transform + " not found"))
    if not os.path.exists(nonlinear_inverse_transform):
        raise(IOError(nonlinear_inverse_transform + " not found"))

    return affine_transform, nonlinear_transform,\
           nonlinear_inverse_transform, output_stem


def WarpImageMultiTransform(source, target, output='',
                            interp='--use-NN', xfm_stem='',
                            affine_transform='', nonlinear_transform='',
                            inverse=False, affine_only=False):
    """
    Use ANTs to transform a source image volume to a target image volume.

    This program uses the ANTs WarpImageMultiTransform function.

    Parameters
    ----------
    source : string
        file name of source image volume
    target : string
        file name of target (reference) image volume
    output : string
        file name of output image volume
    interp : string
        interpolation type ("--use-NN" for nearest neighbor)
    xfm_stem : string
        file name stem for output transform
    affine_transform : string
        file containing affine transform
    nonlinear_transform : string
        file containing nonlinear transform
    inverse : Boolean
        apply inverse transform?
    affine_only : Boolean
        apply only affine transform?

    Returns
    -------
    output : string
        output label file name

    """
    import os
    import sys
    from mindboggle.utils.utils import execute

    if xfm_stem:
        affine_transform = xfm_stem + 'Affine.txt'
        if inverse:
            nonlinear_transform = xfm_stem + 'InverseWarp.nii.gz'
        else:
            nonlinear_transform = xfm_stem + 'Warp.nii.gz'
    elif not affine_transform and not nonlinear_transform:
        sys.exit('Provide either xfm_stem or affine_transform and '
                 'nonlinear_transform.')

    if not output:
        output = os.path.join(os.getcwd(), 'WarpImageMultiTransform.nii.gz')

    if not os.path.exists(nonlinear_transform):
        affine_only = True

    if affine_only:
        if inverse:
            cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
                   target, interp, '-i', affine_transform]
        else:
            cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
                   target, interp, affine_transform]
    else:
        if inverse:
            cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
                   target, interp, '-i', affine_transform, nonlinear_transform]
        else:
            cmd = ['WarpImageMultiTransform', '3', source, output, '-R',
                   target, interp, nonlinear_transform, affine_transform]
    execute(cmd, 'os')

    if not os.path.exists(output):
        raise(IOError(output + " not found"))

    return output


def PropagateLabelsThroughMask(mask, labels, mask_index=None,
                               output_file='', binarize=True, stopvalue=''):
    """
    Use ANTs to fill a binary volume mask with initial labels.

    This program uses ThresholdImage and the ImageMath
    PropagateLabelsThroughMask functions in ANTS.

    ThresholdImage ImageDimension ImageIn.ext outImage.ext
        threshlo threshhi <insideValue> <outsideValue>

    PropagateLabelsThroughMask: Final output is the propagated label image.
        ImageMath ImageDimension Out.ext PropagateLabelsThroughMask
        speed/binaryimagemask.nii.gz initiallabelimage.nii.gz ...

    Parameters
    ----------
    mask : string
        nibabel-readable image volume
    labels : string
        nibabel-readable image volume with integer labels
    mask_index : integer (optional)
        mask with just voxels having this value
    output_file : string
        nibabel-readable labeled image volume
    binarize : Boolean
        binarize mask?
    stopvalue : integer
        stopping value

    Returns
    -------
    output_file : string
        name of labeled output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import PropagateLabelsThroughMask
    >>> from mindboggle.utils.plots import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> labels = os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> mask = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> mask_index = None
    >>> output_file = ''
    >>> binarize = True
    >>> stopvalue = None
    >>> output_file = PropagateLabelsThroughMask(mask, labels, mask_index, output_file, binarize, stopvalue)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os
    from mindboggle.utils.utils import execute

    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'PropagateLabelsThroughMask.nii.gz')
        output_file = os.path.join(os.getcwd(),
                                   os.path.basename(labels) + '_through_' +
                                   os.path.basename(mask))

    print('mask: {0}, labels: {1}'.format(mask, labels))

    # Binarize image volume:
    if binarize:
        temp_file = os.path.join(os.getcwd(),
                                 'PropagateLabelsThroughMask.nii.gz')
        cmd = ['ThresholdImage', '3', mask, temp_file, '0 1 0 1']
        execute(cmd, 'os')
        mask = temp_file

    # Mask with just voxels having mask_index value:
    if mask_index:
        mask2 = os.path.join(os.getcwd(), 'temp.nii.gz')
        cmd = 'ThresholdImage 3 {0} {1} {2} {3} 1 0'.format(mask, mask2,
               mask_index, mask_index)
        execute(cmd)
    else:
        mask2 = mask

    # Propagate labels:
    cmd = ['ImageMath', '3', output_file, 'PropagateLabelsThroughMask',
            mask2, labels]
    if stopvalue:
        cmd.extend(stopvalue)
    execute(cmd, 'os')
    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file


def fill_volume_with_surface_labels(hemi, left_mask, right_mask,
                                    surface_files, mask_index=None,
                                    output_file='', binarize=False):
    """
    Use ANTs to fill a volume mask with surface mesh labels.

    Note ::

        - This uses PropagateLabelsThroughMask in the ANTS ImageMath function.

        - Partial volume information is lost when mapping surface to volume.

    Parameters
    ----------
    hemi : string
        either 'lh' or 'rh', indicating brain's left or right hemisphere
    left_mask : string
        nibabel-readable left brain image mask volume
    right_mask : string
        nibabel-readable left brain image mask volume
    surface_files : string or list of strings
        VTK file(s) containing surface mesh(es) with labels as scalars
    mask_index : integer (optional)
        mask with just voxels having this value
    output_file : string
        name of output file
    binarize : Boolean
        binarize mask?

    Returns
    -------
    output_file : string
        name of labeled output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import fill_volume_with_surface_labels
    >>> from mindboggle.utils.plots import plot_volumes
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> surface_files = [os.path.join(path, 'arno', 'labels',
    >>>     'lh.labels.DKT25.manual.vtk'), os.path.join(path, 'arno', 'labels',
    >>>     'rh.labels.DKT25.manual.vtk')]
    >>> # For a quick test, simply mask with whole brain:
    >>> hemi = 'rh'
    >>> left_mask = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> right_mask = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> mask_index = None
    >>> output_file = ''
    >>> binarize = True
    >>> output_file = fill_volume_with_surface_labels(hemi, left_mask, right_mask, surface_files, mask_index, output_file, binarize)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os

    from mindboggle.utils.io_vtk import transform_to_volume
    from mindboggle.labels.relabel import overwrite_volume_labels
    from mindboggle.utils.ants import PropagateLabelsThroughMask
    from mindboggle.utils.utils import execute

    if isinstance(surface_files, str):
        surface_files = [surface_files]

    if hemi == 'lh':
        mask = left_mask
    elif hemi == 'rh':
        mask = right_mask

    # Transform vtk coordinates to voxel index coordinates in a target
    # volume by using the header transformation:
    surface_in_volume = transform_to_volume(surface_files[0], mask)

    # Do the same for additional vtk surfaces:
    if len(surface_files) == 2:
        surfaces_in_volume = os.path.join(os.getcwd(), 'surfaces.nii.gz')
        surface_in_volume2 = transform_to_volume(surface_files[1], mask)

        overwrite_volume_labels(surface_in_volume, surface_in_volume2,
                                surfaces_in_volume, ignore_labels=[0])
        surface_in_volume = surfaces_in_volume

    # Use ANTs to fill a binary volume mask with initial labels:
    output_file = PropagateLabelsThroughMask(mask, surface_in_volume,
                                             mask_index, output_file,
                                             binarize)
    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file  # surface_in_volume


def thickinthehead(segmented_file, labeled_file, cortex_value=2,
                   noncortex_value=3, labels=[], resize=True, propagate=True,
                   output_dir='', output_table=''):
    """
    Compute a simple thickness measure for each labeled cortex region.

    Note::

      - Cortex, noncortex, and labeled files are the same coregistered brain.
      - Calls ANTs functions: ImageMath, Threshold, ResampleImageBySpacing

    Example preprocessing steps ::

      1. Run Freesurfer and antsCorticalThickness.sh on T1-weighted image.
      2. Convert FreeSurfer volume labels (e.g., wmparc.mgz or aparc+aseg.mgz)
         to cortex (2) and noncortex (3) segments using relabel_volume()
         function [refer to LABELS.py or FreeSurferColorLUT labels file].
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
    resize : Boolean
        resize (2x) segmented_file for more accurate thickness estimates?
    propagate : Boolean
        propagate labels through cortex?
    output_dir : string
        output directory
    output_table : string
        output table file with labels and thickness values
        (if empty, don't save table file)

    Returns
    -------
    label_volume_area51_thickness : list of lists of integers and floats
        label indices, volumes, areas, and thickness values
    table_file : string
        name of output table file

    Examples
    --------
    >>> from mindboggle.utils.ants import thickinthehead
    >>> segmented_file = '/Users/arno/Data/antsCorticalThickness/OASIS-TRT-20-1/tmp23314/tmpBrainSegmentation.nii.gz'
    >>> labeled_file = '/appsdir/freesurfer/subjects/OASIS-TRT-20-1/mri/labels.DKT31.manual.nii.gz'
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
    >>> resize = True
    >>> propagate = False
    >>> output_dir = ''
    >>> output_table = ''
    >>> label_volume_area51_thickness, table_file = thickinthehead(segmented_file, labeled_file, cortex_value, noncortex_value, labels, resize, propagate, output_dir, output_table)

    """
    import os
    import numpy as np
    import nibabel as nb

    from mindboggle.utils.utils import execute

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

    #-------------------------------------------------------------------------
    # Extract noncortex and cortex:
    #-------------------------------------------------------------------------
    cmd = ['ThresholdImage 3', segmented_file,
           noncortex, str(noncortex_value), str(noncortex_value), '1 0']
    execute(cmd)
    cmd = ['ThresholdImage 3', segmented_file,
           cortex, str(cortex_value), str(cortex_value), '1 0']
    execute(cmd)

    #-------------------------------------------------------------------------
    # Either mask labels with cortex or fill cortex with labels:
    #-------------------------------------------------------------------------
    if propagate:
        cmd = ['ImageMath', '3', cortex, 'PropagateLabelsThroughMask',
               cortex, labeled_file]
        execute(cmd)
    else:
        cmd = ['ImageMath 3', cortex, 'm', cortex, labeled_file]
        execute(cmd)

    #-------------------------------------------------------------------------
    # Resample cortex and noncortex files from 1x1x1 to 0.5x0.5x0.5
    # to better represent the contours of the boundaries of the cortex:
    #-------------------------------------------------------------------------
    if resize:
        rescale = 2.0

        dims = ' '.join([str(1/rescale), str(1/rescale), str(1/rescale)])
        cmd = ['ResampleImageBySpacing 3', cortex, cortex, dims, '0 0 1']
        execute(cmd)
        cmd = ['ResampleImageBySpacing 3', noncortex, noncortex, dims, '0 0 1']
        execute(cmd)

    #-------------------------------------------------------------------------
    # Extract outer and inner boundary voxels of the cortex,
    # by eroding 1 (resampled) voxel for cortex voxels (2) bordering
    # the outside of the brain (0) and bordering noncortex (3):
    #-------------------------------------------------------------------------
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
    label_volume_area51_thickness = np.zeros((len(labels), 4))
    label_volume_area51_thickness[:,0] = labels
    for ilabel, label in enumerate(labels):

        #---------------------------------------------------------------------
        # Compute thickness as a ratio of label volume and edge volume:
        #   - Estimate middle cortical surface area by the average volume
        #     of the outer and inner boundary voxels of the cortex.
        #   - Compute the volume of a labeled region of cortex.
        #   - Estimate the thickness of the labeled cortical region as the
        #     volume of the labeled region divided by the surface area.
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
            label_volume_area51_thickness[ilabel,1] = label_cortex_volume
            label_volume_area51_thickness[ilabel,2] = label_area
            label_volume_area51_thickness[ilabel,3] = thickness
            print('label {0} volume: cortex={1:2.2f}, inner={2:2.2f}, '
                  'outer={3:2.2f}, area51={4:2.2f}, thickness={5:2.2f}mm'.
                  format(label, label_cortex_volume, label_inner_edge_volume,
                  label_outer_edge_volume, label_area, thickness))

    if output_table:
        table_file = os.path.join(output_dir, output_table)
        np.savetxt(table_file, label_volume_area51_thickness,
                   fmt='%d %2.4f %2.4f %2.4f', delimiter='\t', newline='\n')
    else:
        table_file = ''

    label_volume_area51_thickness = label_volume_area51_thickness.\
        transpose().tolist()

    return label_volume_area51_thickness, table_file
