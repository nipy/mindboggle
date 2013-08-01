#!/usr/bin/python
"""
Functions that call the legacy ANTs package created by PICSL at UPenn.

Authors:
Arno Klein, 2011-2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


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
    >>> from mindboggle.utils.ants import register_volume
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> source = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> target = os.path.join(path, 'atlases', 'MNI152_T1_1mm_brain.nii.gz')
    >>> iterations = "0"
    >>> output_stem = ""
    >>> #
    >>> ANTS(source, target, iterations, output_stem)

    """
    import os

    if not output_stem:
        src = os.path.basename(source).split('.')[0]
        tgt = os.path.basename(target).split('.')[0]
        output_stem = os.path.join(os.getcwd(), src+'_to_'+tgt)

    args = ['ANTS 3 -m', 'CC[' + target + ',' + source + ',1,2]',
            '-r Gauss[2,0] -t SyN[0.5] -i', iterations,
            '-o', output_stem, '--use-Histogram-Matching',
            '--number-of-affine-iterations 10000x10000x10000x10000x10000']
    cmd = ' '.join(args)
    print(cmd)
    os.system(cmd)  # p = Popen(args);

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

    if xfm_stem:
        affine_transform = xfm_stem + 'Affine.txt'
        if inverse:
            nonlinear_transform = xfm_stem + 'InverseWarp.nii.gz'
        else:
            nonlinear_transform = xfm_stem + 'Warp.nii.gz'
    else:
        xfm_stem = os.path.join(os.getcwd(),
                                os.path.basename(source).split('.')[0] +
                                '_to_' +
                                os.path.basename(target).split('.')[0])

    if not output:
        output = os.path.join(os.getcwd(), 'transformed.nii.gz')

    if not os.path.exists(nonlinear_transform):
        affine_only = True

    if affine_only:
        if inverse:
            args = ['WarpImageMultiTransform 3', source, output, '-R', target,
                    interp, '-i', affine_transform]
        else:
            args = ['WarpImageMultiTransform 3', source, output, '-R', target,
                    interp, affine_transform]
    else:
        if inverse:
            args = ['WarpImageMultiTransform 3', source, output, '-R', target,
                    interp, '-i', affine_transform, nonlinear_transform]
        else:
            args = ['WarpImageMultiTransform 3', source, output, '-R', target,
                    interp, nonlinear_transform, affine_transform]
    cmd = ' '.join(args)
    print(cmd)
    os.system(cmd)  # p = Popen(args);

    if not os.path.exists(output):
        raise(IOError(output + " not found"))

    return output


def PropagateLabelsThroughMask(mask_volume, label_volume, output_file='',
                               binarize=True):
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
    label_volume : string
        nibabel-readable image volume with integer labels
    mask_volume : string
        nibabel-readable image volume
    output_file : string
        nibabel-readable labeled image volume
    binarize : Boolean
        binarize mask_volume?

    Returns
    -------
    output_file : string
        name of labeled output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import PropagateLabelsThroughMask
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> label_volume = os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> mask_volume = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> output_file = ''
    >>> binarize = True
    >>> output_file = PropagateLabelsThroughMask(mask_volume, label_volume,
    >>>                                          output_file, binarize)
    >>> # View
    >>> plot_volumes(output_file)

    """
    import os

    if not output_file:
        output_file = os.path.join(os.getcwd(), 'propagated_labels.nii.gz')

    # Binarize image volume:
    if binarize:
        temp_file = os.path.join(os.getcwd(), 'propagated_labels.nii.gz')
        args = ['ThresholdImage 3', mask_volume, temp_file, '0 1 0 1']
        cmd = ' '.join(args)
        print(cmd)
        os.system(cmd)  # p = Popen(args);
        mask_volume = temp_file

    # Propagate labels:
    args = ['ImageMath 3', output_file, 'PropagateLabelsThroughMask',
            mask_volume, label_volume]
    cmd = ' '.join(args)
    print(cmd)
    os.system(cmd)  # p = Popen(args);

    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file


def fill_volume_with_surface_labels(volume_mask, surface_files,
                                    output_file='', binarize=False):
    """
    Use ANTs to fill a volume mask with surface mesh labels.

    This program uses PropagateLabelsThroughMask in ANTS's ImageMath function.

    Note ::
        Partial volume information is lost when mapping the surface
        to the volume.

    Parameters
    ----------
    volume_mask : string
        nibabel-readable image volume
    surface_files : string or list of strings
        VTK file(s) containing surface mesh(es) with labels as scalars
    output_file : string
        name of output file
    binarize : Boolean
        binarize volume_mask?

    Returns
    -------
    output_file : string
        name of labeled output nibabel-readable image volume

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.ants import fill_volume_with_surface_labels
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> surface_files = [os.path.join(path, 'arno', 'labels',
    >>>     'lh.labels.DKT25.manual.vtk'), os.path.join(path, 'arno', 'labels',
    >>>     'rh.labels.DKT25.manual.vtk')]
    >>> volume_mask = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> output_file = ''
    >>> binarize = True
    >>> fill_volume_with_surface_labels(volume_mask, surface_files,
    >>>                                 output_file, binarize)

    """
    import os

    from mindboggle.utils.io_vtk import transform_to_volume
    from mindboggle.labels.relabel import overwrite_volume_labels
    from mindboggle.utils.ants import PropagateLabelsThroughMask

    if isinstance(surface_files, str):
        surface_files = [surface_files]

    # Transform vtk coordinates to voxel index coordinates in a target
    # volume by using the header transformation:
    surface_in_volume = transform_to_volume(surface_files[0], volume_mask)

    # Do the same for additional vtk surfaces:
    if len(surface_files) == 2:
        surfaces_in_volume = os.path.join(os.getcwd(), 'surfaces.nii.gz')
        surface_in_volume2 = transform_to_volume(surface_files[1], volume_mask)

        overwrite_volume_labels(surface_in_volume, surface_in_volume2,
                                surfaces_in_volume, ignore_labels=[0])
        surface_in_volume = surfaces_in_volume

    # Use ANTs to fill a binary volume mask with initial labels:
    output_file = PropagateLabelsThroughMask(volume_mask, surface_in_volume,
                                             output_file, binarize)

    if not os.path.exists(output_file):
        raise(IOError(output_file + " not found"))

    return output_file  # surface_in_volume
