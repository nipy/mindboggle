#!/usr/bin/python

"""
Save the vertices of a FreeSurfer surface mesh as an image volume.


Authors:  
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def polydata2volume(surface_file, volume_file, output_file,
                    use_freesurfer_surfaces):
    """
    Save the vertices of a FreeSurfer surface mesh as an image volume.
    """

    from os import path
    import numpy as np
    import nibabel as nb
    import pyvtk

    if use_freesurfer_surfaces:
        trans = 128  # translation to middle of FreeSurfer conformed space

    # Check type:
    if type(surface_file) == str:
        pass
    elif type(surface_file) == list:
        surface_file = surface_file[0]
    else:
        import sys
        sys.error("Check format of " + surface_file)

    # Check type:
    if type(volume_file) == str:
        pass
    elif type(volume_file) == list:
        volume_file = volume_file[0]
    else:
        import sys
        sys.error("Check format of " + volume_file)

    # Load image volume
    vol = nb.load(volume_file)
    vol_shape = vol.shape
    xfm = vol.get_affine()

    # Load surface
    VTKReader = pyvtk.VtkData(surface_file)
    xyz =  VTKReader.structure.points
    xyz = np.array(xyz)

    # Create a new volume (permuted and flipped)
    V = np.zeros(vol_shape)
    for vertex in xyz:
        if use_freesurfer_surfaces:
            V[-vertex[0]+trans, -vertex[2]+trans, vertex[1]+trans] = 1
        else:
            V[vertex[0], vertex[1], vertex[2]] = 1
        
    # Save the image with the same affine transform
    img = nb.Nifti1Image(V, xfm)
    img.to_filename(output_file)


    """
    # Create a new volume (permuted and flipped)
    from apply_utils import apply_affine
    xfm2 = np.array([[-1,0,0,128],
                    [0,0,-1,128],
                    [0,1,0,128],
                    [0,0,0,1]],dtype=float)
    xyz = apply_affine(xyz[:,0], xyz[:,1], xyz[:,2], xfm2)

    V = np.zeros(vol_shape)
    for vertex in xyz:
        V[vertex[0], vertex[1], vertex[2]] = 1
    """

def label_volume(output_file, mask_file, input_file):
    """
    Fill (e.g., gray matter) volume with surface labels using ANTS
    (ImageMath's PropagateLabelsThroughMask)

    Brian avants: 
    The initial box labels are propagated through the gray matter with
    gm-probability dependent speed. It uses the fast marching algorithm.
    You can control how tightly the propagation follows the gray matter
    label by adjusting the speed image -- e.g. a binary speed image
    will constrain the propagated label only to the gm.

    """
    from os import path, getcwd
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    print("Fill gray matter volume with surface labels using ANTS...")

    output_file = path.join(getcwd(), output_file)

    # Check type:
    if type(mask_file) == str:
        pass
    elif type(mask_file) == list:
        mask_file = mask_file[0]
    else:
        import sys
        sys.error("Check format of " + mask_file)

    args = ['3',
            output_file,
            'PropagateLabelsThroughMask',
            mask_file,
            input_file]

    cli = CommandLine(command='ImageMath')
    cli.inputs.args = ' '.join(args)
    logger.info(cli.cmdline)
    cli.run()

    return output_file

"""
NB: To fill gray matter with labels using FreeSurfer,
    we would need to save the labels as an .annot file 
    in the subject directory (annot_name).

def label_volume_FS(subject_id, annot_name, output_name):
    ""
    Propagate surface labels through a gray matter volume 
    using FreeSurfer's mri_aparc2aseg
    ""
    from os import path, getcwd
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    print("Fill gray matter volume with surface labels using FreeSurfer...")

    output_file = path.join(getcwd(), output_name)

    args = ['--s', subject_id,
            '--annot', annot_name,
            '--o', output_file]

    cli = CommandLine(command='mri_aparc2aseg')
    cli.inputs.args = ' '.join(args)
    logger.info(cli.cmdline)
    cli.run()

    return output_file
"""

