
##############################################################################
#
#   Use ANTS ImageMath rather than FreeSurfer to fill a volume:
#
##############################################################################

# pipeline: 

    """
    #-------------------------------------------------------------------------
    # Put surface vertices in a volume
    #-------------------------------------------------------------------------
    surf2vol = node(name='Surface_to_volume',
                    interface = fn(function = surface_to_volume,
                                   input_names = ['surface_file',
                                                  'volume_file',
                                                  'use_freesurfer'],
                                   output_names = ['output_file']))
    surf2vol.inputs.use_freesurfer = use_freesurfer
    atlasflow.add_nodes([surf2vol])
    atlasflow.connect([(vote, surf2vol, [('maxlabel_file','surface_file')])])
    if use_freesurfer:
        mbflow.connect([(vol, atlasflow,
                         [('volume_file','Surface_to_volume.volume_file')])])
    else:
        mbflow.connect([(vol, atlasflow,
                         [('volume_file','Surface_to_volume.volume_file')])])
    #-------------------------------------------------------------------------
    # Fill volume mask with surface vertex labels
    #-------------------------------------------------------------------------
    fill_maxlabels = node(name='Fill_volume_maxlabels',
                          interface = fn(function = fill_volume,
                                         input_names = ['command',
                                                        'input_file',
                                                        'mask_file'],
                                         output_names = ['output_file']))
    fill_maxlabels.inputs.command = imagemath
    atlasflow.add_nodes([fill_maxlabels])
    atlasflow.connect([(surf2vol, fill_maxlabels,
                        [('output_file', 'input_file')])])
    if use_freesurfer:
        mbflow.connect([(convertvol, atlasflow,
                         [('out_file','Fill_volume_maxlabels.mask_file')])])
    else:
        mbflow.connect([(vol, atlasflow,
                         [('volume_file','Fill_volume_maxlabels.mask_file')])])
    mbflow.connect([(atlasflow, datasink,
                     [('Fill_volume_maxlabels.output_file',
                       'labels.@maxvolume')])])
    """

# volume_functions.py:

def surface_to_volume(surface_file, volume_file, use_freesurfer):
    """
    Save the vertices of a FreeSurfer surface mesh as an image volume.

    FIX:  labels are incorrect!
    """

    from os import path, getcwd, error
    import numpy as np
    import nibabel as nb
    import vtk

    scalar_name = "Max_(majority_labels)"
    if use_freesurfer:
        trans = 128  # translation to middle of FreeSurfer conformed space

    # Load image volume
    vol = nb.load(volume_file)
    vol_shape = vol.shape
    xfm = vol.get_affine()

    # Load surface
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(surface_file)
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    d = data.GetPointData()
    labels = d.GetArray(scalar_name)

    # Create a new volume (permuted and flipped)
    V = np.zeros(vol_shape)
    npoints = data.GetNumberOfPoints()
    for i in range(npoints):
        point = data.GetPoint(i)
        label = labels.GetValue(i)
        if use_freesurfer:
            V[-point[0]+trans, -point[2]+trans, point[1]+trans] = label
        else:
            V[point[0], point[1], point[2]] = label

    # Save the image with the same affine transform
    output_file = path.join(getcwd(), surface_file.strip('.vtk')+'.nii.gz')
    img = nb.Nifti1Image(V, xfm)
    img.to_filename(output_file)

    return output_file

    """
    # Alternative (NEEDS A FIX):
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

def fill_volume(command, input_file, mask_file):
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
    from os import path, getcwd, error
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    print("Fill gray matter volume with surface labels using ANTS...")

    output_file = path.join(getcwd(), input_file.strip('.nii.gz')+'.fill.nii.gz')

    args = ['3',
            output_file,
            'PropagateLabelsThroughMask',
            mask_file,
            input_file]

    cli = CommandLine(command=command)
    cli.inputs.args = ' '.join(args)
    logger.info(cli.cmdline)
    cli.run()

    return output_file
