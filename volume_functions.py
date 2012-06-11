#!/usr/bin/python

"""
Save the vertices of a FreeSurfer surface mesh as an image volume.


Authors:  
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def surface_to_volume(surface_file, volume_file, use_freesurfer):
    """
    Save the vertices of a FreeSurfer surface mesh as an image volume.
    """

    from os import path, getcwd, error
    import numpy as np
    import nibabel as nb
    import vtk

    scalar_name = "Max_(majority_labels)"
    if use_freesurfer:
        trans = 128  # translation to middle of FreeSurfer conformed space

    # Check type:
    if type(surface_file) == str:
        pass
    elif type(surface_file) == list:
        surface_file = surface_file[0]
    else:
        error("Check format of " + surface_file)

    # Check type:
    if type(volume_file) == str:
        pass
    elif type(volume_file) == list:
        volume_file = volume_file[0]
    else:
        error("Check format of " + volume_file)

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
    output_file = path.join(getcwd(), surface_file.strip('.nii.gz')+'.nii.gz')
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

    # Check type:
    if type(mask_file) == str:
        pass
    elif type(mask_file) == list:
        mask_file = mask_file[0]
    else:
        error("Check format of " + mask_file)

    # Check type:
    if type(input_file) == str:
        pass
    elif type(input_file) == list:
        input_file = input_file[0]
    else:
        error("Check format of " + input_file)

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

"""
NB: To fill gray matter with labels using FreeSurfer,
    we would need to save the labels as an .annot file 
    in the subject directory (annot_name).

def maxlabel_volume_FS(subject, annot_name, output_name):
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

    args = ['--s', subject,
            '--annot', annot_name,
            '--o', output_file]

    cli = CommandLine(command='mri_aparc2aseg')
    cli.inputs.args = ' '.join(args)
    logger.info(cli.cmdline)
    cli.run()

    return output_file
"""

def measure_volume_overlap(subject, labels, input_file, atlases_path, atlases, atlases2):
    """
    Measure overlap between individual label regions in a source and target image.

    Input:
    arg1, arg2:  source, target images, consisting of index-labeled pixels/voxels
    arg3:  list of label indices

    """

    from os import path, getcwd, error
    import nibabel as nb
    import numpy as np
    from nipype import logging
    logger = logging.getLogger('interface')

    # Input files
    if type(input_file) == str:
        pass
    elif type(input_file) == list:
        input_file = input_file[0]
    else:
        error("Check format of " + input_file)

    # Find atlas in atlases list that corresponds to subject (in atlases2 list)
    if subject in atlases2:
        iatlas = atlases2.index(subject)
        atlas = atlases[iatlas]
        atlas_file = path.join(atlases_path, 'atlases', atlas, 'aparcNMMjt+aseg.mgz')
    else:
        import sys
        sys.exit(subject + " not in list")
    input_data = nb.load(input_file).get_data().ravel()
    atlas_data = nb.load(atlas_file).get_data().ravel()

    # Set up the output csv file
    output_table = path.join(getcwd(), input_file.strip('.nii.gz')+'.csv')
    try:
        f = open(output_table,"w")
    except IOError:
        raise

    # For each label, compute Dice and Jaccard coefficients
    avg_dice = 0
    avg_jacc = 0
    f.writelines("Label, Dice, Jaccard\n")
    for label in labels:
        f.writelines(str(label)+", ")
        input_indices = np.where(input_data==label)[0]
        atlas_indices = np.where(atlas_data==label)[0]
        input_label_sum = np.sum(input_indices)
        atlas_label_sum = np.sum(atlas_indices)
        intersect_label_sum = np.sum(np.intersect1d(input_indices, atlas_indices))
        union_label_sum = np.sum(np.union1d(input_indices, atlas_indices))
        #print(str(label), str(input_label_sum), str(atlas_label_sum),
        #      str(intersect_label_sum), str(union_label_sum))

        if intersect_label_sum > 0:
            dice = 2 * intersect_label_sum / (input_label_sum + atlas_label_sum)
            jacc = intersect_label_sum / union_label_sum
            f.writelines(str(dice) + ", " + str(jacc) + "\n")
            avg_dice += dice
            avg_jacc += jacc
        else:
            f.writelines("0, 0\n")
    avg_dice = avg_dice/len(labels)
    avg_jacc = avg_jacc/len(labels)

    f.writelines("Total, " + str(avg_dice) + ", " + str(avg_jacc) + "\n")
    f.close()

    print('Average Dice: ' + str(avg_dice))
    print('Average Jacc: ' + str(avg_jacc))

    return output_table

