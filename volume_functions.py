#!/usr/bin/python

"""
Save the vertices of a FreeSurfer surface mesh as an image volume.


Authors:  
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""
import os

def write_label_file(hemi, surface_file, label_index, label_name):
    """
    Save label file for a given label from the vertices of a labeled VTK surface mesh.

    """

    from os import path, getcwd, error
    import numpy as np
    import vtk

    scalar_name = "Max_(majority_labels)"

    # Check type:
    if type(surface_file) == str:
        pass
    elif type(surface_file) == list:
        surface_file = surface_file[0]
    else:
        error("Check format of " + surface_file)

    # Load surface
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(surface_file)
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    d = data.GetPointData()
    labels = d.GetArray(scalar_name)

    # Write vertex index, coordinates, and 0
    count = 0
    npoints = data.GetNumberOfPoints()
    L = np.zeros((npoints,5))
    for i in range(npoints):
        label = labels.GetValue(i)
        if label == label_index:
            L[count,0] = i
            L[count,1:4] = data.GetPoint(i)
            count += 1
 
    # Save the label file
    if count > 0:
        label_file = path.join(getcwd(), hemi + '.' + label_name + '.label')
        f = open(label_file, 'w')
        f.writelines('#!ascii label\n' + str(count) + '\n')
        for i in range(npoints):
            if any(L[i,:]):
                printline = '{0} {1} {2} {3} 0\n'.format(
                             np.int(L[i,0]), L[i,1], L[i,2], L[i,3])
                f.writelines(printline)
            else:
                break
        f.close()
        return label_file

def label_to_annot_file(hemi, subjects_path, subject, label_files, lookup_table):
    """
    Save label file for a given label from the vertices of a labeled VTK surface mesh.
    """

    from os import path
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')
    
    label_files = [f for f in label_files if f!=None]
    if label_files:
        cli = CommandLine(command='mris_label2annot')
        annot_name = 'labels.max'
        annot_file = hemi + '.' + annot_name + '.annot'
        if path.exists(path.join(subjects_path, subject, 'label', annot_file)):
            pass
        else:
            cli.inputs.args = ' '.join(['--hemi', hemi, '--s', subject, \
                                        '--l', ' --l '.join(label_files), \
                                        '--ctab', lookup_table, \
                                        '--a', annot_name])
            logger.info(cli.cmdline)
            cli.run()
        return annot_name, annot_file

def fill_label_volume(subject, annot_name):
    """
    Propagate surface labels through a gray matter volume 
    using FreeSurfer's mri_aparc2aseg
    """

    from os import path, getcwd
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    print("Fill gray matter volume with surface labels using FreeSurfer...")

    output_file = path.join(getcwd(), annot_name + '.nii.gz')

    args = ['--s', subject,
            '--annot', annot_name,
            '--o', output_file]

    cli = CommandLine(command='mri_aparc2aseg')
    cli.inputs.args = ' '.join(args)
    logger.info(cli.cmdline)
    cli.run()

    return output_file

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

