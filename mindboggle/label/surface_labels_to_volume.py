#!/usr/bin/python

"""
Convert surface mesh labels to volume labels and evaluate.

1. Write surface mesh labels FreeSurfer's .label file.
2. Use FreeSurfer's mris_label2annot and mri_aparc2aseg
   to convert these label files to .annot files and fill
   a gray matter volume with the labels.
3. Measure volume overlap between the labels of two volumes.

Authors:
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def write_label_file(hemi, surface_file, label_number, label_name):
    """
    Save a FreeSurfer .label file for a given label from the vertices
    of a labeled VTK surface mesh.

    """

    import os
    import numpy as np
    import io_file
    import vtk

    scalar_name = "Labels"

    # Check type to make sure the filename is a string
    # (if a list, return the first element)
    surface_file = io_file.string_vs_list_check(surface_file)

    # Check type to make sure the number is an int
    label_number = int(label_number)

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
        if label == label_number:
            L[count,0] = i
            L[count,1:4] = data.GetPoint(i)
            count += 1

    # Save the label file
    if count > 0:
        label_file = os.path.join(os.getcwd(),
                                  hemi + '.' + label_name + '.label')
        f = open(label_file, 'w')
        f.writelines('#!ascii label\n' + str(count) + '\n')
        for i in range(npoints):
            if any(L[i,:]):
                pr = '{0} {1} {2} {3} 0\n'.format(
                     np.int(L[i,0]), L[i,1], L[i,2], L[i,3])
                f.writelines(pr)
            else:
                break
        f.close()

        return label_file

def label_to_annot_file(hemi, subjects_path, subject, label_files, colortable):
    """
    Convert FreeSurfer .label files as a FreeSurfer .annot file
    using FreeSurfer's mris_label2annot.

    """

    import os
    from nipype.interfaces.base import CommandLine

    label_files = [f for f in label_files if f!=None]
    if label_files:
        annot_name = 'labels.max'
        annot_file = hemi + '.' + annot_name + '.annot'
        if os.path.exists(os.path.join(subjects_path, subject, 'label', annot_file)):
            cli = CommandLine(command='rm')
            cli.inputs.args = os.path.join(subjects_path, subject, \
                                           'label', annot_file)
            cli.cmdline
            cli.run()
        cli = CommandLine(command='mris_label2annot')
        cli.inputs.args = ' '.join(['--h', hemi, '--s', subject, \
                                    '--l', ' --l '.join(label_files), \
                                    '--ctab', colortable, \
                                    '--a', annot_name])
        cli.cmdline
        cli.run()
        return annot_name, annot_file

def fill_label_volume(subject, annot_name):
    """
    Propagate surface labels through a gray matter volume
    using FreeSurfer's mri_aparc2aseg

    """

    import os
    from nipype.interfaces.base import CommandLine

    print("Fill gray matter volume with surface labels using FreeSurfer...")

    output_file = os.path.join(os.getcwd(), annot_name + '.nii.gz')

    args = ['--s', subject,
            '--annot', annot_name,
            '--o', output_file]

    cli = CommandLine(command='mri_aparc2aseg')
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    cli.run()

    return output_file
