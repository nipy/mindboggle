#!/usr/bin/env python
"""
Functions for converting to/from FreeSurfer formats.


Authors:
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Oliver Hinds, 2013  (ohinds@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def labels_to_annot(hemi, subjects_path, subject, label_files,
                    colortable, annot_name):
    """
    Convert FreeSurfer .label files to a FreeSurfer .annot file
    using FreeSurfer's mris_label2annot:
    https://surfer.nmr.mgh.harvard.edu/fswiki/mris_label2annot

    The order of the .label files must equal the order
    of the labels in the colortable:
    https://surfer.nmr.mgh.harvard.edu/fswiki/LabelsClutsAnnotationFiles

    NOTE:  The resulting .annot file will have incorrect labels
           if the numbering of the labels is not sequential from 1,2,3...
           For programs like tksurfer, the labels are identified
           by the colortable's RGB values, so to some programs that display
           the label names, the labels could appear correct when not.
    NOTE:  You cannot overwrite a .annot file of the same name,
           so in this script I delete it before running.

    Parameters
    ----------
    hemi :  hemisphere [string]
    subjects_path :  path to file
    subject :  subject name
    label_files :  .label file names [list of strings]
    colortable :  file of label numbers & names (same order as label_files)
    annot_name :  name of the output .annot file (without prepending hemi)

    Returns
    -------
    annot_name :  name of .annot file (without prepend)
    annot_file :  name of .annot file (with prepend)

    """
    import os
    from nipype.interfaces.base import CommandLine

    label_files = [f for f in label_files if f!=None]
    if label_files:
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

def thickness_to_ascii(hemi, subject, subjects_path):
    """
    Convert a FreeSurfer thickness (per-vertex) file
    to an ascii file.

    Note:  Untested function

    Parameters
    ----------
    hemi : string indicating left or right hemisphere
    subject_path: string
        path to subject directory where the binary FreeSurfer
        thickness file is found ("lh.thickness")

    Returns
    -------
    thickness_file : string
        name of output file, where each element is the thickness
        value of a FreeSurfer mesh vertex. Elements are ordered
        by orders of vertices in FreeSurfer surface file.

    """
    import os
    from nipype.interfaces.base import CommandLine

    filename = hemi + 'thickness'
    filename_full = os.path.join(subjects_path, subject, filename)
    thickness_file = os.path.join(os.getcwd(), filename + '.dat')

    cli = CommandLine(command='mri_convert')
    cli.inputs.args = ' '.join([filename_full, '--ascii+crsf', thickness_file])
    cli.cmdline
    cli.run()

    return thickness_file

#=============================================================================
# Functions for converting FreeSurfer files to/from VTK format
#=============================================================================
def surface_to_vtk(surface_file):
    """
    Convert FreeSurfer surface file to VTK format.

    If a file named orig.mgz exists in '../mri', the surface coordinates
    are transformed into scanner RAS space during format conversion
    according to the vox2ras transform in that file.

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_free import surface_to_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial')
    >>> #
    >>> surface_to_vtk(surface_file)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.plots import plot_vtk
    >>> plot_vtk('lh.pial.vtk')

    """
    import os
    import nibabel as nb

    from mindboggle.utils.io_vtk import write_header, write_points, write_faces

    surf = nb.freesurfer.read_geometry(surface_file)
    points = surf[0]
    faces = surf[1]

    # Transform surface coordinates into normal scanner RAS.
    # See example 3 in "Transforms within a subject's anatomical space":
    # https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
    orig_file = os.path.join(os.path.dirname(surface_file),
                             "..", "mri", "orig.mgz")

    if os.path.exists(orig_file):
        import numpy as np
        Norig = nb.load(orig_file).get_affine()
        Torig = np.array([[-1, 0, 0, 128],
                          [0, 0, 1, -128],
                          [0, -1, 0, 128],
                          [0, 0, 0, 1]], dtype=float)
        xfm = np.dot(Norig, np.linalg.inv(Torig))
        points = np.transpose(np.dot(xfm, np.transpose(
            np.concatenate((points, np.ones((np.shape(points)[0],1))),
                           axis=1))))[:,0:3]

    output_vtk = os.path.join(os.getcwd(),
                              os.path.basename(surface_file + '.vtk'))
    Fp = open(output_vtk, 'w')
    write_header(Fp, Title='vtk output from ' + surface_file)
    write_points(Fp, points)
    write_faces(Fp, faces)
    Fp.close()

    return output_vtk

def curvature_to_vtk(surface_file, vtk_file):
    """
    Convert FreeSurfer curvature, thickness, or convexity file to VTK format.

    Parameters
    ----------
    surface_file : string  (name of FreeSurfer surface file)
    vtk_file : string  (name of VTK surface file)

    Returns
    -------
    output_vtk : string
        name of output VTK file, where each vertex is assigned
        the corresponding shape value.

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import curvature_to_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(path, 'arno', 'freesurfer', 'lh.thickness')
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> #
    >>> curvature_to_vtk(surface_file, vtk_file)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.plots import plot_vtk
    >>> plot_vtk('lh.thickness.vtk')

    """
    import os
    import nibabel as nb

    from mindboggle.utils.io_vtk import rewrite_scalars

    output_vtk = os.path.join(os.getcwd(), os.path.basename(surface_file)+'.vtk')

    curvature_values = nb.freesurfer.read_morph_data(surface_file)
    scalar_names = os.path.basename(surface_file)

    rewrite_scalars(vtk_file, output_vtk, curvature_values, scalar_names)

    return output_vtk

def annot_to_vtk(annot_file, vtk_file):
    """
    Load a FreeSurfer .annot file and save as a VTK format file.

    Parameters
    ----------
    annot_file : string
        name of FreeSurfer .annot file
    vtk_file : string
        name of VTK surface file

    Returns
    -------
    labels : list
        integers (one label per vertex)
    output_vtk : string
        name of output VTK file, where each vertex is assigned
        the corresponding shape value.

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import annot_to_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> annot_file = os.path.join(path, 'arno', 'freesurfer', 'lh.aparc.annot')
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> #
    >>> labels, output_vtk = annot_to_vtk(annot_file, vtk_file)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.plots import plot_vtk
    >>> plot_vtk('lh.aparc.vtk')

    """
    import os
    import nibabel as nb
    from mindboggle.utils.io_vtk import rewrite_scalars

    labels, colortable, names = nb.freesurfer.read_annot(annot_file)

    output_vtk = os.path.join(os.getcwd(),
                              os.path.basename(annot_file).strip('.annot') + '.vtk')

    rewrite_scalars(vtk_file, output_vtk, labels, 'Labels')

    return labels, output_vtk

def vtk_to_labels(hemi, surface_file, label_numbers, label_names,
                  RGBs, scalar_name):
    """
    Write FreeSurfer .label files from a labeled VTK surface mesh.

    From https://surfer.nmr.mgh.harvard.edu/fswiki/LabelsClutsAnnotationFiles:

        "A label file is a text file capturing a list of vertices belonging to a region,
        including their spatial positions(using R,A,S coordinates). A label file
        corresponds only to a single label, thus contains only a single list of vertices"::

            1806
            7  -22.796  -66.405  -29.582 0.000000
            89  -22.273  -43.118  -24.069 0.000000
            138  -14.142  -81.495  -30.903 0.000000
            [...]

    Parameters
    ----------
    hemi :  string
        hemisphere
    surface_file :  string
        vtk surface mesh file with labels
    label_numbers :  list of integers
        label numbers
    label_names :  list of strings
        label names
    RGBs :  list of lists of 3-tuples
        label RGB values for later conversion to a .annot file
    scalar_name :  string
        name of scalar values in vtk file

    Returns
    -------
    label_files :  list of strings
        label file names (order must match label list)
    colortable :  string
        file with list of labels and RGB values
        NOTE: labels are identified by the colortable's RGB values

    """
    import os
    import numpy as np
    import vtk

    def string_vs_list_check(var):
        """
        Check type to make sure it is a string.

        (if a list, return the first element)
        """

        # Check type:
        if type(var) == str:
            return var
        elif type(var) == list:
            return var[0]
        else:
            os.error("Check format of " + var)

    # Check type to make sure the filename is a string
    # (if a list, return the first element)
    surface_file = string_vs_list_check(surface_file)

    # Initialize list of label files and output colortable file
    label_files = []
    #relabel_file = os.path.join(os.getcwd(), 'relabel_annot.txt')
    #f_relabel = open(relabel_file, 'w')
    colortable = os.path.join(os.getcwd(), 'colortable.ctab')
    f_rgb = open(colortable, 'w')

    # Loop through labels
    irgb = 0
    for ilabel, label_number in enumerate(label_numbers):

        # Check type to make sure the number is an int
        label_number = int(label_number)
        label_name = label_names[ilabel]

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
            irgb += 1

            # Write to relabel_file
            #if irgb != label_number:
            #    f_relabel.writelines('{0} {1}\n'.format(irgb, label_number))

            # Write to colortable
            f_rgb.writelines('{0} {1} {2}\n'.format(
                irgb, label_name, "0 0 0 0")) # ".join(RGBs[ilabel])))

            # Store in list of .label files
            label_file = hemi + '.' + label_name + '.label'
            label_file = os.path.join(os.getcwd(), label_file)
            label_files.append(label_file)

            # Write to .label file
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
    f_rgb.close()
    #f_relabel.close()

    return label_files, colortable
