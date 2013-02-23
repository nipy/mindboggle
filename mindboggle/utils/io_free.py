#!/usr/bin/env python
"""
Functions for reading surfaces and converting between FreeSurfer formats.

This Python library reads and writes different file types.
The function read_surface() reads in surface files,
while the function read_curvature() reads in both
curvature (.curv) and convexity (.sulc) files.


Authors:
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def read_surface(filename):
    """
    Read in a FreeSurfer triangle surface mesh in binary format.

    Parameters
    ----------
    filename : string
        A binary FreeSurfer Triangle Surface file

    Returns
    -------
    Vertex : list of 3-tuples of floats
        Each element is a 3-tuple (list) of floats, which are the X, Y and Z coordinates of a vertex, respectively.
        A 3-tuple's index in the list *Vertex* is the ID of a vertex.

    Face : list of 3-tuples of integers
        Each element is a 3-tuple (list) of integers, which are the IDs of 3 vertices that form one face

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_free import read_surface
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(data_path, 'arno', 'freesurfer', 'lh.pial')
    >>> Vrtx, Face = read_surface(surface_file)
    >>> len(Vrtx)
      145069
    >>> len(Face)
      290134
    >>> Vrtx[10]
      [-16.585412979125977, -101.663330078125, 3.4243862628936768]
    >>> Face[10]
      [2, 50, 54]

    """
    import os
    import struct

    f = open(filename, "rb")
    f.seek(3)  # skip the first 3 Bytes "Magic" number

    s = f.read(500)   # the second field is a string of variable length
    End2 = s.find('\n\n',0)  # end of the second field is a '\n\n'

    f.seek(3+End2+2)  # jump to immediate Byte after the creating information

    s = f.read(8)
    VertexCount, FaceCount = struct.unpack(">ii", s)
    # print("This hemisphere has", VertexCount, "Vertexes and", FaceCount, "Faces")

    Vertex, Face = [], []

    for i in xrange(0, VertexCount):
        s = f.read(8)
        R, A = struct.unpack(">ff", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        A, S = struct.unpack(">ff", s)
        Vertex.append([R,A,S]) # R, A, S are the coordinates of vertices

    for i in xrange(0, FaceCount):
        s = f.read(8)
        V0, V1 = struct.unpack(">ii", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        V1, V2 = struct.unpack(">ii", s)
        Face.append([V0, V1, V2])

    return Vertex, Face

def read_curvature(filename):
    """
    Read in a FreeSurfer curvature, convexity, or thickness file.

    Parameters
    ----------
    filename : string
        A binary FreeSurfer (per-vertex) curvature file

    Returns
    -------
    Curvature : list of floats
        Each element is the curvature value of a FreeSurfer mesh vertex.

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_free import read_curvature
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(data_path, 'arno', 'freesurfer', 'lh.curv')
    >>> Curv = read_curvature(surface_file)
    >>> len(Curv)
      145069
    >>> Curv[10]
      -0.27393144369125366

    """
    import os
    import struct

    f = open(filename, "rb")

    f.seek(3) # skip the first 3 Bytes "Magic" number

    s = f.read(8)  # get the VertexCount and FaceCount
    VertexCount, FaceCount = struct.unpack(">ii", s)
    # print("# of Vertexes:", VertexCount, ", # of Faces:", FaceCount)

    Curvature = [0.0]

    s = f.read(8)
    ValsPerVertex, Curvature[0] = struct.unpack(">if", s)

    VertexCount -= 1  # because the first curvature value has been loaded

    while VertexCount > 1:
        s = f.read(8)
        VertexVal1, VertexVal2  =  struct.unpack(">ff", s)
        Curvature += [VertexVal1, VertexVal2]
        VertexCount -= 2

    if VertexCount != 0:  # number of vertices is even (NOT ODD!!!)
        f.seek(-4, os.SEEK_CUR)  # backward 4 Bytes from current position
        s = f.read(8)
        VertexVal1, VertexVal2 = struct.unpack(">ff", s)
        Curvature.append(VertexVal2)

    f.close()

    return Curvature

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

def labels_to_volume(subject, annot_name, original_space, reference):
    """
    Propagate surface labels through a gray matter volume
    and convert label volume from FreeSurfer 'unconformed' to original space
    using FreeSurfer's mri_aparc2aseg and mri_vol2vol.

    Note
    ----
    From the mri_aparc2aseg documentation:
    The volumes of the cortical labels will be different than
    reported by mris_anatomical_stats because partial volume information
    is lost when mapping the surface to the volume. The values reported by
    mris_anatomical_stats will be more accurate than the volumes from the
    aparc+aseg volume.

    Parameters
    ----------
    subject : string
        subject name
    annot_name: string
        FreeSurfer annot filename without the hemisphere prepend or .annot append
    original_space: Boolean
        convert from FreeSurfer unconformed to original space?
    reference : string
        file in original space

    Returns
    -------
    output_file : string
        name of output file

    """
    import os
    from nipype.interfaces.base import CommandLine

    # Fill gray matter volume with surface labels using FreeSurfer
    print("Fill gray matter volume with surface labels using FreeSurfer...")

    output_file1 = os.path.join(os.getcwd(), annot_name + '.nii.gz')

    args = ['--s', subject,
            '--annot', annot_name,
            '--o', output_file1]

    cli = CommandLine(command='mri_aparc2aseg')
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    cli.run()

    # Convert label volume from FreeSurfer to original space
    if original_space:
        print("Convert label volume from FreeSurfer 'unconformed' to original space...")

        output_file2 = os.path.join(os.getcwd(), annot_name + '.native.nii.gz')

        interp = 'nearest'
        args = ['--mov', output_file1,
                '--targ', reference,
                '--interp', interp,
                '--regheader --o', output_file2]

        cli = CommandLine(command='mri_vol2vol')
        cli.inputs.args = ' '.join(args)
        cli.cmdline
        cli.run()
        output_file = output_file2
    else:
        output_file = output_file1

    return output_file

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

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import surface_to_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial')
    >>> #
    >>> surface_to_vtk(surface_file)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('lh.pial.vtk')

    """
    import os
    from mindboggle.utils.io_free import read_surface
    from mindboggle.utils.io_vtk import write_header, write_points, write_faces

    points, faces = read_surface(surface_file)

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
    >>> vtk_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> #
    >>> curvature_to_vtk(surface_file, vtk_file)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('lh.thickness.vtk')

    """
    import os
    from mindboggle.utils.io_free import read_curvature
    from mindboggle.utils.io_vtk import rewrite_scalars

    output_vtk = os.path.join(os.getcwd(), os.path.basename(surface_file)+'.vtk')
    curvature_values = read_curvature(surface_file)
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
    >>> vtk_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> #
    >>> labels, output_vtk = annot_to_vtk(annot_file, vtk_file)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.mesh import plot_vtk
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
    hemi :  hemisphere [string]
    surface_file :  vtk surface mesh file with labels [string]
    label_numbers :  label numbers [list of strings]
    label_names :  label names [list of strings]
    RGBs :  list of label RGB values for later conversion to a .annot file
    scalar_name :  name of scalar values in vtk file [string]

    Returns
    -------
    label_files :  list of .label file names (order must match label list)
    colortable :  file with list of labels and RGB values
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
                irgb, label_name, RGBs[ilabel]))

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

    return label_files, colortable  #relabel_file

