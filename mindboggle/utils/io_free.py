#!/usr/bin/python

"""
Functions related to reading and writing FreeSurfer files.

This Python library reads and writes different file types.
In particular, it has functions to read some FreeSurfer files,
including surface, curvature, and convexity files.
The function read_surface() reads in surface files,
while the function read_curvature() reads in both
curvature (.curv) and convexity (.sulc) files.

1. Functions for reading surfaces and converting between FreeSurfer formats
2. Functions for converting to VTK format
3. Functions specific to Mindboggle that call the read_surface() function


Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
# Functions for reading surfaces and converting between FreeSurfer formats
#=============================================================================

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
    >>> import readFreeSurfer as rfs
    >>> Vrtx, Face = rfs.read_surface('lh.pial')
    >>> len(Vrtx)
    130412
    >>> len(Face)
    260820
    >>> Vrtx[10]
    [-7.902474880218506, -95.6839370727539, -21.856534957885742]
    >>> Face[10]
    [2, 39, 3]

    """

    import os
    import struct

    f = open(filename, "rb")
    f.seek(3)  # skip the first 3 Bytes "Magic" number

    s = f.read(50)   # the second field is a string of variable length
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
    >>> Curv = read_curvature('lh.curv')
    >>> len(Curv)
    130412
    >>> Curv[10]
    -0.37290969491004944

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

def labels_to_volume(subject, annot_name):
    """
    Propagate surface labels through a gray matter volume
    using FreeSurfer's mri_aparc2aseg.

    Note
    ----
    From the mri_aparc2aseg documentation:
    The volumes of the cortical labels will be different than
    reported by mris_anatomical_stats because partial volume information
    is lost when mapping the surface to the volume. The values reported by
    mris_anatomical_stats will be more accurate than the volumes from the
    aparc+aseg volume.

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

def thickness_to_ascii(hemi, subject, subjects_path):
    """
    Convert a FreeSurfer thickness (per-vertex) file
    to an ascii file.

.. note::
    Untested function

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
# Functions for converting to VTK format
#=============================================================================

def surf_to_vtk(surface_file):
    """
    Convert FreeSurfer surface file to VTK format.
    """
    import os
    from utils.io_vtk import write_vtk_header, write_vtk_points, \
        write_vtk_faces
    from utils.io_free import read_surface

    Vertex, Face = read_surface(surface_file)

    vtk_file = os.path.join(os.getcwd(),
                            os.path.basename(surface_file + '.vtk'))
    Fp = open(vtk_file, 'w')
    write_vtk_header(Fp, Title='vtk output from ' + surface_file)
    write_vtk_points(Fp, Vertex)
    write_vtk_faces(Fp, Face)
    Fp.close()

    return vtk_file

def curv_to_vtk(file_string, surface_file, hemi, subject, subjects_path):
    """
    Convert FreeSurfer curvature, thickness, or convexity file to VTK format.

    Parameters
    ----------
    file_string : string
        string for FreeSurfer file: 'curv', 'thickness', 'sulc.pial'
    surface_file : string  (name of VTK surface file)
    hemi : string indicating left or right hemisphere
    subject : string
        name of subject directory
    subjects_path: string
        path to subject directory

    Returns
    -------
    vtk_file : string
        name of output VTK file, where each vertex is assigned
        the corresponding shape value.

    """
    import os
    from utils.io_free import read_curvature
    from utils.io_vtk import load_scalar, write_scalars

    filename = os.path.join(subjects_path, subject, 'surf',
                            hemi + '.' + file_string)
    vtk_file = os.path.join(os.getcwd(), file_string + '.vtk')

    curvature_values = read_curvature(filename)

    # Load VTK surface
    Points, Faces, Scalars = load_scalar(surface_file, return_arrays=0)
    Vertices =  range(1, len(Points) + 1)

    LUTs = [curvature_values]
    LUT_names = [file_string]
    write_scalars(vtk_file, Points, Vertices, Faces, LUTs, LUT_names)

    return vtk_file

def annot_to_vtk(surface_file, hemi, subject, subjects_path, annot_name):
    """
    Load a FreeSurfer .annot file and save as a VTK format file.

    Parameters
    ----------
    surface_file : string  (name of VTK surface file)
    annot_file : strings  (name of FreeSurfer .annot file)

    Returns
    -------
    labels : list of integers (one label per vertex)
    vtk_file : output VTK file

    """

    import os
    import nibabel as nb
    from utils.io_vtk import load_scalar, write_scalars

    annot_file = os.path.join(subjects_path, subject, 'label',
                              hemi + '.' + annot_name)

    labels, colortable, names = nb.freesurfer.read_annot(annot_file)

    # Load FreeSurfer surface
    #from utils.io_file import read_surface
    #Points, Faces = read_surface(surface_file)

    # Load VTK surface
    Points, Faces, Scalars = load_scalar(surface_file, return_arrays=0)
    Vertices =  range(1, len(Points) + 1)

    output_stem = os.path.join(os.getcwd(),
                  os.path.basename(surface_file.strip('.vtk')))
    vtk_file = output_stem + '.' + annot_name.strip('.annot') + '.vtk'

    LUTs = [labels.tolist()]
    LUT_names = ['Labels']
    write_scalars(vtk_file, Points, Vertices, Faces, LUTs, LUT_names)

    return labels, vtk_file

def vtk_to_label_files(hemi, surface_file, label_numbers, label_names,
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
    from utils import io_file
    import vtk

    # Check type to make sure the filename is a string
    # (if a list, return the first element)
    surface_file = io_file.string_vs_list_check(surface_file)

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


#=============================================================================
# Functions specific to Mindboggle that call the read_surface() function
#=============================================================================
def vertex_list_to_vtk(vtk_file, surface_file, index_pair_file,
                       LUT=[], LUTname=[]):
    """
    Load a vertex list file and a surface file to map vertices onto the surface
    and save the result into VTK format.

    Parameters
    ----------
    LUT : list of LUTs
        LUT[i] is the i-th LUT, and consists of a list of float numbers,
        representing for example distance transformation values, curvature, depth, etc.
        More than one LUT may be written to the final VTK file

    LUTname : list of strings
        LUTname[i] is the name of the lookup table to be inserted into VTK file

    """

    import os
    from utils import io_vtk, io_free

    vtk_file = os.path.join(os.getcwd(), vtk_file)

    Fp = open(vtk_file,'w')
    Vertex, Face = io_free.read_surface(surface_file)
    index_pair_list = io_vtk.load_fundi_list(index_pair_file)
    io_vtk.write_vertices_to_fundi(Fp, Vertex, index_pair_list)
    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            if i == 0:
                io_vtk.write_vtk_LUT(Fp, LUT[i], LUTname[i])
            else:
                io_vtk.write_vtk_LUT(Fp, LUT[i], LUTname[i],
                                        at_LUT_begin=False)
    Fp.close()

def face_list_to_vtk(vtk_file, surface_file, index_pair_file,
                     LUT=[], LUTname=[]):
    """
    Load a face list file and a surface file to map faces onto the surface
    and save the result into VTK format.

    This function is called by libbasin.getBasin()

    """

    import os
    from utils import io_vtk, io_free

    vtk_file = os.path.join(os.getcwd(), vtk_file)

    Fp = open(vtk_file,'w')
    Vertex, Face = io_free.read_surface(surface_file)
    index_pair_list = io_vtk.load_fundi_list(index_pair_file)
    io_vtk.write_feature_to_face(Fp, Vertex, Face, index_pair_list)

    if LUT!=[] :
        for i in xrange(0, len(LUT)):
            if i == 0:
                io_vtk.write_vtk_LUT(Fp, LUT[i], LUTname[i])
            else:
                io_vtk.write_vtk_LUT(Fp, LUT[i], LUTname[i],
                                        at_LUT_begin=False)
    Fp.close()

def line_segments_to_vtk(vtk_file, surface_file, index_pair_file, LUT=[], LUTname=[]):
    """
    Load a fundus curve segment list file and a surface file
    to map curve segments onto the surface and save the result into VTK format.

    Parameters
    ----------
    LUT : list of LUTs
        LUT[i] is the i-th LUT, e.g., distance transformation values,
                                      curvature, convexity values.
        LUT[i] is a list of float numbers
        So more than one LUTs may be written to the final VTK file

    LUTname : list of strings
        LUTname[i] is the name of the lookup table to be inserted into VTK file

    """

    import os
    from utils import io_vtk, io_free

    vtk_file = os.path.join(os.getcwd(), vtk_file)

    Fp = open(vtk_file,'w')
    Vertex, Face = io_free.read_surface(surface_file)
    index_pair_list = io_vtk.load_segmented_fundi(index_pair_file)
    io_vtk.write_line_segments_to_fundi(Fp, Vertex, index_pair_list)

    if LUT!=[]:
        for i in xrange(0, len(LUT)):
            if i == 0:
                io_vtk.write_vtk_LUT(Fp, LUT[i], LUTname[i])
            else:
                io_vtk.write_vtk_LUT(Fp, LUT[i], LUTname[i],
                                        at_LUT_begin=False)

    Fp.close()

