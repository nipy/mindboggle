#!/usr/bin/env python
"""
Functions for reading surfaces and converting between FreeSurfer formats.

This Python library reads and writes different file types.
The function read_surface() reads in surface files,
while the function read_curvature() reads in both
curvature (.curv) and convexity (.sulc) files.


Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import os
#import struct
#from nipype.interfaces.base import CommandLine


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
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                        'surf', 'lh.pial')
    >>> Vrtx, Face = rfs.read_surface(surface_file)
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
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> surface_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                        'surf', 'lh.curv')
    >>> Curv = read_curvature(surface_file)
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
