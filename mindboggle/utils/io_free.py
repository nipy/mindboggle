#!/usr/bin/python

"""
FreeSurfer functions

This Python library reads and writes different file types.
In particular, it has functions to read some FreeSurfer files,
including surface, curvature, and convexity files.
The function read_surface reads in surface files,
while the function read_curvature reads in both
curvature (.curv) and convexity (.sulc) files.

Authors:
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com
Forrest Sheng Bao  .  http://fsbao.net

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def read_surface(filename):
    """
    Read in a FreeSurfer triangle surface mesh in binary format.

    Parameters
    ===========
    filename : string
        A binary FreeSurfer Triangle Surface file

    Outputs
    =======
    Vertex : list of 3-tuples of floats
        Each element is a 3-tuple (list) of floats, which are the X, Y and Z coordinates of a vertex, respectively.
        A 3-tuple's index in the list *Vertex* is the ID of a vertex.

    Face : list of 3-tuples of integers
        Each element is a 3-tuple (list) of integers, which are the IDs of 3 vertexes that form one face

    Example
    ========

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
        Vertex.append([R,A,S]) # R, A, S are the coordinates of vertexes

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
    Read in a FreeSurfer curvature (per-vertex) file.

    Parameters
    ==========

    filename : string
        A binary FreeSurfer curvature (pre-vertex) file

    Outputs
    ========

    Curvature : list of floats
        Each element is the curvature value of a FreeSurfer mesh vertex.
        Elements are ordered by orders of vertexes in FreeSurfer surface file.

    Example
    ========

    >>> import readFreeSurfer as rfs
    >>> Curv = rfs.read_curvature('lh.curv')
    >>> len(Curv)
    130412
    >>> Curv[10]
    -0.37290969491004944

    """

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

    if VertexCount != 0:  # number of vertexes is even (NOT ODD!!!)
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

    Inputs:
    ------
    hemi:  hemisphere [string]
    subjects_path:  path to file
    subject:  subject name
    label_files:  .label file names [list of strings]
    colortable:  file of label numbers & names (same order as label_files)
    annot_name:  name of the output .annot file (without prepending hemi)

    Output:
    ------
    annot_name, annot_file:  name of .annot file (without & with prepend)

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
