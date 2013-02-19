#!/usr/bin/env python
"""
Atlas-based functions for surface registration-based labeling:

Register to template
    Register surface to template with FreeSurfer's mris_register.
    Transform the labels from multiple atlases via a template
    (using FreeSurfer's mri_surf2surf).

Transform atlas labels
    For each brain hemisphere (left and right) in a given subject,
    read in FreeSurfer *.annot files (multiple labelings) and output one VTK file
    of majority vote labels, representing a "maximum probability" labeling.
    The main function is majority_vote() and calls vote_labels().


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

##############################################################################
#   Template-based, multi-atlas registration
##############################################################################

def register_template(hemi, sphere_file, transform,
                      templates_path, template):
    """
    Register surface to template with FreeSurfer's mris_register.
    """
    from os import path
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    template_file = path.join(templates_path, hemi + '.' + template)
    output_file = hemi + '.' + transform
    cli = CommandLine(command='mris_register')
    cli.inputs.args = ' '.join(['-curv', sphere_file,
                                template_file, output_file])
    logger.info(cli.cmdline)
    cli.run()

    return transform

def transform_atlas_labels(hemi, subject, transform,
                           subjects_path, atlas, atlas_string):
    """
    Transform the labels from a surface atlas via a template
    using FreeSurfer's mri_surf2surf (wrapped in NiPype).

    nipype.workflows.smri.freesurfer.utils.fs.SurfaceTransform
    wraps command ``mri_surf2surf``:

        "Transform a surface file from one subject to another via a spherical registration.
        Both the source and target subject must reside in your Subjects Directory,
        and they must have been processed with recon-all, unless you are transforming
        to one of the icosahedron meshes."

    Parameters
    ----------
    hemi : ``string``: hemisphere
    subject : ``string``: subject, corresponding to FreeSurfer subject directory
    transform : ``string``: FreeSurfer spherical surface registration transform
    subjects_path : ``string``: FreeSurfer subjects directory
    atlas : ``string``: name of atlas
    atlas_string : ``string``: name of atlas labeling protocol

    """
    from os import path, getcwd
    from nipype.interfaces.freesurfer import SurfaceTransform

    sxfm = SurfaceTransform()
    sxfm.inputs.hemi = hemi
    sxfm.inputs.target_subject = subject
    sxfm.inputs.source_subject = atlas

    # Source file
    sxfm.inputs.source_annot_file = path.join(subjects_path,
                                    atlas, 'label',
                                    hemi + '.' + atlas_string + '.annot')
    # Output annotation file
    output_file = path.join(getcwd(), hemi + '.' + atlas + '.' + atlas_string + \
                                      '_to_' + subject + '.annot')
    sxfm.inputs.out_file = output_file

    # Arguments: strings within registered files
    args = ['--srcsurfreg', transform,
            '--trgsurfreg', transform]
    sxfm.inputs.args = ' '.join(args)

    sxfm.run()

    return output_file

##############################################################################
#   Multi-atlas labeling
##############################################################################

def vote_labels(label_lists):
    """
    For each vertex, vote on the majority label.

    Parameters
    ----------
    label_lists : list of lists of integers  (vertex labels assigned by each atlas)
    n_atlases : integer  (number of atlases / lists of labels)
    npoints : integer  (number of vertices / elements in each list)

    Returns
    -------
    labels_max : list of integers  (majority labels for vertices)
    label_counts : list of integers  (number of different labels for vertices)
    label_votes : list of integers  (number of votes for the majority labels)

    Examples
    --------
    >>> from collections import Counter
    >>> X = [1,1,2,3,4,2,1,2,1,2,1,2]
    >>> Votes = Counter(X)
    >>> Votes
    Counter({1: 5, 2: 5, 3: 1, 4: 1})
    >>> Votes.most_common(1)
    [(1, 5)]
    >>> Votes.most_common(2)
    [(1, 5), (2, 5)]
    >>> len(Votes)
    4

    """

    from collections import Counter

    print("Begin voting...")
    n_atlases = len(label_lists)  # number of atlases used to label subject
    npoints = len(label_lists[0])
    labels_max = [-1 for i in xrange(npoints)]
    label_counts = [1 for i in xrange(npoints)]
    label_votes = [n_atlases for i in xrange(npoints)]

    consensus_vertices = []
    for vertex in xrange(npoints):
        votes = Counter([label_lists[i][vertex] for i in xrange(n_atlases)])

        labels_max[vertex] = votes.most_common(1)[0][0]
        label_votes[vertex] = votes.most_common(1)[0][1]
        label_counts[vertex] = len(votes)
        if len(votes) == n_atlases:
            consensus_vertices.append(vertex)

    print("Voting done.")

    return labels_max, label_votes, label_counts, consensus_vertices

def majority_vote_label(surface_file, annot_files):
    """
    Load a VTK surface and corresponding FreeSurfer annot files.
    Write majority vote labels, and label counts and votes as VTK files.

    Parameters
    ----------
    surface_file : string  (name of VTK surface file)
    annot_files : list of strings  (names of FreeSurfer annot files)

    Returns
    -------
    labels_max : list of integers  (majority labels for vertices)
    label_counts : list of integers  (number of different labels for vertices)
    label_votes : list of integers  (number of votes for the majority labels)
    consensus_vertices : list of integers  (indicating which are consensus labels)
    maxlabel_file : VTK file containing majority vote labels
    labelcounts_file : VTK file containing number of different label counts
    labelvotes_file : VTK file containing number of votes per majority label

    """
    from os import path, getcwd
    import nibabel as nb
    import pyvtk
    from mindboggle.label.multiatlas_labeling import vote_labels
    from mindboggle.utils.io_file import string_vs_list_check

    # Load multiple label sets
    print("Load annotation files...")
    label_lists = []
    for annot_file in annot_files:
        labels, colortable, names = nb.freesurfer.read_annot(annot_file)
        label_lists.append(labels)
    print("Annotations loaded.")

    # Vote on labels for each vertex
    labels_max, label_votes, label_counts, \
    consensus_vertices = vote_labels(label_lists)

    # Check type to make sure the filename is a string
    # (if a list, return the first element)
    surface_file = string_vs_list_check(surface_file)

    # Save files
    VTKReader = pyvtk.VtkData(surface_file)
    Vertices =  VTKReader.structure.points
    Faces =     VTKReader.structure.polygons

    output_stem = path.join(getcwd(), path.basename(surface_file.strip('.vtk')))
    maxlabel_file = output_stem + '.labels.max.vtk'
    labelcounts_file = output_stem + '.labelcounts.vtk'
    labelvotes_file = output_stem + '.labelvotes.vtk'

    pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(labels_max,\
                                  name='Max (majority labels)'))).\
          tofile(maxlabel_file, 'ascii')

    pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(label_counts,\
                                  name='Counts (number of different labels)'))).\
          tofile(labelcounts_file, 'ascii')

    pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(label_votes,\
                                  name='Votes (number of votes for majority labels)'))).\
          tofile(labelvotes_file, 'ascii')

    return labels_max, label_counts, label_votes, consensus_vertices, \
           maxlabel_file, labelcounts_file, labelvotes_file
