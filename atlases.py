#!/usr/bin/python

"""
Atlas-based functions for surface registration-based labeling:

1. Multi-atlas registration
Register surface to template with FreeSurfer's mris_register.
Transform the labels from multiple atlases via a template 
(using FreeSurfer's mri_surf2surf).

2. Multi-atlas labeling
For each brain hemisphere (left and right) in a given subject, 
read in FreeSurfer *.annot files (multiple labelings) and output one VTK file 
of majority vote labels, representing a "maximum probability" labeling.
The main function is multilabel() and calls: 
load_labels, vote_labels, combine_labels, and read_annot (0-index labels).


Authors:  
Forrest Bao  .  forrest.bao@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

##############################################################################
#   Multi-atlas registration
##############################################################################

def register_template(hemi, sph_surface_file, 
                      template_reg_name, template_name, templates_path):
    """
    Register surface to template with FreeSurfer's mris_register

    Example: lh
             bert
             /Applications/freesurfer/subjects
             sphere_to_template.reg
             KKI_2.tif
             ./templates_freesurfer
             
    """
    from os import path, getcwd
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    template_file = path.join(templates_path, hemi + '.' + template_name)
    output_file = hemi + '.' + template_reg_name
    cli = CommandLine(command='mris_register')
    cli.inputs.args = ' '.join(['-curv', sph_surface_file, 
                                template_file, output_file])
    logger.info(cli.cmdline)
    cli.run()
    
    return template_reg_name

"""
def register_atlas(hemi, subject_id, template_reg_name,
                   atlas_name, atlases_path, atlas_annot_name):
    ""
    Transform the labels from multiple atlases via a template
    (using FreeSurfer's mri_surf2surf)
    ""
    from os import system, path, getcwd

    source_annot_file = path.join(atlases_path, atlas_name, 'label',
                                  hemi + '.' + atlas_annot_name) 
    output_file = path.join(getcwd(), hemi + '.' + atlas_name + '_to_' + \
                            subject_id + '_' + atlas_annot_name)
    args = ['mri_surf2surf',
            '--hemi', hemi,
            '--srcsubject', atlas_name,
            '--trgsubject', subject_id,
            '--sval-annot', source_annot_file,
            '--tval', output_file,
            '--srcsurfreg', template_reg_name,
            '--trgsurfreg', template_reg_name]
    print(' '.join(args)); system(' '.join(args))
    return output_file
"""

def register_atlas(hemi, subject_id, template_reg_name,
                   atlas_name, atlases_path, atlas_annot_name):
    """
    Transform the labels from multiple atlases via a template
    using FreeSurfer's mri_surf2surf (wrapped in NiPype)

    nipype.workflows.smri.freesurfer.utils.fs.SurfaceTransform
    wraps command **mri_surf2surf**:
    "Transform a surface file from one subject to another via a spherical registration.
    Both the source and target subject must reside in your Subjects Directory,
    and they must have been processed with recon-all, unless you are transforming
    to one of the icosahedron meshes."
    """
    from os import path, getcwd
    from nipype.interfaces.freesurfer import SurfaceTransform

    sxfm = SurfaceTransform()
    sxfm.inputs.hemi = hemi
    sxfm.inputs.target_subject = subject_id
    sxfm.inputs.source_subject = atlas_name

    # Source file
    sxfm.inputs.source_annot_file = path.join(atlases_path, 
                                    atlas_name, 'label',
                                    hemi + '.' + atlas_annot_name) 
    # Output annotation file
    output_file = path.join(getcwd(), hemi + '.' + atlas_name + '_to_' + \
                            subject_id + '_' + atlas_annot_name)
    sxfm.inputs.out_file = output_file

    # Arguments: strings within registered files
    args = ['--srcsurfreg', template_reg_name,
            '--trgsurfreg', template_reg_name]
    sxfm.inputs.args = ' '.join(args)

    sxfm.run()

    return output_file

##############################################################################
#   Multi-atlas labeling
##############################################################################

def combine_labels(label_index):
    """
    Replace some labels by combined labels. 

    Parameters
    ==========
    label_index: integer label of a vertex

    Returns
    =======
    No variable. Direct return. 

    Notes
    ===== 
    label_index combinations:

    2  = 2, 10, 23, 26
    3  = 3, 27  
    18 = 18, 19, 20
    
    Temporal (33) and frontal (32) poles, and bankstss (1) regions eliminated, 
    corresponding cortex absorbed by adjacent regions.

    Caudal (2), isthmus (10), posterior (23), and rostral anterior (26) cingulate 
    combined to form single cingulate region (2)
 
    Caudal (3) and rostral (27) middle frontal regions combined to form 
    single middle frontal region (3)
    
    Opercular (18), orbital (19), and triangular (20) inferior frontal regions 
    combined to form a single inferior frontal region (18)
    
    """
    
    if label_index == 10 or label_index == 23 or label_index == 26:
        return 2
    elif label_index == 27:
        return 3
    elif label_index == 19 or label_index == 20:
        return 18
    else:
        return label_index

def read_annot(filepath, orig_ids=False):
    """
    Read in a Freesurfer annotation from a .annot file.
    From https://github.com/nipy/PySurfer/blob/master/surfer/io.py

    Parameters
    ==========
    filepath : str  (path to annotation file)
    orig_ids : bool  (?return the vertex ids as stored in the annotation file 
                       or colortable ids)
    
    Returns
    =======
    labels : n_vtx numpy array  (annotation id at each vertex)
    ctab : numpy array  (RGBA + label id colortable array)
    names : numpy array  (array of region names as stored in the annot file)   

    """
    import numpy as np
    
    with open(filepath, "rb") as fobj:
        dt = ">i4"
        vnum = np.fromfile(fobj, dt, 1)[0]
        data = np.fromfile(fobj, dt, vnum * 2).reshape(vnum, 2)
        labels = data[:, 1]
        ctab_exists = np.fromfile(fobj, dt, 1)[0]
        if not ctab_exists:
            raise Exception('Color table not found in annotation file.')
        n_entries = np.fromfile(fobj, dt, 1)[0]
        if n_entries > 0:
            length = np.fromfile(fobj, dt, 1)[0]
            orig_tab = np.fromfile(fobj, '>c', length)
            orig_tab = orig_tab[:-1]

            names = list()
            ctab = np.zeros((n_entries, 5), np.int)
            for i in xrange(n_entries):
                name_length = np.fromfile(fobj, dt, 1)[0]
                name = np.fromfile(fobj, "|S%d" % name_length, 1)[0]
                names.append(name)
                ctab[i, :4] = np.fromfile(fobj, dt, 4)
                ctab[i, 4] = (ctab[i, 0] + ctab[i, 1] * (2 ** 8) +
                              ctab[i, 2] * (2 ** 16) +
                              ctab[i, 3] * (2 ** 24))
        else:
            ctab_version = -n_entries
            if ctab_version != 2:
                raise Exception('Color table version not supported.')
            n_entries = np.fromfile(fobj, dt, 1)[0]
            ctab = np.zeros((n_entries, 5), np.int)
            length = np.fromfile(fobj, dt, 1)[0]
            _ = np.fromfile(fobj, "|S%d" % length, 1)[0] # Orig table path
            entries_to_read = np.fromfile(fobj, dt, 1)[0]
            names = list()
            for i in xrange(entries_to_read):
                _ = np.fromfile(fobj, dt, 1)[0] # Structure
                name_length = np.fromfile(fobj, dt, 1)[0]
                name = np.fromfile(fobj, "|S%d" % name_length, 1)[0]
                names.append(name)
                ctab[i, :4] = np.fromfile(fobj, dt, 4)
                ctab[i, 4] = (ctab[i, 0] + ctab[i, 1] * (2 ** 8) +
                                ctab[i, 2] * (2 ** 16))
        ctab[:, 3] = 255
    if not orig_ids:
        ord = np.argsort(ctab[:, -1])
        labels = ord[np.searchsorted(ctab[ord, -1], labels)]
    return labels, ctab, names 

def load_labels_from_annot(annot_file):
    """
    Load annotation file

    Parameters
    ==========
    annot_file: string
    atlas_annot_name: string  (identifies annot file)

    Returns 
    =======
    labels: list of lists of integers  (vertex labels)

    """

    from os import listdir, path

    print("Loading multiple annotations...")

    labels, colortable, names = read_annot(annot_file)

    print("Annotation loaded.")

    return labels

def vote_labels(labels):
    """
    For each vertex, vote on the majority label.

    Parameters 
    ==========
    left_labels: list of lists of integers  (labels for left vertices)          
    right_labels: list of lists of integers  (labels for right vertices)
    n_vertices_left: integer  (number of left vertices)
    n_vertices_right: integer  (number of right vertices)
        
    Returns
    =======
    left_max: list of integers  (majority labels for left vertices)  
    right_max: list of integers  (majority labels for right vertices)
    left_counts: list of integers  (number of different labels for left vertices)
    right_counts: list of integers  (number of different labels for right vertices)
    left_votes: list of integers  (number of votes for the left majority labels) 
    right_votes: list of integers  (number of votes for the right majority labels)

    Example of Counter:
    In [1]: from collections import Counter
    In [2]: X = [1,1,2,3,4,2,1,2,1,2,1,2]
    In [4]: Votes = Counter(X)
    In [7]: Votes
    Out[7]: Counter({1: 5, 2: 5, 3: 1, 4: 1})
    In [5]: Votes.most_common(1)
    Out[5]: [(1, 5)]
    In [6]: Votes.most_common(2)
    Out[6]: [(1, 5), (2, 5)]
    In [8]: len(Votes)
    Out[8]: 4

    """

    from collections import Counter
    
    print("Begin voting...")
    n_atlases = len(labels)  # number of atlases used to label subject
    n_vertices = len(labels[0])
    labels_max = [-1 for i in xrange(n_vertices)]  
    label_counts = [1 for i in xrange(n_vertices)]
    label_votes = [n_atlases for i in xrange(n_vertices)]

    for vertex in xrange(n_vertices):
        votes = Counter([labels[i][vertex] for i in xrange(n_atlases)])

        labels_max[vertex] = votes.most_common(1)[0][0]
        label_votes[vertex] = votes.most_common(1)[0][1]
        label_counts[vertex] = len(votes)

    print("Voting done.")

    return labels_max, label_votes, label_counts
 
def multilabel(hemi, subject_id, subjects_path, in_files, atlases_path, atlas_annot_name):
    """
    Load VTK surfaces and write majority vote labels as VTK files, 
    according to multiple labelings (atlas_annot_name).
    """
    from os import path, getcwd
    import pyvtk
    from atlases import load_labels, vote_labels

#    subject_surf_path = path.join(subjects_path, subject_id, 'surf')
#    atlas_annot_path = path.join(atlases_path, subject_id, 'label')

    # Load multiple label sets
#    labels = []
    labels = load_labels(atlas_annot_path, atlas_annot_name)
#    labels.append(map(combine_labels,Labels))

    # Vote on labels for each vertex
    left_max, right_max, left_votes, right_votes,\
    left_counts, right_counts = vote_labels(left_labels, right_labels)

        input_files = ['lh.pial.vtk', 'rh.pial.vtk']
        output_files = [['lh.pial.labels.max.vtk', 
                         'lh.pial.labelcounts.vtk', 
                         'lh.pial.labelvotes.vtk'], 
                        ['rh.pial.labels.max.vtk', 
                         'rh.pial.labelcounts.vtk', 
                         'rh.pial.labelvotes.vtk']]
        dv = [[left_max, left_counts, left_votes],
              [right_max, right_counts, right_votes]]

    for i, input_file in enumerate(input_files):
        VTKReader = pyvtk.VtkData(path.join(subject_surf_path, input_file))
        Vertices =  VTKReader.structure.points
        Faces =     VTKReader.structure.polygons
        pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
              pyvtk.PointData(pyvtk.Scalars(dv[i][0],\
                    name='Max (majority labels)'))).\
              tofile(path.join(getcwd(), output_files[i][0]), 'ascii')
        pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
              pyvtk.PointData(pyvtk.Scalars(dv[i][1],\
                    name='Counts (number of different labels)'))).\
              tofile(path.join(getcwd(), output_files[i][1]), 'ascii')
        pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
              pyvtk.PointData(pyvtk.Scalars(dv[i][2],\
                    name='Votes (number of votes for majority labels)'))).\
              tofile(path.join(getcwd(), output_files[i][2]), 'ascii')

    return output_files  
    
