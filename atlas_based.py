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

def register_template(subject_id, subjects_path, 
                      template_name, templates_path, reg_name):
    """
    Register surface to template with FreeSurfer's mris_register

    Example: bert
             /Applications/freesurfer/subjects
             ./templates_freesurfer
             KKI_2.tif
             sphere_to_template.reg
    """
    from os import path
    from nipype.interfaces.base import CommandLine

    for hemi in ['lh','rh']:
        input_file = path.join(subjects_path, subject_id, 'surf', hemi + '.sphere')
        output_file = path.join(subjects_path, subject_id, 'surf', hemi + '.' + reg_name)
        template_file = path.join(templates_path, hemi + '.' + template_name)
        cli = CommandLine(command='mris_register')
        cli.inputs.args = ' '.join(['-curv', input_file, template_file, output_file])
        cli.cmdline
    return reg_name

"""
def register_atlases(subject_id, subjects_path, atlas_list, 
                     annot_name, reg_name):
    ""
    Transform the labels from multiple atlases via a template
    (using FreeSurfer's mri_surf2surf)
    ""
    from os import path
    from nipype.interfaces.base import CommandLine

    # For each atlas
    for atlas_name in atlas_list:
        # For each hemisphere
        for hemi in ['lh','rh']:        
            annot_file = path.join(subjects_path, atlas_name, 'label',
                                   hemi + '.' + annot_name) 
            output_annot = hemi + '.' + atlas_name + '_to_' + \
                           subject_id + '_' + annot_name
            args = ['--hemi', hemi,
                    '--srcsubject', atlas_name,
                    '--trgsubject', subject_id,
                    '--sval-annot', annot_file,
                    '--tval', output_annot,
                    '--srcsurfreg', reg_name,
                    '--trgsurfreg', reg_name]
            cli = CommandLine(command='mri_surf2surf')
            cli.inputs.args = ' '.join(args)
            cli.cmdline
    return annot_name

"""
def register_atlases(subject_id, subjects_path, atlas_list, 
                     annot_name, reg_name):
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
    from os import path
    from nipype.interfaces.freesurfer import SurfaceTransform

    sxfm = SurfaceTransform()
    sxfm.inputs.target_subject = subject_id

    # For each atlas
    for atlas_name in atlas_list:
        sxfm.inputs.source_subject = atlas_name

        # For each hemisphere
        for hemi in ['lh','rh']:        
            sxfm.inputs.hemi = hemi

            # Source file
            source_file = path.join(subjects_path, atlas_name, 'label',
                                    hemi + '.' + annot_name) 
            sxfm.inputs.source_file = source_file
            sxfm.inputs.svalannot = source_file

            # Output annotation file
            sxfm.inputs.tval = hemi + '.' + atlas_name + '_to_' + \
                               subject_id + '_' + annot_name

            # Arguments: strings within registered files
            args = ['--srcsurfreg', reg_name,
                    '--trgsurfreg', reg_name]
            sxfm.inputs.args = ' '.join(args)

            sxfm.run()

    return annot_name

##############################################################################
#   Multi-atlas labeling
##############################################################################

def combine_labels(label_index):
    """
    Replace some labels by combined labels. 

    Parameters
    ============
    label_index: integer label of a vertex

    Returns
    =========
    No variable. Direct return. 

    Notes
    ======= 
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
    ----------
    filepath : str  (path to annotation file)
    orig_ids : bool  (?return the vertex ids as stored in the annotation file 
                       or colortable ids)
    
    Returns
    -------
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

def load_labels(annot_path, annot_name):
    """
    Load multiple annotation files for each of the hemispheres of a subject
    
    Parameters
    ===========
    annot_name: string  (identifies annot files by their file name)
    left_files: list of strings  (annotation files for the left hemisphere)
    right_files: list of strings  (annotation files for the right hemisphere)
    AnnotPath: string  (path where all annotation files are saved)
    
    Returns 
    ========
    left_labels: list of lists of integers  (labels for left vertices)      
    right_labels: list of lists of integers  (labels for right vertices) 
    
    Notes
    ======
    The labels from read_annot range from 1 to 35 

    """
    
    print("Loading multiple annotations...")
    
    from os import listdir, path

    all_files  = listdir(annot_path)
    left_files, right_files = [], []
    for file1 in all_files:
        if file1.find(annot_name) > -1 and file1.find('_to_') > -1:
            if file1[0] == 'l':
                left_files.append(file1)
            elif file1[0] == 'r':        
                right_files.append(file1)
            else:
                print("Unable to match any file names.")

    left_labels, right_labels = [], []
    for file1 in left_files:
        Labels, ColorTable, Names = read_annot(path.join(annot_path, file1))
        left_labels.append(map(combine_labels,Labels))

    for file1 in right_files:
        Labels, ColorTable, Names = read_annot(path.join(annot_path, file1))
        right_labels.append(map(combine_labels,Labels))
    
    print("Multiple annotations loaded.")
    
    return left_labels, right_labels
    
def vote_labels(left_labels, right_labels):
    """
    For each vertex, vote on the majority label.

    Parameters 
    ========
    left_labels: list of lists of integers  (labels for left vertices)          
    right_labels: list of lists of integers  (labels for right vertices)
    n_vertices_left: integer  (number of left vertices)
    n_vertices_right: integer  (number of right vertices)
        
    Returns
    ========
    left_max: list of integers  (majority labels for left vertices)  
    right_max: list of integers  (majority labels for right vertices)
    left_counts: list of integers  (number of different labels for left vertices)
    right_counts: list of integers  (number of different labels for right vertices)
    left_votes: list of integers  (number of votes for the left majority labels) 
    right_votes: list of integers  (number of votes for the right majority labels)
  
    """
    from collections import Counter
    
    print("Begin voting...")
    
    n_atlases = len(left_labels)  # number of atlases used to label subject
    n_vertices_left, n_vertices_right = len(left_labels[0]), len(right_labels[0])
    left_max = [-1 for i in xrange(n_vertices_left)]  
    right_max = [-1 for i in xrange(n_vertices_right)]
    left_counts = [1 for i in xrange(n_vertices_left)]   
    right_counts = [1 for i in xrange(n_vertices_right)]
    left_votes = [n_atlases for i in xrange(n_vertices_left)]   
    right_votes = [n_atlases for i in xrange(n_vertices_right)]
    
    """
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
    for vertex in xrange(n_vertices_left):
        votes = Counter([left_labels[i][vertex] for i in xrange(n_atlases)])

        left_max[vertex] = votes.most_common(1)[0][0]
        left_votes[vertex] = votes.most_common(1)[0][1]
        left_counts[vertex] = len(votes)
            
    for vertex in xrange(n_vertices_right):
        votes = Counter([right_labels[i][vertex] for i in xrange(n_atlases)])
        
        right_max[vertex] = votes.most_common(1)[0][0]
        right_votes[vertex] = votes.most_common(1)[0][1]
        right_counts[vertex] = len(votes)
        
    print("Voting done.")
    return left_max, right_max,\
           left_votes, right_votes,\
           left_counts, right_counts 
 
def multilabel(subject_id, subjects_path, annot_name, use_inflated_surfaces=1):
    """
    Load VTK surfaces and write majority vote labels as VTK files, 
    according to multiple labelings (annot_name).
    """
    from os import path
    import pyvtk
    from atlas_based import load_labels, vote_labels

    subject_surf_path = path.join(subjects_path, subject_id, 'surf')
    annot_path = path.join(subjects_path, subject_id, 'label')

    # Load multiple label sets
    left_labels, right_labels = load_labels(annot_path, annot_name)
 
    # Vote on labels for each vertex
    left_max, right_max, left_votes, right_votes,\
    left_counts, right_counts = vote_labels(left_labels, right_labels)
    
    if use_inflated_surfaces:
        input_files = ['lh.pial.vtk', 'rh.pial.vtk', 
                       'lh.inflated.vtk', 'rh.inflated.vtk']
        output_files = [['lh.pial.labels.max.vtk',
                         'lh.pial.labelcounts.vtk', 
                         'lh.pial.labelvotes.vtk'], 
                        ['rh.pial.labels.max.vtk', 
                         'rh.pial.labelcounts.vtk', 
                         'rh.pial.labelvotes.vtk'], 
                        ['lh.inflated.labels.max.vtk', 
                         'lh.inflated.labelcounts.vtk', 
                         'lh.inflated.labelvotes.vtk'], 
                        ['rh.inflated.labels.max.vtk', 
                         'rh.inflated.labelcounts.vtk', 
                         'rh.inflated.labelvotes.vtk']]
        dv = [[left_max, left_counts, left_votes],
              [right_max, right_counts, right_votes],
              [left_max, left_counts, left_votes],
              [right_max, right_counts, right_votes]]
    else:
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
              tofile(path.join(annot_path, output_files[i][0]), 'ascii')
        pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
              pyvtk.PointData(pyvtk.Scalars(dv[i][1],\
                    name='Counts (number of different labels)'))).\
              tofile(path.join(annot_path, output_files[i][1]), 'ascii')
        pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
              pyvtk.PointData(pyvtk.Scalars(dv[i][2],\
                    name='Votes (number of votes for majority labels)'))).\
              tofile(path.join(annot_path, output_files[i][2]), 'ascii')

    return annot_name  
    
