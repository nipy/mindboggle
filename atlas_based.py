#!/usr/bin/python

"""
Atlas-based functions for surface registration-based labeling:

1. Surface map conversion (from FreeSurfer to VTK format)

2. Multi-atlas registration
Register surface to template with FreeSurfer's mris_register.
Transform the labels from multiple atlases via a template 
(using FreeSurfer's mri_surf2surf).

3. Multi-atlas labeling
For each brain hemisphere (left and right) in a given subject, 
read in FreeSurfer *.annot files (multiple labelings) and output one VTK file 
of majority vote labels, representing a "maximum probability" labeling.
The main function is labeling() and calls: 
load_labels, vote_labels, combine_labels, and read_annot (0-index labels).


Authors:  
Forrest Bao  .  forrest.bao@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

##############################################################################
#   Surface map conversion
##############################################################################

def convert_to_vtk(fs_surface_files):
    """
    Measure

    # import nipype.interfaces.freesurfer as fs
    # mris = fs.MRIsConvert()
    # mris.inputs(in_file = '', out_file = '', subjects_dir = '')
    """
    #import subprocess as sp
    import os
    surface_files = []
    for fs_surface_file in fs_surface_files:
        surface_file = fs_surface_file + '.vtk'
        surface_files.append(surface_file)
        args = ['mris_convert', fs_surface_file, surface_file]
        print(' '.join(args)); os.system(' '.join(args))
    #proc = sp.Popen(' '.join(args))
    #o, e = proc.communicate()
    #if proc.returncode > 0 :
    #    raise Exception('\n'.join([args + ' failed', o, e]))
    return surface_files

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
    import os

    for hemi in ['lh','rh']:
        input_file = os.path.join(subjects_path, subject_id, 'surf', hemi + '.sphere')
        output_file = os.path.join(subjects_path, subject_id, 'surf', hemi + '.' + reg_name)
        template_file = os.path.join(templates_path, hemi + '.' + template_name)
        args = ['mris_register -curv', input_file, template_file, output_file]
        print(' '.join(args)); os.system(' '.join(args))
    return reg_name

"""
def register_atlases(subject_id, subjects_path, atlas_list_file, 
                     annot_name, reg_name):
    ""
    Transform the labels from multiple atlases via a template
    (using FreeSurfer's mri_surf2surf)
    ""
    import os

    # Get list of atlas subjects from a file
    f = open(atlas_list_file)
    atlas_list = f.readlines()
    for atlas_line in atlas_list:
        # For each atlas
        atlas_name = atlas_line.strip("\n")
        # For each hemisphere
        for hemi in ['lh','rh']:        
            annot_file = os.path.join(subjects_path, atlas_name, 'label',
                                      hemi + '.' + annot_name) 
            output_annot = hemi + '.' + atlas_name + '_to_' + \
                           subject_id + '_' + annot_name
            args = ['mri_surf2surf',
                    '--hemi', hemi,
                    '--srcsubject', atlas_name,
                    '--trgsubject', subject_id,
                    '--sval-annot', annot_file,
                    '--tval', output_annot,
                    '--srcsurfreg', reg_name,
                    '--trgsurfreg', reg_name]
            print(' '.join(args)); os.system(' '.join(args))
    return annot_name
"""

def register_atlases(subject_id, subjects_path, atlas_list_file, 
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
    import os
    from nipype.interfaces.freesurfer import SurfaceTransform

    sxfm = SurfaceTransform()
    sxfm.inputs.target_subject = subject_id

    # Get list of atlas subjects from a file
    f = open(atlas_list_file)
    atlas_list = f.readlines()
    for atlas_line in atlas_list:
        # For each atlas
        atlas_name = atlas_line.strip("\n")
        sxfm.inputs.source_subject = atlas_name
        # For each hemisphere
        for hemi in ['lh','rh']:        
            sxfm.inputs.hemi = hemi
            source_file = os.path.join(subjects_path, atlas_name, 'label',
                                       hemi + '.' + annot_name) 
            sxfm.inputs.source_file = source_file
            # Specify annotation file
            output_annot = os.path.join(subjects_path, subject_id, 'label',
                                        hemi + '.' + atlas_name + '_to_' + \
                                        subject_id + '_' + annot_name)
            sxfm.inputs.out_file = output_annot
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
        2= 2,10,23,26
        3= 3, 27
        18=18,19,20

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
    
    print "Loading multiple annotations..."
    
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
                print "Unable to match any file names."

    left_labels, right_labels = [], []
    for file1 in left_files:
        Labels, ColorTable, Names = read_annot(path.join(annot_path, file1))
        left_labels.append(map(combine_labels,Labels))

    for file1 in right_files:
        Labels, ColorTable, Names = read_annot(path.join(annot_path, file1))
        right_labels.append(map(combine_labels,Labels))
    
    print "Multiple annotations loaded."
    
    return left_labels, right_labels
    
def vote_labels(left_labels, right_labels):
    """
    For each vertex, vote on the majority label.

    Parameters 
    ========
    left_labels: list of lists of integers  (labels for left vertices)          
    right_labels: list of lists of integers  (labels for right vertices)
    left_n_labels: integer  (number of left vertices)
    right_n_labels: integer  (number of right vertices)
        
    Returns
    ========
    left_assign: list of integers  (winning labels for left vertices)  
    right_assign: list of integers  (winning labels for right vertices)
    left_consensus: list of integers  (number of left consensus labels) 
    right_consensus: list of integers  (number of right consensus labels)
  
    """
    from collections import Counter
    
    compute_diff_consensus = 0
    
    print "Begin voting..."
    
    n_labelings = len(left_labels)  # number of atlases used to label subject
    left_n_labels, right_n_labels = len(left_labels[0]), len(right_labels[0])
    left_assign = [-1 for i in xrange(left_n_labels)]  
                  # if no consistent vote, the label is -1 (vertex is unlabeled)
    right_assign = [-1 for i in xrange(right_n_labels)]

    if compute_diff_consensus:
        left_diff = [1 for i in xrange(left_n_labels)]   
        right_diff = [1 for i in xrange(right_n_labels)]
        left_consensus = [n_labelings for i in xrange(left_n_labels)]   
        right_consensus = [n_labelings for i in xrange(right_n_labels)]
    
    for vertex in xrange(left_n_labels):
        votes = Counter([left_labels[i][vertex] for i in xrange(n_labelings)])
        
        left_assign[vertex] = votes.most_common(1)[0][0]
        if compute_diff_consensus:
            left_consensus[vertex] = votes.most_common(1)[0][1]
            left_diff[vertex] = len(votes)
            
    for vertex in xrange(right_n_labels):
        votes = Counter([right_labels[i][vertex] for i in xrange(n_labelings)])
        
        right_assign[vertex] = votes.most_common(1)[0][0]
        if compute_diff_consensus:
            right_consensus[vertex] = votes.most_common(1)[0][1]
            right_diff[vertex] = len(votes)
        
    print "Voting done."
    if compute_diff_consensus:
        return left_assign, right_assign,\
               left_consensus, right_consensus, left_diff, right_diff 
    else:
        return left_assign, right_assign
 
def labeling(subject_id, subjects_path, annot_name):
    """
    Load VTK surfaces and write majority vote labels as VTK files, 
    according to multiple labelings (annot_name).
    """
    from os import path
    import pyvtk
    from atlas_based import load_labels, vote_labels

    compute_diff_consensus = 0
    save_inflated = 0

    subject_surf_path = path.join(subjects_path, subject_id, 'surf')
    annot_path = path.join(subjects_path, subject_id, 'label')
    
    # Load multiple label sets
    left_labels, right_labels = load_labels(annot_path, annot_name)
 
    # Vote on labels for each vertex
    if compute_diff_consensus:
        left_assign, right_assign, left_consensus, right_consensus,\
        left_diff, right_diff = vote_labels(left_labels, right_labels)
    else:
        left_assign, right_assign = vote_labels(left_labels, right_labels)
    
    if save_inflated:
        input_files = ['lh.pial.vtk', 'rh.pial.vtk', 
                       'lh.inflated.vtk', 'rh.inflated.vtk']
        output_files = ['lh.pial.labels.vtk', 'rh.pial.labels.vtk', 
                        'lh.inflated.labels.vtk', 'rh.inflated.labels.vtk']
        if compute_diff_consensus:
            dv = [[left_assign, left_diff, left_consensus],
                  [right_assign, right_diff, right_consensus],
                  [left_assign, left_diff, left_consensus],
                  [right_assign, right_diff, right_consensus]]
        else:
            dv = [[left_assign], [right_assign],
                  [left_assign], [right_assign]]
    else:
        input_files = ['lh.pial.vtk', 'rh.pial.vtk']
        output_files = ['lh.pial.labels.vtk', 'rh.pial.labels.vtk']
        if compute_diff_consensus:
            dv = [[left_assign, left_diff, left_consensus],
                  [right_assign, right_diff, right_consensus]]
        else:
            dv = [[left_assign], [right_assign]]

    for i, input_file in enumerate(input_files):
        VTKReader = pyvtk.VtkData(path.join(subject_surf_path, input_file))
        Vertices =  VTKReader.structure.points
        Faces =     VTKReader.structure.polygons
        
        if compute_diff_consensus:
            pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(dv[i][0], name='Assigned'),\
                                  pyvtk.Scalars(dv[i][1], name='Different'),\
                                  pyvtk.Scalars(dv[i][2], name='Common'))).\
                  tofile(path.join(annot_path, output_files[i]), 'ascii')
        else:
            pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(dv[i][0], name='Assigned'))).\
                  tofile(path.join(annot_path, output_files[i]), 'ascii')
    return annot_name  
    
