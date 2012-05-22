#!/usr/bin/python

"""
Combine labels

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
#   Combine labels
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
    
def relabel(subject_id, subjects_path, annot_name):
    """
    Load VTK surfaces and combine labels.
    """
    from os import path
    import pyvtk
    from combine_labels import load_labels

    save_inflated = 0

    subject_surf_path = path.join(subjects_path, subject_id, 'surf')
    annot_path = path.join(subjects_path, subject_id, 'label')
    
    # Load multiple label sets
    left_labels, right_labels = load_labels(annot_path, annot_name)
 
    if save_inflated:
        input_files = ['lh.pial.vtk', 'rh.pial.vtk', 
                       'lh.inflated.vtk', 'rh.inflated.vtk']
        output_files = ['lh.pial.labels.vtk', 'rh.pial.labels.vtk', 
                        'lh.inflated.labels.vtk', 'rh.inflated.labels.vtk']
        dv = [[left_assign], [right_assign],
              [left_assign], [right_assign]]
    else:
        input_files = ['lh.pial.vtk', 'rh.pial.vtk']
        output_files = ['lh.pial.labels.vtk', 'rh.pial.labels.vtk']
        dv = [[left_assign], [right_assign]]

    for i, input_file in enumerate(input_files):
        VTKReader = pyvtk.VtkData(path.join(subject_surf_path, input_file))
        Vertices =  VTKReader.structure.points
        Faces =     VTKReader.structure.polygons
        
        pyvtk.VtkData(pyvtk.PolyData(points=Vertices, polygons=Faces),\
              pyvtk.PointData(pyvtk.Scalars(dv[i][0], name='Assigned'))).\
              tofile(path.join(annot_path, output_files[i]), 'ascii')
    
"""
def labelMap(label_index):
    # Given a label as in http://surfer.nmr.mgh.harvard.edu/fswiki/AnnotFiles, 
    # return its index
    
    Map = {1639705:1, 2647065:2, 10511485:3, 6500:4, 3294840:5,\
           6558940:6, 660700:7, 9231540:8, 14433500:9, 7874740:10,\
           9180300:11, 9182740:12, 3296035:13, 9211105:14, 4924360:15,\
           3302560:16, 3988500:17, 3988540:18, 9221340:19, 3302420:20,\
           1326300:21, 3957880:22, 1316060:23, 14464220:24, 14423100:25,\
           11832480:26, 9180240:27, 8204875:28, 10542100:29, 9221140:30,\
           14474380:31, 1351760:32, 6553700:33, 11146310:34, 13145750:35, 2146559:36, 0:0}

    return Map[label_index]
"""

