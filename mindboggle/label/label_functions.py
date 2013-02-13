#!/usr/bin/env python

"""
Functions for labeling brains.

Notes regarding creation and use of a FreeSurfer Gaussian classifier atlas:

Create the DKT classifier atlas (?h.DKTatlas40.gcs) --NEED TO VERIFY THIS:
$ mris_ca_train -t $FREESURFERHOME/average/colortable_desikan_killiany.txt \
                $hemi sphere.reg aparcNMMjt.annot $SCANS ./$hemi.DKTatlas40.gcs

Label a brain with the DKT atlas (surface annotation file ?h.DKTatlas40.annot):
$ mris_ca_label -l ./$x/label/$hemi.cortex.label $x/ $hemi sphere.reg \
                ./$hemi.DKTatlas40.gcs ./$x/label/$hemi.DKTatlas40.annot

Label the cortex of a subject's segmented volume according
to the edited surface labels (?h.aparcNMMjt.annot):
$ mri_aparc2aseg --s ./x --volmask --annot aparcNMMjt

Label a brain with the DKT atlas using FreeSurfer's mris_ca_label:
$ mris_ca_label MMRR-21-1 lh lh.sphere.reg ../lh.DKTatlas40.gcs ../out.annot


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com
    - Jason Tourville  (jtour@bu.edu)
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#===============================================================================
# Label with Gaussian classifier atlas
#===============================================================================
def label_with_classifier(hemi, subject, subjects_path, sphere_file,
                          classifier_path, classifier_atlas):
    """
    Label a brain with the DKT atlas using FreeSurfer's mris_ca_label::

        SYNOPSIS
        mris_ca_label [options] <subject> <hemi> <canonsurf> <classifier> <outputfile>

        DESCRIPTION
        For a single subject, produces an annotation file, in which each
        cortical surface vertex is assigned a neuroanatomical label.This
        automatic procedure employs data from a previously-prepared atlas
        file. An atlas file is created from a training set, capturing region
        data manually drawn by neuroanatomists combined with statistics on
        variability correlated to geometric information derived from the
        cortical model (sulcus and curvature). Besides the atlases provided
        with FreeSurfer, new ones can be prepared using mris_ca_train).

    Parameters
    ----------
    hemi : ``string``: hemisphere
    subject : ``string``: subject, corresponding to FreeSurfer subject directory
    subjects_path : ``string``: FreeSurfer subjects directory
    sphere_file : ``string``: FreeSurfer spherical surface
    classifier_path : ``string``: FreeSurfer classifier atlas parent directory
    classifier_atlas : ``string``: name of FreeSurfer classifier atlas (no hemi)

    Returns
    -------
    annot_file : ``string``: name of output .annot file

    Examples
    --------
    $ mris_ca_label MMRR-21-1 lh sphere ../lh.DKTatlas40.gcs ../out.annot
    """
    import os
    from nipype.interfaces.base import CommandLine
    from nipype import logging
    logger = logging.getLogger('interface')

    classifier_file = os.path.join(classifier_path, hemi + '.' + classifier_atlas)
    annot_name = classifier_atlas
    annot_file = os.path.join(subjects_path, subject, 'label',
                              hemi + '.' + annot_name + '.annot')
    cli = CommandLine(command='mris_ca_label')
    cli.inputs.args = ' '.join([subject, hemi, sphere_file,
                                classifier_file, annot_file])
    logger.info(cli.cmdline)
    cli.run()

    return annot_name, annot_file

#===============================================================================
# Extract label borders
#===============================================================================
def extract_borders(labels_file, mask_file='', values_file=''):
    """
    Extract borders (between labels) on a surface.
    Options: Mask out values; extract border values on a second surface.

    Parameters
    ----------
    labels_file : string
        file name for surface mesh with labels
    mask_file : string
        file name for surface mesh with mask (>-1) values
    values_file : string
        file name for surface mesh with values to extract along borders

    Returns
    -------
    border_file : string
        file name for surface mesh with label borders (-1 background values)
    border_values : numpy array
        values for all vertices (-1 for vertices not along label borders)

    Examples
    --------
    >>> import os
    >>> from mindboggle.extract.extract_folds import extract_borders
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>>
    >>> # Load labels, folds, neighbor lists, and sulcus names and label pairs
    >>> labels_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                            'labels', 'lh.labels.DKT25.manual.vtk')
    >>> mask_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                     'features', 'lh.sulci.vtk')
    >>> values_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                       'measures', 'lh.pial.depth.vtk')
    >>> border_file, border_values = extract_borders(labels_file,
    >>>                                              mask_file, values_file)
    >>> os.system('mayavi2 -m Surface -d ' + border_file + ' &')

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    from mindboggle.utils.mesh_operations import find_neighbors, detect_boundary_indices

    # Load labeled surface file
    points, faces, labels, n_vertices = load_scalars(labels_file, True)

    # Detect boundaries
    neighbor_lists = find_neighbors(faces, n_vertices)
    indices_boundaries = detect_boundary_indices(range(len(points)),
        labels, neighbor_lists)

    # Filter values with label boundaries
    border_values = -1 * np.ones(len(points))
    if values_file:
        points, faces, values, n_vertices = load_scalars(values_file, True)
        border_values[indices_boundaries] = values[indices_boundaries]
    else:
        border_values[indices_boundaries] = 1

    # Mask values (for mask >-1)
    if mask_file:
        points, faces, mask_values, n_vertices = load_scalars(mask_file, True)
    else:
        mask_values = []

    # Write out label boundary vtk file
    border_file = os.path.join(os.getcwd(),
                               'borders_' + os.path.basename(labels_file))
    rewrite_scalar_lists(labels_file, border_file,
                         border_values, 'label_borders_in_mask', mask_values)

    return border_file, border_values


def find_superset_subset_lists(labels, label_lists):
    """
    Find *label_lists* that are supersets or subsets of *labels*.

    Parameters
    ----------
    labels : list of integers
        label numbers
    label_lists : list of lists of integers
        each list contains label numbers

    Returns
    -------
    superset_indices : list of integers
        indices to label_lists that are a superset of labels
    subset_indices : list of integers
        indices to label_lists that are a subset of labels

    Example
    -------
    >>> find_superset_subset_lists([1,2],[[1,2],[3,4]])
    >>> [0]
    >>> find_superset_subset_lists([1,2],[[1,2,3],[1,2,5]])
    >>> [0, 1]
    >>> find_superset_subset_lists([1,2],[[2,3],[1,2,5]])
    >>> [1]

    """

    labels = set(labels)
    superset_indices = []
    subset_indices = []
    for Id, label_list in enumerate(label_lists):
        if labels.issubset(set(label_list)):
            superset_indices.append(Id)
        if set(label_list).issubset(labels):
            subset_indices.append(Id)

    return superset_indices, subset_indices
