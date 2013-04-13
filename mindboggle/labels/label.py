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
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
    hemi : string
        hemisphere ['lh' or 'rh']
    subject : string
        subject corresponding to FreeSurfer subject directory
    subjects_path : string
        name of FreeSurfer subjects directory
    sphere_file : string
        name of FreeSurfer spherical surface file
    classifier_path : string
        name of FreeSurfer classifier atlas parent directory
    classifier_atlas : string
        name of FreeSurfer classifier atlas (no hemi)

    Returns
    -------
    annot_file : string
        name of output .annot file

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

#------------------------------------------------------------------------------
# Extract borders between labels
#------------------------------------------------------------------------------
def extract_borders(indices, labels, neighbor_lists,
                    ignore_values=[], return_label_pairs=False):
    """
    Detect the label boundaries in a collection of vertices such as a region.

    Label boundaries are the set of all vertices
    whose neighbors do not share the same label.

    Parameters
    ----------
    indices : list of integers
        indices to (a subset of) vertices
    labels : numpy array of integers
        label numbers for all vertices (-1 for unlabeled vertices)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    ignore_values : list of integers
        integers to ignore (e.g., background)

    Returns
    -------
    boundary_indices : list of integers
        indices to label boundary vertices
    boundary_label_pairs : list of lists of sorted pairs of integers
        sorted label pairs
    unique_boundary_label_pairs : list of sorted pairs of integers
        unique, sorted label pairs

    Examples
    --------
    >>> # Small example:
    >>> from mindboggle.labels.label import extract_borders
    >>> indices = [0,1,2,4,5,8,9]
    >>> labels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, -1, -1]
    >>> neighbor_lists = [[1,2,3], [1,2], [2,3], [2], [4,7], [3,2,3]]
    >>> extract_borders(indices, labels, neighbor_lists, [], True)
        ([1, 2, 4, 5],
         [[20, 30], [30, 40], [50, 80], [30, 40]],
         [[20, 30], [30, 40], [50, 80]])
    >>> # Real example -- extract sulcus label boundaries:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.label import extract_borders
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(labels_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> indices_boundaries, label_pairs, foo = extract_borders(range(npoints),
    >>>     labels, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> IDs = -1 * np.ones(npoints)
    >>> IDs[indices_boundaries] = 1
    >>> rewrite_scalars(labels_file, 'extract_borders.vtk',
    >>>                 IDs, 'boundaries', IDs)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('extract_borders.vtk')

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if not isinstance(labels, np.ndarray):
        labels = np.array(labels)

    # Construct a list of labels corresponding to the neighbor lists
    label_lists = [list(set(labels[lst])) for lst in neighbor_lists]

    # Find indices to sets of two labels
    boundary_indices = [i for i,x in enumerate(label_lists)
                        if len(set(x)) == 2
                        if i in indices]

    if return_label_pairs:
        boundary_label_pairs = [np.sort(x).tolist()
                                for i,x in enumerate(label_lists)
                                if len(set(x)) == 2
                                if i in indices]
    else:
        boundary_label_pairs = []

    if ignore_values:
        Ikeep = [i for i,x in enumerate(boundary_label_pairs)
                 if not len(frozenset(x).intersection(ignore_values))]
        boundary_indices = [x for i,x in enumerate(boundary_indices)
                            if i in Ikeep]
        if return_label_pairs:
            boundary_label_pairs = [np.sort(x).tolist()
                                    for i,x in enumerate(boundary_label_pairs)
                                    if i in Ikeep]

    if return_label_pairs:
        unique_boundary_label_pairs = []
        for pair in boundary_label_pairs:
            if pair not in unique_boundary_label_pairs:
                unique_boundary_label_pairs.append(pair)
    else:
        unique_boundary_label_pairs = []

    return boundary_indices, boundary_label_pairs, unique_boundary_label_pairs

#===============================================================================
# Extract label borders
#===============================================================================
def extract_border_values(labels_file, mask_file='', values_file=''):
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
    >>> # Extract depth values along label boundaries in sulci (mask):
    >>> import os
    >>> from mindboggle.labels.label import extract_borders
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> mask_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> values_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> #
    >>> border_file, border_values = extract_borders(labels_file, mask_file, values_file)
    >>> #
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk(border_file)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.labels.label import extract_borders

    # Load labeled surface file
    faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(labels_file,
                                                                    return_first=True, return_array=True)

    # Detect boundaries
    neighbor_lists = find_neighbors(faces, npoints)
    indices_boundaries, label_pairs, foo = extract_borders(range(npoints),
                                                           labels, neighbor_lists)

    # Filter values with label boundaries
    border_values = -1 * np.ones(npoints)
    if values_file:
        values, name = read_scalars(values_file, return_first=True, return_array=True)
        border_values[indices_boundaries] = values[indices_boundaries]
    else:
        border_values[indices_boundaries] = 1

    # Mask values (for mask >-1)
    if mask_file:
        mask_values, name = read_scalars(mask_file)
    else:
        mask_values = []

    # Write out label boundary vtk file
    border_file = os.path.join(os.getcwd(), 'borders_' + os.path.basename(labels_file))
    rewrite_scalars(labels_file, border_file, border_values, \
                    'label_borders_in_mask', mask_values)

    return border_file, border_values
