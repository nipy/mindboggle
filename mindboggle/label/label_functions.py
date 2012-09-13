#!/usr/bin/python

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

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

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
    annot_name = classifier_atlas + '.annot'
    annot_file = os.path.join(subjects_path, subject, 'label', hemi + '.' + annot_name)
    cli = CommandLine(command='mris_ca_label')
    cli.inputs.args = ' '.join([subject, hemi, sphere_file,
                                classifier_file, annot_file])
    logger.info(cli.cmdline)
    cli.run()

    return annot_name, annot_file