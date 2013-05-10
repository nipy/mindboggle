#!/usr/bin/env python
"""
Atlas-based functions for FreeSurfer surface registration-based labeling:

1. Gaussian classifier atlas-based labeling
2. Template registration, multi-atlas-based labeling


Authors:
- Arno Klein, 2012-2013 (arno@mindboggle.info) http://binarybottle.com

Copyright 2013, Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
# Gaussian classifier atlas-based labeling
#=============================================================================
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

    Notes regarding creation and use of FreeSurfer Gaussian classifier atlas:

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

#=============================================================================
# Template registration, multi-atlas-based labeling
#=============================================================================

def register_template(hemi, sphere_file, transform,
                      templates_path, template):
    """
    Register surface to template with FreeSurfer's mris_register.

    Transform the labels from multiple atlases via a template
    (using FreeSurfer's mris_register).

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
    Transform atlas labels.

    For each brain hemisphere (left and right) in a given subject,
    read in FreeSurfer *.annot files (multiple labelings), output one VTK file
    of majority vote labels, representing a "maximum probability" labeling.
    The main function is majority_vote() and calls vote_labels().

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
    hemi : string
        hemisphere ['lh' or 'rh']
    subject : string
        subject corresponding to FreeSurfer subject directory
    transform : string
        name of FreeSurfer spherical surface registration transform file
    subjects_path : string
        name of FreeSurfer subjects directory
    atlas : string
        name of atlas
    atlas_string : string
        name of atlas labeling protocol

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
