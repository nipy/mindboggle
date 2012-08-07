==========
Mindboggle
==========

This package enables a user to automate anatomical labeling and shape analysis
of MRI brain image data.

Website
=======

Current information can always be found at the Mindboggle website::

    http://www.mindboggle.info

Documentation
=============

#Mailing Lists
#=============
#
#Please see the developer's list here::
#
#    http://mail.scipy.org/mailman/listinfo/nipy-devel

Code
====

You can find our sources and single-click downloads:

* `Main repository`_ on Github.
* Documentation_ for all releases and current development tree.
* Download as a tar/zip file the `current trunk`_.
* Downloads of all `available releases`_.

.. _main repository: http://github.com/binarybottle/mindboggle
.. _Documentation: http://www.mindboggle.info/documentation
.. _current trunk: http://github.com/binarybottle/mindboggle/archives/master
.. _available releases: http://github.com/binarybottle/mindboggle/downloads

License
=======

Mindboggle is licensed under the terms of the Apache 2.0 license.
Please the COPYING file in the mindboggle distribution.

Structure
=========
In the directory tree below, ? stands for multiple entries,
such as ``?h'' for left and right hemispheres: [lh, rh]

CHANGES      # version information
COPYING      # copyright and license information
INSTALL.rst  # installation instructions
MANIFEST.in  # used for installation
THANKS       # acknowledgments
README.rst   # this file
mindboggle/
    pipeline.py          # main program (nipype pipeline)
    atlas_functions.py   # functions for multi-atlas surface registration/labeling
    surface_functions.py # functions for computing values on surface meshes
    volume_functions.py  # functions for label filling a volume
    utils/
data/
    templates/
        templates_freesurfer/  # FreeSurfer-ready templates
    atlases/
        labels.txt  # list of atlas labels
        info/  # atlas scanning and participant information
        subjects/
            <subject>/
                originalT1.nii.gz
                ?h.pial.[depth, curv…].vtk
                ?.inflated.vtk
                labels.manual.vtk
                labels.manual.nii.gz
                labels.ants.subcortex.nii.gz
                labels.freesurfer.subcortex.nii.gz
        subjects_freesurfer/ # (see below)
            <subject>/
                label/
                    ?h.labels.manual.annot  # manual label files
                surf/
                    ?h.[pial, inflated]  # surfaces
                    ?h.sphere_to_?_template.reg  # template transforms
