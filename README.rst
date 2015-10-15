.. _README:

==============================================================================
Software -- NOTE: MINDBOGGLE IS STILL IN PRE-RELEASE AND IS UNDERGOING CHANGES
==============================================================================
| 1. `Introduction and help`_
| 2. `Set up Mindboggle`_
| 3. `Example Mindboggle commands`_
| 4. `FreeSurfer and ANTs preprocessing`_
| 5. `Mindboggle processing steps`_
| 6. `Mindboggle output files`_

------------------------------------------------------------------------------
_`Introduction and help`
------------------------------------------------------------------------------
The Mindboggle software package automates shape analysis of anatomical labels
and features extracted from human brain MR image data (read the
`story <http://mindboggle.info/faq/why_mindboggle.html>`_).
Mindboggle can be run as a single command, and can be
easily installed as a cross-platform virtual machine for convenience and
reproducibility of results. Behind the scenes, open source
Python and C++ code run within a Nipype pipeline framework.

- For help in a terminal window (see below for inputs and outputs)::

    mindboggle -h

- `FAQs <http://www.mindboggle.info/faq/>`_
- `Example output data <http://media.mindboggle.info/data/examples/arno.tar.gz>`_
- `Documentation <http://mindboggle.info/documentation.html>`_
- `Installation <http://mindboggle.info/users/INSTALL.html>`_
- `GitHub <http://github.com/binarybottle/mindboggle>`_
- `License <http://mindboggle.info/users/LICENSE.html>`_
- `Contributors <http://mindboggle.info/users/THANKS.html>`_

------------------------------------------------------------------------------
_`Set up Mindboggle`
------------------------------------------------------------------------------
If running Mindboggle in a virtual machine (recommended),
type the following two commands in a terminal window,
in the same directory as the Vagrantfile
you generated (see `INSTALL <http://mindboggle.info/users/INSTALL.html>`_). This will launch and log into
the Mindboggle virtual machine (requires an active Internet connection)::

    vagrant up
    vagrant ssh

------------------------------------------------------------------------------
_`Example Mindboggle commands`
------------------------------------------------------------------------------
**Example 1:**
The following bare-bones command runs Mindboggle (and its dependencies)
on data processed by FreeSurfer but not ANTs
(replace ``SUBJECT`` with the name of a FreeSurfer subject directory, such as ``bert``)::

    mindboggle SUBJECT

**Example 2:**
To generate only volume data (no surface labels or measures),
this command uses ANTs output files
(replace ``SEGMENTS`` with an ANTs segmented file, such as
``ants_output/subject1/antsBrainSegmentation.nii.gz``)::

    mindboggle SUBJECT --ants SEGMENTS --no_surfaces

**Example 3:**
To compute all shape measures on all labels and features using 8 processors
(type ``mindboggle --help`` for more options)::

    mindboggle SUBJECT --ants SEGMENTS --all -p 8

------------------------------------------------------------------------------
_`FreeSurfer and ANTs preprocessing`
------------------------------------------------------------------------------
Mindboggle currently takes output from `FreeSurfer <http://surfer.nmr.mgh.harvard.edu>`_
(preferably v5.3 and above) and the latest `ANTs <http://stnava.github.io/ANTs/>`_ software packages.

**FreeSurfer** generates labeled cortical surfaces, and labeled cortical and
noncortical volumes. Run ``recon-all`` on a T1-weighted ``IMAGE`` file
(e.g., ``subject1.nii.gz``) and set the output ``SUBJECT`` name (e.g., to ``subject1``),
while calling the DKT40 cortical surface atlas to aid in cortical labeling::

    recon-all -all -i IMAGE -s SUBJECT -gcs DKTatlas40.gcs

..
    - mri/orig/001.mgz
    - mri/[wmparc,aparc+aseg].mgz
    - surf/[lh,rh].pial
    - label/[lh,rh].aparc.annot

**ANTs** provides brain volume extraction, segmentation, and registration-based labeling.
To generate the ANTs transforms and segmentation files used by
Mindboggle, run the ``antsCorticalThickness.sh`` script on the same ``IMAGE`` file,
set an output ``PREFIX``, and provide paths to the
`OASIS-30 Atropos template <http://mindboggle.info/data/templates/atropos/OASIS-30_Atropos_template.tar.gz>`_
files (backslash denotes a line return)::
    antsCorticalThickness.sh -d 3 -a IMAGE -o PREFIX \
      -e OASIS-30_Atropos_template/T_template0.nii.gz \
      -t OASIS-30_Atropos_template/T_template0_BrainCerebellum.nii.gz \
      -m OASIS-30_Atropos_template/T_template0_BrainCerebellumProbabilityMask.nii.gz \
      -f OASIS-30_Atropos_template/T_template0_BrainCerebellumExtractionMask.nii.gz \
      -p OASIS-30_Atropos_template/Priors2/priors%d.nii.gz

------------------------------------------------------------------------------
_`Mindboggle processing steps`
------------------------------------------------------------------------------
    1. Create hybrid gray/white segmentation from FreeSurfer and ANTs output (`combine_2labels_in_2volumes <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/guts/segment.py>`_).
    2. Fill hybrid segmentation with FreeSurfer- or ANTs-registered labels.
    3. Compute volume shape measures for each labeled region:

        - volume (`volume_per_label <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/guts/compute.py>`_)
        - thickness of cortical labels (`thickinthehead <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/guts/ants.py>`_)

    4. Compute surface shape measures for every cortical mesh vertex:

        - `surface area <https://github.com/binarybottle/mindboggle/blob/master/surface_cpp_tools/PointAreaComputer.cpp>`_
        - `travel depth <https://github.com/binarybottle/mindboggle/blob/master/surface_cpp_tools/TravelDepth.cpp>`_
        - `geodesic depth <https://github.com/binarybottle/mindboggle/blob/master/surface_cpp_tools/geodesic_depth/GeodesicDepthMain.cpp>`_
        - `mean curvature <https://github.com/binarybottle/mindboggle/blob/master/surface_cpp_tools/curvature/CurvatureMain.cpp>`_
        - convexity (from FreeSurfer)
        - thickness (from FreeSurfer)

    5. Extract cortical surface features:

        - `folds <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/features/folds.py>`_
        - `sulci <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/features/sulci.py>`_
        - `fundi <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/features/fundi.py>`_

    6. For each cortical surface label/sulcus, compute:

        - area
        - mean coordinates
        - mean coordinates in MNI152 space
        - `Laplace-Beltrami spectrum <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/shapes/laplace_beltrami.py>`_
        - `Zernike moments <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/shapes/zernike/zernike.py>`_

    7. Compute statistics (``stats_per_label`` in `compute.py <https://github.com/binarybottle/mindboggle/blob/master/mindboggle/guts/compute.py>`_) for each shape measure in #4 for each label/feature:

        - median
        - median absolute deviation
        - mean
        - standard deviation
        - skew
        - kurtosis
        - lower quartile
        - upper quartile

------------------------------------------------------------------------------
_`Mindboggle output files`
------------------------------------------------------------------------------
Example output data can be downloaded from 
`here <http://media.mindboggle.info/data/examples/arno.tar.gz>`_.
By default, output files are saved in HOME/mindboggled/SUBJECT,
where HOME is the home directory and SUBJECT is the name of the subject.
Volume files are in `Nifti <http://nifti.nimh.nih.gov>`_ format,
surface meshes in `VTK <http://www.vtk.org/>`_ format,
and tables are comma-delimited.
Each file contains integers that correspond to anatomical
`labels <http://mindboggle.info/faq/labels.html>`_
or features (e.g., 0-24 for sulci).
All output data are in the original subject's space.
The following include outputs from most, but not all, optional arguments.

+-+---------------+----------------------------------------------------+--------------+
| |  **Folder**   | **Contents**                                       | **Format**   |
+-+---------------+----------------------------------------------------+--------------+
| |   labels/     |  number-labeled surfaces and volumes               | .vtk, .nii.gz|
+-+---------------+----------------------------------------------------+--------------+
| |   features/   |  surfaces with features:  sulci, fundi             | .vtk         |
+-+---------------+----------------------------------------------------+--------------+
| |   shapes/     |  surfaces with shape measures (per vertex)         | .vtk         |
+-+---------------+----------------------------------------------------+--------------+
| |   tables/     |tables of shape measures (per label/feature/vertex) | .csv         |
+-+---------------+----------------------------------------------------+--------------+

**mindboggled** / SUBJECT /

    **labels** /

        **FreeSurfer_wmparc_filled_labels.nii.gz**:  *hybrid segmentation filled with FS labels*

        **ANTs_filled_labels.nii.gz**:  *hybrid segmentation filled with ANTs + FS cerebellar labels*

        [left,right]_surface / **FreeSurfer_cortex_labels.vtk**:  *FS or* `DKT <http://mindboggle.info/data/>`_ *cortical surface labels*

    **features** / [left,right]_surface /

            **sulci.vtk**:  *sulci defined by* `DKT <http://mindboggle.info/data/>`_ *label pairs in depth-based folds*

            **fundus_per_sulcus.vtk**:  *fundus curve per sulcus*  **-- UNDER EVALUATION --**

    **shapes** / [left,right]_surface /

            **area.vtk**:  *per-vertex surface area*

            **mean_curvature.vtk**:  *per-vertex mean curvature*

            **geodesic_depth.vtk**:  *per-vertex geodesic depth*

            **travel_depth.vtk**:  *per-vertex travel depth*

            **FreeSurfer_convexity.vtk**:  *FS sulc files converted to VTK*

            **FreeSurfer_thickness.vtk**:  *FS thickness files converted to VTK*

    **tables** /

        **volumes_FreeSurfer_labels.csv**:  *volume per FS-filled label*

        **volumes_ANTs_labels.csv**:  *volume per ANTs-filled label*

        **thickinthehead_FreeSurfer_labels.csv**:  *thickness measure per FS-filled cortical label*

        **thickinthehead_ANTs_labels.csv**:  *thickness measure per ANTs-filled cortical label*

        [left,right]_surface /

            **label_shapes.csv**:  *per-label surface shape statistics*

            **sulcus_shapes.csv**:  *per-sulcus surface shape statistics*

            **fundus_shapes.csv**:  *per-fundus surface shape statistics*  **-- UNDER EVALUATION --**

            **vertices.csv**:  *per-vertex surface shape statistics*
