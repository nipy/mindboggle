==============================================================================
Software
==============================================================================
.. role:: red

The Mindboggle software package automates shape analysis of anatomical labels
and features extracted from human brain magnetic resonance image data.
Mindboggle can be run as a single command, and can be installed as a
cross-platform Docker container for convenience and reproducibility of results.
Behind the scenes, open source Python 3 and C++ code run within a modular
Nipype pipeline framework on Linux (tested with Python 3.5.1 on Ubuntu 14.04).

..
    1. Installing Mindboggle
    2. Running Mindboggle
    3. Preprocessing
    4. Processing steps
    5. Output

:Release: |version|
:Date: |today|

Links:

.. toctree::
    :maxdepth: 1

    FAQ <faq.rst>
    license

- `GitHub <http://github.com/nipy/mindboggle>`_
- `Contributors <http://mindboggle.info/people.html>`_

* :ref:`modindex`
* :ref:`genindex`

------------------------------------------------------------------------------
Installing Mindboggle
------------------------------------------------------------------------------
We recommend installing Mindboggle as a cross-platform Docker container 
for greater convenience and reproducibility of results.
Alternatively, Mindboggle can be installed from scratch on a Linux machine 
using the installation script `install_mindboggle.sh <https://raw.githubusercontent.com/nipy/mindboggle/master/install/install_mindboggle.sh>`_.

1. Install and run Docker on your (macOS, Linux, or Windows) host machine:  

    https://docs.docker.com/engine/installation/

2. Clone the Mindboggle Docker app and change to its directory::

    git clone https://github.com/BIDS-Apps/mindboggle
    cd mindboggle

3. Set the path on your host machine for the Docker container to access
inputs and outputs (the root directory in this example), and enter the
container's bash shell (backslash denotes line return)::

    PATH_ON_HOST=/
    docker run --rm -ti -v $PATH_ON_HOST:/root/data \
        --entrypoint /bin/bash bids/mindboggle 

------------------------------------------------------------------------------
Running Mindboggle
------------------------------------------------------------------------------
Mindboggle takes as its input preprocessed brain MR image data
(see `Preprocessing`_ below). To get up and running with the following
examples, download and unzip the
`mindboggle_input_example.zip <https://osf.io/3xfb8/?action=download&version=1>`_
directory (455 MB), which contains some example preprocessed data.
More example input and output data can be found
on Mindboggle's `examples <https://osf.io/8cf5z>`_ osf.io site.

**Set paths:**

For brevity in the commands below, we set path variables to the example data
and Mindboggle working and output directories::

    HOST=/root/data/Users/arno
    FREESURFER_SUBJECT=$HOST/mindboggle_input_example/freesurfer/subjects/arno
    MINDBOGGLING=$HOST/mindboggle_working
    MINDBOGGLED=$HOST/mindboggled

Additional paths optionally access ANTs commands (if ANTs is installed)
and ANTs preprocessed data::

    ANTSPATH=$HOST/ants/bin
    ANTS_SUBJECT=$HOST/mindboggle_input_example/ants/subjects/arno

**Help and options:**

For help with the Mindboggle command and to learn about all of its
command-line options, type the following in a terminal window::

    mindboggle -h

**Example 1:**
The following bare-bones command runs Mindboggle
on data processed by FreeSurfer but not ANTs (backslash denotes line return)::

    mindboggle $FREESURFER_SUBJECT \
        --working $MINDBOGGLING --out $MINDBOGGLED

**Example 2:**
Same as #1, but takes advantage of ANTs output::

    mindboggle $FREESURFER_SUBJECT \
        --working $MINDBOGGLING --out $MINDBOGGLED \
        --ants $ANTS_SUBJECT/antsBrainSegmentation.nii.gz

**Example 3:**
To generate only volume (and not surface) labels and shape measures
from FreeSurfer and ANTs data, using 8 processors::

    mindboggle $FREESURFER_SUBJECT/arno --no_surfaces -p 8 \
        --working $MINDBOGGLING --out $MINDBOGGLED \
        --ants $ANTS_SUBJECT/antsBrainSegmentation.nii.gz

------------------------------------------------------------------------------
Preprocessing
------------------------------------------------------------------------------
As you may have inferred from the "Running Mindboggle" examples and example
data above, Mindboggle currently takes output from
`FreeSurfer <http://surfer.nmr.mgh.harvard.edu>`_ (v6 or higher recommended)
and optionally from `ANTs <http://stnava.github.io/ANTs/>`_
(v2.1.0rc3 or higher recommended).
Example input data can be found
on Mindboggle's `examples <https://osf.io/8cf5z>`_ site on osf.io.

**FreeSurfer** generates labeled cortical surfaces, and labeled cortical and
noncortical volumes. Run ``recon-all`` on a T1-weighted $IMAGE file
(e.g., subject1.nii.gz) and set the output $SUBJECT name (e.g., subject1)::

    recon-all -all -i $IMAGE -s $SUBJECT

**ANTs** provides brain volume extraction, segmentation, and
registration-based labeling. To generate the ANTs transforms and segmentation
files used by Mindboggle, run the ``antsCorticalThickness.sh`` script on the
same $IMAGE file, set an output $PREFIX, and provide paths to the
`OASIS-30 Atropos template <https://osf.io/rh9km/?action=download&version=1>`_
files (backslash denotes a line return)::

    antsCorticalThickness.sh -d 3 -a $IMAGE -o $PREFIX \
      -e OASIS-30_Atropos_template/T_template0.nii.gz \
      -t OASIS-30_Atropos_template/T_template0_BrainCerebellum.nii.gz \
      -m OASIS-30_Atropos_template/T_template0_BrainCerebellumProbabilityMask.nii.gz \
      -f OASIS-30_Atropos_template/T_template0_BrainCerebellumExtractionMask.nii.gz \
      -p OASIS-30_Atropos_template/Priors2/priors%d.nii.gz

------------------------------------------------------------------------------
Processing steps
------------------------------------------------------------------------------
The following steps are performed by Mindboggle (with links to code on GitHub):

1. Create hybrid gray/white segmentation from FreeSurfer and ANTs output (`combine_2labels_in_2volumes <https://github.com/nipy/mindboggle/blob/master/mindboggle/guts/segment.py>`_).
2. Fill hybrid segmentation with FreeSurfer- or ANTs-registered labels.
3. Compute volume shape measures for each labeled region:

    - volume (`volume_per_brain_region <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/volume_shapes.py>`_)
    - thickness of cortical labels (`thickinthehead <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/volume_shapes.py>`_)

4. Compute surface shape measures for every cortical mesh vertex:

    - `surface area <https://github.com/nipy/mindboggle/blob/master/vtk_cpp_tools/PointAreaComputer.cpp>`_
    - `travel depth <https://github.com/nipy/mindboggle/blob/master/vtk_cpp_tools/TravelDepth.cpp>`_
    - `geodesic depth <https://github.com/nipy/mindboggle/blob/master/vtk_cpp_tools/geodesic_depth/GeodesicDepthMain.cpp>`_
    - `mean curvature <https://github.com/nipy/mindboggle/blob/master/vtk_cpp_tools/curvature/CurvatureMain.cpp>`_
    - convexity (from FreeSurfer)
    - thickness (from FreeSurfer)

5. Extract cortical surface features:

    - `folds <https://github.com/nipy/mindboggle/blob/master/mindboggle/features/folds.py>`_
    - `sulci <https://github.com/nipy/mindboggle/blob/master/mindboggle/features/sulci.py>`_
    - `fundi <https://github.com/nipy/mindboggle/blob/master/mindboggle/features/fundi.py>`_

6. For each cortical surface label/sulcus, compute:

    - `area <https://github.com/nipy/mindboggle/blob/master/vtk_cpp_tools/area/PointAreaMain.cpp>`_
    - mean coordinates: `means_per_label <https://github.com/nipy/mindboggle/blob/master/mindboggle/guts/compute.py>`_
    - mean coordinates in MNI152 space
    - `Laplace-Beltrami spectrum <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/laplace_beltrami.py>`_
    - `Zernike moments <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/zernike/zernike.py>`_

7. Compute statistics (``stats_per_label`` in `compute.py <https://github.com/nipy/mindboggle/blob/master/mindboggle/guts/compute.py>`_) for each shape measure in #4 for each label/feature:

    - median
    - median absolute deviation
    - mean
    - standard deviation
    - skew
    - kurtosis
    - lower quartile
    - upper quartile

------------------------------------------------------------------------------
Output
------------------------------------------------------------------------------
Example output data can be found
on Mindboggle's `examples <https://osf.io/8cf5z>`_ site on osf.io.
By default, output files are saved in $HOME/mindboggled/SUBJECT, where $HOME
is the home directory and SUBJECT is a name representing the person's
brain that has been scanned.
Volume files are in `NIfTI <http://nifti.nimh.nih.gov>`_ format,
surface meshes in `VTK <http://www.vtk.org/>`_ format,
and tables are comma-delimited.
Each file contains integers that correspond to anatomical :doc:`labels <labels>`
or features (0-24 for sulci).
All output data are in the original subject's space.
The following include outputs from most, but not all, optional arguments.

+----------------+----------------------------------------------------+--------------+
|   **Folder**   | **Contents**                                       | **Format**   |
+----------------+----------------------------------------------------+--------------+
|    labels/     |  number-labeled surfaces and volumes               | .vtk, .nii.gz|
+----------------+----------------------------------------------------+--------------+
|    features/   |  surfaces with features:  sulci, fundi             | .vtk         |
+----------------+----------------------------------------------------+--------------+
|    shapes/     |  surfaces with shape measures (per vertex)         | .vtk         |
+----------------+----------------------------------------------------+--------------+
|    tables/     |tables of shape measures (per label/feature/vertex) | .csv         |
+----------------+----------------------------------------------------+--------------+

**mindboggled** / $SUBJECT /

    **labels** /

        **freesurfer_wmparc_labels_in_hybrid_graywhite.nii.gz**:  *hybrid segmentation filled with FS labels*

        **ants_labels_in_hybrid_graywhite.nii.gz**:  *hybrid segmentation filled with ANTs + FS cerebellar labels*

        [left,right]_cortical_surface / **freesurfer_cortex_labels.vtk**: `DKT <http://mindboggle.info/data.html>`_ *cortical surface labels*

    **features** / [left,right]_cortical_surface /

            **folds.vtk**:  *(unidentified) depth-based folds*

            **sulci.vtk**:  *sulci defined by* `DKT <http://mindboggle.info/data.html>`_ *label pairs in depth-based folds*

            **fundus_per_sulcus.vtk**:  *fundus curve per sulcus*  **-- UNDER EVALUATION --**

            **cortex_in_MNI152_space.vtk**:  *cortical surfaces aligned to an MNI152 template*

    **shapes** / [left,right]_cortical_surface /

            **area.vtk**:  *per-vertex surface area*

            **mean_curvature.vtk**:  *per-vertex mean curvature*

            **geodesic_depth.vtk**:  *per-vertex geodesic depth*

            **travel_depth.vtk**:  *per-vertex travel depth*

            **freesurfer_curvature.vtk**:  *FS curvature files converted to VTK*

            **freesurfer_sulc.vtk**:  *FS sulc (convexity) files converted to VTK*

            **freesurfer_thickness.vtk**:  *FS thickness files converted to VTK*

    **tables** /

        **volume_per_freesurfer_label.csv**:  *volume per FS label*

        **volumes_per_ants_label.csv**:  *volume per ANTs label*

        **thickinthehead_per_freesurfer_cortex_label.csv**:  *FS cortex label thickness*

        **thickinthehead_per_ants_cortex_label.csv**:  *ANTs cortex label thickness*

        [left,right]_cortical_surface /

            **label_shapes.csv**:  *per-label surface shape statistics*

            **sulcus_shapes.csv**:  *per-sulcus surface shape statistics*

            **fundus_shapes.csv**:  *per-fundus surface shape statistics*  **-- UNDER EVALUATION --**

            **vertices.csv**:  *per-vertex surface shape statistics*

