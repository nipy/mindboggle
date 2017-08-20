==============================================================================
Software
==============================================================================
.. role:: red

Mindboggle's open source brain morphometry platform takes in preprocessed T1-weighted
MRI data, and outputs volume, surface, and tabular data containing label, feature, and shape
information for further analysis. Mindboggle can be run on the command line as "mindboggle"
and also exists as a cross-platform Docker container for convenience and reproducibility
of results. The software runs on Linux and is written in Python 3 and Python-wrapped C++ code 
called within a `Nipype <http://nipy.org/nipype>`_ pipeline framework. 
We have tested the software most extensively with Python 3.5.1 on Ubuntu Linux 14.04.

..
    1. Reference
    2. Getting help
    3. Installation
    4. Run one command
    5. Run separate commands
    6. Visualize output
    7. Appendix: processing
    8. Appendix: output

:Release: |version|
:Date: |today|

Links:

.. toctree::
    :maxdepth: 1

    FAQ <faq.rst>
    license

- `GitHub <http://github.com/nipy/mindboggle>`_ and `Circleci tests <https://circleci.com/gh/nipy/mindboggle>`_
- `Contributors <http://mindboggle.info/people.html>`_

* :ref:`modindex`
* :ref:`genindex`

------------------------------------------------------------------------------
Reference
------------------------------------------------------------------------------
A Klein, SS Ghosh, FS Bao, J Giard, Y Hame, E Stavsky, N Lee, B Rossa,
M Reuter, EC Neto, A Keshavan. 2017.
**Mindboggling morphometry of human brains**.
*PLoS Computational Biology* 13(3): e1005350.
`doi:10.1371/journal.pcbi.1005350 <http://dx.doi.org/10.1371/journal.pcbi.1005350>`_

------------------------------------------------------------------------------
Getting help
------------------------------------------------------------------------------
If you have any questions about Mindboggle, please post to
`NeuroStars <https://neurostars.org/tags/mindboggle/>`_ with the tag
"mindboggle". If you have found a bug, big or small, please
`submit an issue <https://github.com/nipy/mindboggle/issues>`_ on GitHub.

To run the mindboggle jupyter notebook tutorial, first install the mindboggle
docker container (see below) and enter the bash shell of the container
from your $HOST (e.g., /Users/arno)::

    docker run --rm -ti -v $HOST:/home/jovyan/work -p 8888:8888 --entrypoint /bin/bash nipy/mindboggle

Then run the notebook from within the container::

    jupyter notebook /opt/mindboggle/docs/mindboggle_tutorial.ipynb

------------------------------------------------------------------------------
Installation
------------------------------------------------------------------------------
We recommend installing Mindboggle and its dependencies as a cross-platform
Docker container for greater convenience and reproducibility of results.
All the examples below assume you are using this Docker container,
with the path /home/jovyan/work/ pointing to your host machine.
(Alternatively, one can `create a Singularity image <http://mindboggle.readthedocs.io/en/latest/faq/singularity.html>`_,
or Mindboggle can be installed from scratch on a Linux machine using
`this script <https://raw.githubusercontent.com/nipy/mindboggle/master/install_mindboggle.sh>`_).

1. `Install and run Docker <https://docs.docker.com/engine/installation/>`_
on your (macOS, Linux, or Windows) host machine.

2. Download the Mindboggle Docker container (copy/paste the following in a
terminal window)::

    docker pull nipy/mindboggle

*Note 1: This contains FreeSurfer, ANTs, and Mindboggle, so it is currently
over 6GB.*

*Note 2: You may need to increase memory allocated by Docker to at least 5GB.
For example: By default, Docker for `Mac is set to use 2 GB runtime memory <https://docs.docker.com/docker-for-mac/>`_.*

3. Optionally download sample data. To try out the ``mindboggle`` examples
below, download and unzip the directory of example input data
`mindboggle_input_example.zip <https://osf.io/3xfb8/?action=download&version=1>`_ (455 MB).
For example MRI data to preprocess with FreeSurfer and ANTs software,
download and unzip
`example_mri_data.zip <https://osf.io/k3m94/?action=download&version=1>`_ (29 MB).

4. Optionally set environment variables for clarity in the commands below
(modify accordingly, except for DOCK)::

    HOST=/Users/binarybottle  # path on host to access input/output
    DOCK=/home/jovyan/work  # path to HOST from Docker container
    IMAGE=$DOCK/example_mri_data/T1.nii.gz  # input image on HOST
    ID=arno  # ID for brain image

------------------------------------------------------------------------------
Run one command
------------------------------------------------------------------------------
The Mindboggle Docker container can be run as a single command to process
a T1-weighted MR brain image through FreeSurfer, ANTs, and Mindboggle.
Skip to the next section if you wish to run ``recon-all``,
``antsCorticalThickness.sh``, and ``mindboggle`` differently::

    docker run --rm -ti -v $HOST:$DOCK nipy/mindboggle $IMAGE --id $ID

Outputs are stored in $DOCK/mindboggle123_output/ by default,
but you can set a different output path with ``--out $OUT``.

------------------------------------------------------------------------------
Run separate commands
------------------------------------------------------------------------------
If finer control is needed over the software in the Docker container,
the following instructions outline how to run each command separately.
Mindboggle currently takes output from FreeSurfer and optionally from ANTs.
*FreeSurfer version 6 or higher is recommended because by default it uses
Mindboggle’s DKT-100 surface-based atlas to generate corresponding labels
on the cortical surfaces and in the cortical and non-cortical volumes
(v5.3 generates these surface labels by default; older versions require
"-gcs DKTatlas40.gcs" to generate these surface labels).*

1. Enter the Docker container's bash shell to run ``recon-all``,
``antsCorticalThickness.sh``, and ``mindboggle`` commands::

    docker run --rm -ti -v $HOST:$DOCK --entrypoint /bin/bash nipy/mindboggle

2. `FreeSurfer <http://surfer.nmr.mgh.harvard.edu>`_ generates labeled
cortical surfaces, and labeled cortical and noncortical volumes.
Run ``recon-all`` on a T1-weighted IMAGE file (and optionally a T2-weighted
image), and set the output ID name as well as the $FREESURFER_OUT output
directory::

    FREESURFER_OUT=$DOCK/freesurfer_subjects

    recon-all -all -i $IMAGE -s $ID -sd $FREESURFER_OUT

3. `ANTs <http://stnava.github.io/ANTs/>`_ provides brain volume extraction,
segmentation, and registration-based labeling. ``antsCorticalThickness.sh``
generates transforms and segmentation files used by Mindboggle, and is run
on the same IMAGE file and ID as above, with $ANTS_OUT output directory.
TEMPLATE points to the `OASIS-30_Atropos_template <https://osf.io/rh9km/>`_ folder
already installed in the Docker container ("\\" splits the command for readability)::

    ANTS_OUT=$DOCK/ants_subjects
    TEMPLATE=/opt/data/OASIS-30_Atropos_template

    antsCorticalThickness.sh -d 3 -a $IMAGE -o $ANTS_OUT/$ID/ants \
      -e $TEMPLATE/T_template0.nii.gz \
      -t $TEMPLATE/T_template0_BrainCerebellum.nii.gz \
      -m $TEMPLATE/T_template0_BrainCerebellumProbabilityMask.nii.gz \
      -f $TEMPLATE/T_template0_BrainCerebellumExtractionMask.nii.gz \
      -p $TEMPLATE/Priors2/priors%d.nii.gz

4. **Mindboggle** can be run on data preprocessed by ``recon-all`` and
``antsCorticalThickness.sh`` as above by setting::

    FREESURFER_SUBJECT=$FREESURFER_OUT/$ID
    ANTS_SUBJECT=$ANTS_OUT/$ID
    OUT=$DOCK/mindboggled  # output folder

Or it can be run on the
`mindboggle_input_example <https://osf.io/3xfb8/?action=download&version=1>`_
preprocessed data by setting::

    EXAMPLE=$DOCK/mindboggle_input_example
    FREESURFER_SUBJECT=$EXAMPLE/freesurfer/subjects/arno
    ANTS_SUBJECT=$EXAMPLE/ants/subjects/arno
    OUT=$DOCK/mindboggled  # output folder

**Example Mindboggle commands:**

To learn about Mindboggle's command options, type this in a terminal window::

    mindboggle -h

**Example 1:**
Run Mindboggle on data processed by FreeSurfer but not ANTs::

    mindboggle $FREESURFER_SUBJECT --out $OUT

**Example 2:**
Same as Example 1 with output to visualize surface data with roygbiv::

    mindboggle $FREESURFER_SUBJECT --out $OUT --roygbiv

**Example 3:**
Take advantage of ANTs output as well ("\\" splits for readability)::

    mindboggle $FREESURFER_SUBJECT --out $OUT --roygbiv \
        --ants $ANTS_SUBJECT/antsBrainSegmentation.nii.gz

**Example 4:**
Generate only volume (no surface) labels and shapes::

    mindboggle $FREESURFER_SUBJECT --out $OUT \
        --ants $ANTS_SUBJECT/antsBrainSegmentation.nii.gz \
        --no_surfaces

------------------------------------------------------------------------------
Visualize output
------------------------------------------------------------------------------
To visualize Mindboggle output with roygbiv, start the Docker image with::

    docker run --rm -ti -v $HOST:$DOCK -p 5000:5000 --entrypoint /bin/bash nipy/mindboggle

and then inside the image, run roygbiv on an output directory::

    roygbiv $OUT/$ID

and open a browser to `localhost:5000`.

Right now, roygbiv only shows summarized data, but Anisha Keshavan is working
on by-vertex visualizations (for the latter, try `Paraview <https://www.paraview.org/2>`_).

------------------------------------------------------------------------------
Appendix: processing
------------------------------------------------------------------------------
The following steps are performed by Mindboggle (with links to code on GitHub):

1. Create hybrid gray/white segmentation from FreeSurfer and ANTs output (`combine_2labels_in_2volumes <https://github.com/nipy/mindboggle/blob/master/mindboggle/guts/segment.py#L1660>`_).
2. Fill hybrid segmentation with FreeSurfer- or ANTs-registered labels.
3. Compute volume shape measures for each labeled region:

    - volume (`volume_per_brain_region <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/volume_shapes.py#L14>`_)
    - thickness of cortical labels (`thickinthehead <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/volume_shapes.py#L132>`_)

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
    - mean coordinates: `means_per_label <https://github.com/nipy/mindboggle/blob/master/mindboggle/guts/compute.py#L512>`_
    - mean coordinates in MNI152 space
    - `Laplace-Beltrami spectrum <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/laplace_beltrami.py>`_
    - `Zernike moments <https://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/zernike/zernike.py>`_

7. Compute statistics (``stats_per_label`` in `compute.py <https://github.com/nipy/mindboggle/blob/master/mindboggle/guts/compute.py#L716>`_) for each shape measure in #4 for each label/feature:

    - median
    - median absolute deviation
    - mean
    - standard deviation
    - skewhttps://github.com/nipy/mindboggle/blob/master/mindboggle/shapes/volume_shapes.py#L132
    - kurtosis
    - lower quartile
    - upper quartile

------------------------------------------------------------------------------
Appendix: output
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

