====
Data
====

Except where noted, all data are licensed under a Creative Commons License: |CC_license|_

Manually labeled brain surfaces and volumes
-------------------------------------------

**COMING VERY SOON!** (article under review)

|  The Mindboggle-101 dataset includes manually labeled anatomical regions for `101 healthy subjects`_.
|  The manually edited cortical labels_ follow sulcus landmarks_ according to the Desikan-Killiany-Tourville
|  (DKT) protocol.  A paper detailing this protocol is currently under review and will be posted here soon.
|  Here we will provide labeled nifti volumes (nii), vtk surfaces (vtk), and FreeSurfer files (mgh).

..
  - **MMRR-21** [`nii <http://mindboggle.info/data/atlases/MMRR-21_nii.tar.gz>`_, `vtk <http://mindboggle.info/data/atlases/MMRR-21_vtk.tar.gz>`_, `mgh <http://mindboggle.info/data/atlases/MMRR-21_mgh.tar.gz>`_]:
    All 21 subjects in the Multi-Modal MRI Reproducibility Resource |MMRR www|_
  - **NKI-RS-22** [`nii <http://mindboggle.info/data/atlases/NKI-RS-22_nii.tar.gz>`_, `vtk <http://mindboggle.info/data/atlases/NKI-RS-22_vtk.tar.gz>`_,  `mgh <http://mindboggle.info/data/atlases/NKI-RS-22_mgh.tar.gz>`_]:
    22 subjects from the Nathan Klein Institute / Rockland Sample |NKI-RS www|_
  - **NKI-TRT-20** [`nii <http://mindboggle.info/data/atlases/NKI-TRT-20_nii.tar.gz>`_, `vtk <http://mindboggle.info/data/atlases/NKI-TRT-20_vtk.tar.gz>`_, `mgh <http://mindboggle.info/data/atlases/NKI-TRT-20_mgh.tar.gz>`_]:
    20 subjects from the Nathan Klein Institute / Test-Retest Sample |NKI-TRT www|_
  - **OASIS-TRT-20** [`nii <http://mindboggle.info/data/atlases/OASIS-TRT-20_nii.tar.gz>`_, `vtk <http://mindboggle.info/data/atlases/OASIS-TRT-20_vtk.tar.gz>`_, `mgh <http://mindboggle.info/data/atlases/OASIS-TRT-20_mgh.tar.gz>`_]:
    All 20 subjects from the OASIS Test-Retest sample |OASIS-TRT www|_
      - Subcortex_ (nii):  Separate subcortical regions manually labeled by Neuromorphometrics_: |CC_license_nond|_
  - **Extra-18** [`nii <http://mindboggle.info/data/atlases/HLN_MMRR-3T7T_Colin27_Twins_Afterthought_nii.tar.gz>`_, `vtk <http://mindboggle.info/data/atlases/HLN_MMRR-3T7T_Colin27_Twins_Afterthought_vtk.tar.gz>`_, `mgh <http://mindboggle.info/data/atlases/HLN_MMRR-3T7T_Colin27_Twins_Afterthought_mgh.tar.gz>`_]:
    - **HLN-12**:  All 12 subjects from the Human Language Network study
    - **MMRR-3T7T-2**:  2 subjects acquired like MMRR-21
      - multimodal and 7T scans for subject `1 <data/mindboggle101/MMRR-3T7T-2-1_multimodal.tar.gz>`_ and `2 <data/mindboggle101/MMRR-3T7T-2-2_multimodal.tar.gz>`_ (0.4gb each)
    - **Colin27-1**:  Colin Holmes template (average of 27 scans)
    - **Twins-2**:  2 identical twins, including AK
    - **Afterthought-1**:  1 brain imager, SG

| 
|
.. image:: http://media.mindboggle.info/images/data/DKT_labels_width800px.png
|
|
|
.. image:: http://media.mindboggle.info/images/data/DKT_sulci_width800px.png
|

.. _CC_license: http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US
.. |CC_license| image:: http://i.creativecommons.org/l/by-nc-sa/3.0/80x15.png
.. _`101 healthy subjects`: http://media.mindboggle.info/images/data/Mindboggle101_table.pdf
.. _labels: http://media.mindboggle.info/images/data/DKT_label_table.pdf
.. _landmarks: http://media.mindboggle.info/images/data/DKT_sulci_table.pdf
.. _`MMRR www`: http://www.nitrc.org/projects/multimodal
.. _`NKI-RS www`: http://fcon_1000.projects.nitrc.org/indi/pro/nki.html
.. _`NKI-TRT www`: http://fcon_1000.projects.nitrc.org/indi/pro/eNKI_RS_TRT/FrontPage.html
.. _`OASIS-TRT www`: http://www.oasis-brains.org/app/action/BundleAction/bundle/OAS1_RELIABILITY
.. |MMRR www| image:: images/link-brown-12x12.png
.. |NKI-RS www| image:: images/link-brown-12x12.png
.. |NKI-TRT www| image:: images/link-brown-12x12.png
.. |OASIS-TRT www| image:: images/link-brown-12x12.png
.. _Subcortex: http://mindboggle.info/data/atlases/OASIS-TRT-20_subcortex.tar.gz
.. _Neuromorphometrics: http://neuromorphometrics.com
.. _CC_license_nond: http://creativecommons.org/licenses/by-nc-nd/3.0/deed.en_US
.. |CC_license_nond| image:: http://i.creativecommons.org/l/by-nc-nd/3.0/80x15.png
.. _Extra-18: http://mindboggle.info/data/atlases/HLN_MMRR-3T7T_Colin27_Twins_Afterthought.tar.gz

Mindboggle-101 atlases
----------------------

|  Each of the 101 individually labeled brain surfaces and volumes above is an atlas,
|  a labeled or annotated brain image used for transferring labels to unlabeled brains. 
|  We have combined their labels to create aggregate atlases here as well. 
|  The purpose of registering to atlases is to help give a **rough** anatomical labeling,
|  or to initialize labels for further refinement, as is done by the Mindboggle software.

    - `DKT classifier atlas 101`: FreeSurfer atlas (.gcs) from all 101 Mindboggle-101 participants
    - `DKT classifier atlas 40`_: FreeSurfer atlas (.gcs) from 40 of the Mindboggle-101 participants

.. _`DKT classifier atlas 101`: http://mindboggle.info/data/atlases/classifiers/DKTatlas101.tar.gz
.. _`DKT classifier atlas 40`: http://mindboggle.info/data/atlases/classifiers/DKTatlas40.tar.gz


Mindboggle-101 templates
------------------------

|  A template is an unlabeled image used as a reference or standard, often for registering other images to each other. 
|  Each one of the image volumes and surfaces below was constructed by combining the images from multiple subjects. 
|  ANTS templates were made with buildtemplateparallel.sh_ and FreeSurfer templates with make_freesurfer_template.py_.

  **Brain volumes**: ANTS nonlinear optimal average templates (.nii.gz)

  - `MMRR-21 brain`_ template from 21 brains (2012) 
  - `MMRR-21 to MNI152`_: MMRR-21 template `affine`_ transformed to `MNI152`_ (2012) 
  - `OASIS-TRT-20 brain`_ template from 20 brains (2012)
    
  **Head volumes**: ANTS nonlinear optimal average templates (.nii.gz)

  - `HLN-12 head`_ template from 12 heads (`bet brain <http://mindboggle.info/data/templates/ants/HLN-12_head_template_bet.nii.gz>`_) (2012) 
  - `MMRR-21 head`_ template from 21 heads (2012) 
  - `NKI-RS-22 head`_ template from 22 heads (`bet brain <http://mindboggle.info/data/templates/ants/NKI-RS-22_head_template_bet.nii.gz>`_) (2012) 
  - `NKI-TRT-20 head`_ template from 20 heads (`bet brain <http://mindboggle.info/data/templates/ants/NKI-TRT-20_head_template_bet.nii.gz>`_) (2012) 
  - `OASIS-TRT-20 head`_ template from 20 heads (2012)

  **Cortical surfaces**: FreeSurfer nonlinear optimal average templates (.tif)
    
  - `HLN-12 surface`_ template from 12 brains (2012) 
  - `MMRR-21 surface`_ template from 21 brains (2012) 
  - `NKI-RS-22 surface`_ template from 22 brains (2012) 
  - `NKI-TRT-20 surface`_ template from 20 brains (2012) 
  - `OASIS-TRT-20 surface`_ template from 20 brains (2012)


.. _buildtemplateparallel.sh: data/templates/buildtemplateparallel.sh
.. _make_freesurfer_template.py: data/templates/make_freesurfer_template.txt
.. _`MMRR-21 brain`: http://mindboggle.info/data/templates/ants/MMRR-21_template.nii.gz
.. _`MMRR-21 to MNI152`: http://mindboggle.info/data/templates/ants/MMRR-21_template_to_MNI152.nii.gz
.. _`affine`: http://mindboggle.info/data/templates/ants/MMRR-21_template_to_MNI152_affine.txt
.. _`MNI152`: http://mindboggle.info/data/templates/MNI152_T1_1mm_brain.nii.gz
.. _`OASIS-TRT-20 brain`: http://mindboggle.info/data/templates/ants/OASIS-TRT-20_template.nii.gz
.. _`HLN-12 head`: http://mindboggle.info/data/templates/ants/HLN-12_head_template.nii.gz
.. _`MMRR-21 head`: http://mindboggle.info/data/templates/ants/MMRR-21_head_template.nii.gz
.. _`NKI-RS-22 head`: http://mindboggle.info/data/templates/ants/NKI-RS-22_head_template.nii.gz
.. _`NKI-TRT-20 head`: http://mindboggle.info/data/templates/ants/NKI-TRT-20_head_template.nii.gz
.. _`OASIS-TRT-20 head`: http://mindboggle.info/data/templates/ants/OASIS-20_head_template.nii.gz
.. _`HLN-12 surface`: http://mindboggle.info/data/templates/freesurfer/HLN-12_surface_template.nii.gz
.. _`MMRR-21 surface`: http://mindboggle.info/data/templates/freesurfer/MMRR-21_surface_template.nii.gz
.. _`NKI-RS-22 surface`: http://mindboggle.info/data/templates/freesurfer/NKI-RS-22_surface_template.nii.gz
.. _`NKI-TRT-20 surface`: http://mindboggle.info/data/templates/freesurfer/NKI-TRT-20_surface_template.nii.gz
.. _`OASIS-TRT-20 surface`: http://mindboggle.info/data/templates/freesurfer/OASIS-TRT-20_surface_template.nii.gz


Other templates and manually labeled brains
-------------------------------------------

| The Mindboggle-101 templates and manually labeled brains above benefit from the application
| of a consistent labeling protocol by the same labelers, to reduce variability in label assignments.
| The following manually labeled image volumes used different labeling protocols,
| but have been evaluated for use with registration-based labeling and brain extraction
| (see `2009 evaluation`_ and Atropos_ articles), as have the templates.

  **Atlases**: manually labeled volumes (.nii.gz)

  - CUMC-12_: 12 labeled brains
  - IBSR-18_: 18 labeled brains
  - MGH-10_: 10 labeled brains
  - Atropos-18_: 8-class labeled templates for brain extraction, from 18 subjects

  **Templates**: ANTS nonlinear optimal average templates (.nii.gz)

  - `CUMC12 brain`_ template from 12 brains (2010)
  - `LPBA40 brain`_ template from 40 brains (2011)
  - See Satrajit Ghosh's `pediatric template`_ of 31 brains (2011) 

|

.. image:: http://media.mindboggle.info/images/data/evaluation2009_80atlases.png

.. _`2009 evaluation`: http://www.mindboggle.info/papers/evaluation_NeuroImage2009.php
.. _Atropos: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3297199/
.. _CUMC-12: http://mindboggle.info/papers/evaluation_NeuroImage2009/data/CUMC12.tar.gz
.. _IBSR-18: http://mindboggle.info/papers/evaluation_NeuroImage2009/data/IBSR18.tar.gz
.. _MGH-10: http://mindboggle.info/papers/evaluation_NeuroImage2009/data/MGH10.tar.gz
.. _Atropos-18: http://mindboggle.info/data/templates/Atropos_brain_extraction_template.tar.gz
.. _`CUMC12 brain`: http://mindboggle.info/data/templates/ants/CUMC12_template.nii.gz
.. _`LPBA40 brain`: http://mindboggle.info/data/templates/ants/LPBA40_template.nii.gz
.. _`pediatric template`: http://www.mit.edu/~satra/research/pubdata/index.html
