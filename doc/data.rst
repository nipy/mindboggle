==============================================================================
Data
==============================================================================
| Welcome to the world's largest collection of freel, manually labeled human brain image data!
| Please cite the following article and this website when making use of Mindboggle-101 data:
| `101 labeled brain images and a consistent human cortical labeling protocol`_
| Arno Klein, Jason Tourville. Frontiers in Brain Imaging Methods. 6:171. DOI: 10.3389/fnins.2012.00171
|
| See the `README <http://mindboggle.info/data/mindboggle101/README.txt>`_, `labels <http://mindboggle.info/faq/labels.html>`_, the `CHANGELOG <http://mindboggle.info/data/CHANGELOG.txt>`_, and `MD5SUMS <http://mindboggle.info/data/MD5SUMS>`_,
| which describe the labeled nifti volumes (nii), vtk surfaces (vtk), and FreeSurfer files (mgh, etc.).
| Except where noted, all data are licensed under a Creative Commons License: |CC_license|_

|
|
.. image:: http://media.mindboggle.info/images/data/DKT_labels_width400px.png
|

------------------------------------------------------------------------------
Mindboggle-101 atlases
------------------------------------------------------------------------------
| Each of the 101 manually labeled individual brain surfaces and volumes is an atlas,
| a labeled or annotated brain image that can be used in registration-based labeling.
| We have combined the individual labels to create aggregate atlases here as well.

  **Volume atlases**

  - OASIS-TRT-20 joint fusion atlas `in OASIS-30 <http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_OASIS-30.nii.gz>`_ and `in MNI152 <http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152.nii.gz>`_ space [`2mm <http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_2mm.nii.gz>`_ version] (2013)
  - Corresponding label probabilities `in OASIS-30 <http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_OASIS-30.nii.gz>`_ and `in MNI152 <http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities_in_MNI152.nii.gz>`_ space

      Probabilistic labels of the 20 OASIS-TRT brains using joint fusion (Hongzhi Wang, 2013; distributed with ANTs),
      including a single volume of probabilities corresponding to the winning labels.
      Joint fusion was performed on the 20 brains after antsRegistration warped them
      to the brain and cerebellum `OASIS-30 Atropos template`_. Results are also given after
      `affine <http://mindboggle.info/data/templates/atropos/OASIS-30_Atropos_template_to_MNI152_affine.txt>`_
      transformation to the `OASIS-30 Atropos template in MNI152`_ space [`2mm <http://mindboggle.info/data/templates/atropos/OASIS-30_Atropos_template_in_MNI152_2mm.nii.gz>`_ version]

  - Software for `preparing labels <http://mindboggle.info/data/mindboggle101_extras/prep_OASIS-TRT-20_DKT31_CMA_labels.txt>`_, `running jointfusion <http://mindboggle.info/data/atlases/jointfusion/make_jointfusion_atlas.txt>`_, and `resampling to 2mm <http://mindboggle.info/data/resample2mm.txt>`_

  **Cortical surface atlases**

  - `DKT40 classifier atlas`_: FreeSurfer atlas (.gcs) from 40 of the Mindboggle-101 participants (2012)

------------------------------------------------------------------------------
Mindboggle-101 templates
------------------------------------------------------------------------------
|  A template is an unlabeled image used as a reference or standard, often for
|  registering other images to each other. Each one of the image volumes and
|  surfaces below was constructed by combining the images from multiple subjects.
|  ANTS brain templates were made with buildtemplateparallel.sh_,
|  head templates with antsMultivariateTemplateConstruction.sh_,
|  and FreeSurfer templates with make_freesurfer_template.py_.

  **Brain volumes**: ANTS nonlinear average templates (nii)

  - `OASIS-TRT-20 brain`_ template from 20 brains [`2mm <http://mindboggle.info/data/templates/ants/OASIS-TRT-20_template_in_MNI152_2mm.nii.gz>`_ version] (2013)
  - `OASIS-TRT-20 MNI152`_: OASIS-TRT-20 template `affine <http://mindboggle.info/data/templates/ants/OASIS-TRT-20_template_in_MNI152_affine.txt>`_ transformed to `MNI152`_ (2013)
  - `MMRR-21 brain`_ template from 21 brains (2012)
  - `MMRR-21 MNI152`_: MMRR-21 template `affine <http://mindboggle.info/data/templates/ants/MMRR-21_template_in_MNI152_affine.txt>`_ transformed to `MNI152`_ (2012)

  **Head volumes**: ANTS nonlinear average templates (nii)

  - `OASIS-TRT-20 head`_ template from 20 heads (2013) (`FSL's bet2 brain <http://mindboggle.info/data/templates/ants/OASIS-21_head_template_bet.nii.gz>`_)
  - `NKI-RS-22 head`_ template from 22 heads (2013)
  - `NKI-TRT-20 head`_ template from 20 heads (2013)
  - `MMRR-21 head`_ template from 21 heads (2013)
  - `HLN-12 head`_ template from 12 heads (2013)

  **Cortical surfaces**: FreeSurfer nonlinear average templates (tif)

  - `OASIS-TRT-20 surface`_ template from 20 brains (2012)
  - `NKI-RS-22 surface`_ template from 22 brains (2012)
  - `NKI-TRT-20 surface`_ template from 20 brains (2012)
  - `MMRR-21 surface`_ template from 21 brains (2012)
  - `HLN-12 surface`_ template from 12 brains (2012)

------------------------------------------------------------------------------
Mindboggle-101: Individual, manually labeled brain surfaces and volumes
------------------------------------------------------------------------------
|  The Mindboggle-101 dataset includes labeled anatomical regions for `101 healthy subjects`_.
|  The manually edited cortical labels follow sulcus landmarks according to the
|  Desikan-Killiany-Tourville (DKT) protocol (reference at top).

  - **OASIS-TRT-20** cortical labels [`nii <http://mindboggle.info/data/mindboggle101/OASIS-TRT-20_volumes.tar.gz>`_, `nii (MNI152) <http://mindboggle.info/data/mindboggle101/OASIS-TRT-20_volumes_in_MNI152.tar.gz>`_, `vtk <http://mindboggle.info/data/mindboggle101/OASIS-TRT-20_surfaces.tar.gz>`_, `mgh <http://mindboggle.info/data/mindboggle101/OASIS-TRT-20_freesurfer.tar.gz>`_]:
      All 20 subjects from the OASIS Test-Retest sample |OASIS-TRT www|_
  - **OASIS-TRT-20 whole-brain** labels [`nii <http://mindboggle.info/data/mindboggle101_extras/OASIS-TRT-20_DKT31_CMA_labels.tar.gz>`_, `nii (MNI152) <http://mindboggle.info/data/mindboggle101_extras/OASIS-TRT-20_DKT31_CMA_labels_in_MNI152.tar.gz>`_] by Neuromorphometrics_ |CC_license_nond|_
  - **NKI-RS-22** cortical labels [`nii <http://mindboggle.info/data/mindboggle101/NKI-RS-22_volumes.tar.gz>`_, `nii (MNI152) <http://mindboggle.info/data/mindboggle101/NKI-RS-22_volumes_in_MNI152.tar.gz>`_, `vtk <http://mindboggle.info/data/mindboggle101/NKI-RS-22_surfaces.tar.gz>`_,  `mgh <http://mindboggle.info/data/mindboggle101/NKI-RS-22_freesurfer.tar.gz>`_]:
      22 subjects from the Nathan Klein Institute / Rockland Sample |NKI-RS www|_
  - **NKI-TRT-20** cortical labels [`nii <http://mindboggle.info/data/mindboggle101/NKI-TRT-20_volumes.tar.gz>`_, `nii (MNI152) <http://mindboggle.info/data/mindboggle101/NKI-TRT-20_volumes_in_MNI152.tar.gz>`_, `vtk <http://mindboggle.info/data/mindboggle101/NKI-TRT-20_surfaces.tar.gz>`_, `mgh <http://mindboggle.info/data/mindboggle101/NKI-TRT-20_freesurfer.tar.gz>`_]:
      20 subjects from the Nathan Klein Institute / Test-Retest Sample |NKI-TRT www|_
  - **MMRR-21** cortical labels [`nii <http://mindboggle.info/data/mindboggle101/MMRR-21_volumes.tar.gz>`_, `nii (MNI152) <http://mindboggle.info/data/mindboggle101/MMRR-21_volumes_in_MNI152.tar.gz>`_, `vtk <http://mindboggle.info/data/mindboggle101/MMRR-21_surfaces.tar.gz>`_, `mgh <http://mindboggle.info/data/mindboggle101/MMRR-21_freesurfer.tar.gz>`_]:
      All 21 subjects in the Multi-Modal MRI Reproducibility Resource |MMRR www|_
  - **Extra-18** cortical labels [`nii <http://mindboggle.info/data/mindboggle101/Extra-18_volumes.tar.gz>`_, `nii (MNI152) <http://mindboggle.info/data/mindboggle101/Extra-18_volumes_in_MNI152.tar.gz>`_, `vtk <http://mindboggle.info/data/mindboggle101/Extra-18_surfaces.tar.gz>`_, `mgh <http://mindboggle.info/data/mindboggle101/Extra-18_freesurfer.tar.gz>`_]:
      - **HLN-12**:  All 12 subjects from the Human Language Network study
      - **MMRR-3T7T-2**:  2 subjects acquired like MMRR-21 (multimodal + 7T scans: |MMRR www|_)
      - **Colin27-1**:  Colin Holmes template (average of 27 scans)
      - **Twins-2**:  2 identical twins, including AK
      - **Afterthought-1**:  1 brain imager, SG
  - **fsaverage** [`nii and mgh <http://mindboggle.info/data/atlases/fsaverage.tar.gz>`_]:
      The figures above show the DKT cortical labeling protocol_ with `sulcus landmarks`_
      on FreeSurfer's fsaverage surface.


.. _`101 labeled brain images and a consistent human cortical labeling protocol`: http://www.frontiersin.org/Brain_Imaging_Methods/10.3389/fnins.2012.00171/full
.. _`OASIS-30 Atropos template`: http://mindboggle.info/data/templates/atropos/OASIS-30_Atropos_template.nii.gz
.. _`OASIS-30 Atropos template in MNI152`: http://mindboggle.info/data/templates/atropos/OASIS-30_Atropos_template_in_MNI152.nii.gz
.. _`OASIS-TRT-20 joint fusion atlas`: http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_labels.nii.gz
.. _`OASIS-TRT-20 joint fusion atlas in MNI152`: http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152.nii.gz
.. _`label probabilities in OASIS-30`: http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities.nii.gz
.. _`label probabilities in MNI152`: http://mindboggle.info/data/atlases/jointfusion/OASIS-TRT-20_jointfusion_DKT31_CMA_label_probabilities.nii.gz
.. _`DKT40 classifier atlas`: http://mindboggle.info/data/atlases/classifiers/DKTatlas40.tar.gz


.. _MD5SUMS: http://mindboggle.info/data/MD5SUMS
.. _CC_license: http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US
.. |CC_license| image:: http://i.creativecommons.org/l/by-nc-sa/3.0/80x15.png
.. _`101 healthy subjects`: http://media.mindboggle.info/images/data/Mindboggle101_table.pdf
.. _labels: http://mindboggle.info/data/mindboggle101/protocol.txt
.. _protocol: http://mindboggle.info/data/mindboggle101/protocol.txt
.. _`sulcus landmarks`: http://media.mindboggle.info/images/data/DKT_sulci_table.pdf
.. _`MMRR www`: http://www.nitrc.org/projects/multimodal
.. _`NKI-RS www`: http://fcon_1000.projects.nitrc.org/indi/pro/nki.html
.. _`NKI-TRT www`: http://fcon_1000.projects.nitrc.org/indi/pro/eNKI_RS_TRT/FrontPage.html
.. _`OASIS-TRT www`: http://www.oasis-brains.org/app/action/BundleAction/bundle/OAS1_RELIABILITY
.. |MMRR www| image:: images/link-brown-12x12.png
.. |NKI-RS www| image:: images/link-brown-12x12.png
.. |NKI-TRT www| image:: images/link-brown-12x12.png
.. |OASIS-TRT www| image:: images/link-brown-12x12.png
.. _Neuromorphometrics: http://neuromorphometrics.com
.. _CC_license_nond: http://creativecommons.org/licenses/by-nc-nd/3.0/deed.en_US
.. |CC_license_nond| image:: http://i.creativecommons.org/l/by-nc-nd/3.0/80x15.png


.. _numbers: http://media.mindboggle.info/images/data/DKT_label_table.pdf
.. _buildtemplateparallel.sh: data/templates/ants/buildtemplateparallel.sh
.. _antsMultivariateTemplateConstruction.sh: data/templates/ants/antsMultivariateTemplateConstruction.sh
.. _make_freesurfer_template.py: data/templates/freesurfer/make_freesurfer_template.txt
.. _`MMRR-21 brain`: http://mindboggle.info/data/templates/ants/MMRR-21_template.nii.gz
.. _`MMRR-21 MNI152`: http://mindboggle.info/data/templates/ants/MMRR-21_template_in_MNI152.nii.gz
.. _`MNI152`: http://mindboggle.info/data/templates/MNI152_T1_1mm_brain.nii.gz
.. _`OASIS-TRT-20 brain`: http://mindboggle.info/data/templates/ants/OASIS-TRT-20_template.nii.gz
.. _`OASIS-TRT-20 MNI152`: http://mindboggle.info/data/templates/ants/OASIS-TRT-20_template_in_MNI152.nii.gz
.. _`affine`: http://mindboggle.info/data/templates/ants/OASIS-TRT-20_template_in_MNI152_affine.txt
.. _`HLN-12 head`: http://mindboggle.info/data/templates/ants/HLN-12_head_template.nii.gz
.. _`MMRR-21 head`: http://mindboggle.info/data/templates/ants/MMRR-21_head_template.nii.gz
.. _`NKI-RS-22 head`: http://mindboggle.info/data/templates/ants/NKI-RS-22_head_template.nii.gz
.. _`NKI-TRT-20 head`: http://mindboggle.info/data/templates/ants/NKI-TRT-20_head_template.nii.gz
.. _`OASIS-TRT-20 head`: http://mindboggle.info/data/templates/ants/OASIS-TRT-20_head_template.nii.gz
.. _`HLN-12 surface`: http://mindboggle.info/data/templates/freesurfer/HLN-12_surface_template.tar.gz
.. _`MMRR-21 surface`: http://mindboggle.info/data/templates/freesurfer/MMRR-21_surface_template.tar.gz
.. _`NKI-RS-22 surface`: http://mindboggle.info/data/templates/freesurfer/NKI-RS-22_surface_template.tar.gz
.. _`NKI-TRT-20 surface`: http://mindboggle.info/data/templates/freesurfer/NKI-TRT-20_surface_template.tar.gz
.. _`OASIS-TRT-20 surface`: http://mindboggle.info/data/templates/freesurfer/OASIS-TRT-20_surface_template.tar.gz


------------------------------------------------------------------------------
Other templates and manually labeled brains
------------------------------------------------------------------------------
| The following images are not from the Mindboggle-101 data above, and the manual labels are not the same
| as those of the DKT labeling protocol used for the Mindboggle-101 data above:

  **Tissue-segmented templates**: created by Nicholas Tustison for use with antsAtroposN4.sh

  - `OASIS-30 <http://mindboggle.info/data/templates/atropos/OASIS-30_Atropos_template.tar.gz>`_: from 30 `MICCAI challenge <http://mindboggle.info/data/templates/atropos/MICCAI_2012_Workshop_v2.pdf>`_ OASIS images (2013)
  - `NKI-30 <http://mindboggle.info/data/templates/atropos/NKI-30_Atropos_template.tar.gz>`_: from 30 NKI images (2013)
  - `MMRR-41 <http://mindboggle.info/data/templates/atropos/MMRR-41_Atropos_template.tar.gz>`_: from 41 MMRR images (2013)
  - `IXI <http://mindboggle.info/data/templates/atropos/IXI_Atropos_template.tar.gz>`_: from IXI images (2013)

  **Templates**: built with buildtemplateparallel.sh_ with images from `2009 evaluation`_ (nii)

  - `CUMC12 brain`_ template from 12 brains (2010)
  - `LPBA40 brain`_ template from 40 brains (2011)
  - See Satrajit Ghosh's `pediatric template`_ of 31 brains (2011)

  **Atlases**: manually labeled volumes from `2009 evaluation`_ (nii)

  - CUMC-12_: 12 labeled brains (2009)
  - IBSR-18_: 18 labeled brains (2009)
  - MGH-10_: 10 labeled brains (2009)

|

.. image:: http://media.mindboggle.info/images/data/evaluation2009_80atlases.png

.. _`2009 evaluation`: http://www.mindboggle.info/papers/evaluation_NeuroImage2009.php
.. _Atropos: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3297199/
.. _CUMC-12: http://mindboggle.info/papers/evaluation_NeuroImage2009/data/CUMC12.tar.gz
.. _IBSR-18: http://mindboggle.info/papers/evaluation_NeuroImage2009/data/IBSR18.tar.gz
.. _MGH-10: http://mindboggle.info/papers/evaluation_NeuroImage2009/data/MGH10.tar.gz
.. _`CUMC12 brain`: http://mindboggle.info/data/templates/ants/CUMC-12_template.nii.gz
.. _`LPBA40 brain`: http://mindboggle.info/data/templates/ants/LPBA-40_template.nii.gz
.. _`pediatric template`: http://www.mit.edu/~satra/research/pubdata/index.html

.. include:: ./links.txt
