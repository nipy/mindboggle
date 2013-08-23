.. _mindboggle-faq:

=======
History
=======

.. _purpose:

Purpose
-------

The purpose of Mindboggle is to improve the accuracy, precision, and
consistency of labeling and morphometry of brain imaging data,
and to promote open science by making all data, software, and documentation
freely and openly available.


.. _motivation:

Motivation
----------

There is a strong tendency for neuroscientists who are not software developers
to assume that questions of brain image "preprocessing", including registration
and labeling, have been satisfactorily resolved or are not even worthy of mention
([carp2012]_).  We believe that conventional methods do not establish
reasonable correspondences across brain images and threaten to undermine the
credibility of the science and the scientists who use them in the service of
answering the neuroscientific and clinical questions that interest them.

.. [carp2012]
   Carp, J. (2012) *The secret lives of experiments:
   Methods reporting in the fMRI literature*. NeuroImage, 63(1), 289–300.
   doi:10.1016/j.neuroimage.2012.07.004


.. _history:

History
-------

Back in 1998, `Arno Klein`_ left Caltech to visit Cornell Medical School for a summer,
to learn a bit about human brain imaging.  Within a day of arriving in New York,
he witnessed two medical students arguing over what regions of a patient's brain
were active according to fMRI BOLD data printed on some sheets of paper.
To do this, they compared locations of BOLD signal in the young man's
brain with what appeared to them to be corresponding topographical locations
in the 2-D slices of a dead, 60-year French woman's half cerebrum from 20 years prior,
the venerable Talairach atlas.  That clinched it.
Arno stayed on to create the initial Mindboggle software package in Matlab
as part of his doctoral dissertation.

In 2007, Arno was invited to join the faculty at Columbia University
to continue its development, whereupon he started writing proposals.

In 2009, with generous funding from the National Institute of Mental Health
(3-year NIMH R01 #MH084029), he assembled a team to tackle the problem anew.
People who have contributed to the codebase during this period include:

    - Forrest Sheng Bao
    - Satrajit Ghosh
    - Joachim Giard
    - Yrjö Häme
    - Eliezer Stavsky

In 2012, the NIMH funded a 1-year U01 supplement to expand this work to
integrate multiple imaging modalities (fMRI and dMRI) with the structural MRI
shape analysis.  Joining this effort will be:

    - Nolan Nichols

For more information about the Mindboggle team, see http://mindboggle.info/people/.

.. include:: ../links.txt
