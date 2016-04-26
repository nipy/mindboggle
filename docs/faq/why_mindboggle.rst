.. _why_mindboggle:

------------------------------------------------------------------------------
 Why Mindboggle?
------------------------------------------------------------------------------

Purpose
.......

The purpose of Mindboggle is to improve the accuracy, precision, and
consistency of labeling and morphometry of brain imaging data,
and to promote open science by making all data, software, and documentation
freely and openly available.

Motivation
..........

There is a strong tendency for neuroscientists who are not software developers
to assume that questions of brain image "preprocessing", including registration
and labeling, have been satisfactorily resolved or are not even worthy of mention
[carp2012]_.  We believe that conventional methods do not establish
reasonable correspondences across brain images and threaten to undermine the
credibility of the science and the scientists who use them in the service of
answering the neuroscientific and clinical questions that interest them.

.. [carp2012]
   Carp, J. (2012) *The secret lives of experiments:
   Methods reporting in the fMRI literature*. NeuroImage, 63(1), 289–300.
   doi:10.1016/j.neuroimage.2012.07.004

History
.......

Back in 1998, `Arno Klein <http://binarybottle.com>`_
left Caltech to visit Cornell Medical School for a summer,
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
People who contributed to the codebase during this period include:

    - `Forrest Bao <https://sites.google.com/site/forrestbao/>`_
    - `Satrajit Ghosh <http://mit.edu/~satra>`_
    - Joachim Giard
    - Yrjö Häme
    - Eliezer Stavsky

In 2012, the NIMH funded a 1-year U01 supplement to expand Mindboggle's
capabilities, and to apply Mindboggle's structural shape analysis to multiple
imaging modalities such as fMRI, dMRI, and genetic data.
Tank Think Labs joined this U01 effort:

    - `Satrajit Ghosh <http://mit.edu/~satra>`_
    - Nolan Nichols
    - Brian Rossa
    - Oliver Hinds

In 2013, Arno, Rich Stoner, and Jason Bohland teamed up at the Human Brain
Mapping conference to win a hackathon challenge, having coded up an online
interactive visualization of Allen Brain Institute human brain data
in two days. Pointing to a cortical region displaying Mindboggle shape
information revealed corresponding genetic expression data.

In 2015, Mindboggle was used to process ADNI and AddNeuroMed data for
an international `Alzheimer's disease challenge <https://www.synapse.org/#!Synapse:syn2290704/wiki/60828>`_.
Teams performed statistical analyses on Mindboggle shape measures to
try and determine which brains had Alzheimer's, mild cognitive impairment, or
were healthy, and to try and estimate a cognitive measure
(mini-mental state exam score).

In 2016, Arno prepared Mindboggle for broader public use:

    - moved its GitHub repository to the nipy.org community's GitHub account
    - ported Mindboggle from Python 2 to Python 3
    - moved its documentation to readthedocs.org (generated every time a commit is made to the GitHub repository)
    - wrote docstring tests for almost every function
    - set up continuous testing on circleci.com (tests are run every time a commit is made to the GitHub repository)


Please see our updated Mindboggle `team page <http://mindboggle.info/people.html>`_.
