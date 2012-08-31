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

    http://www.mindboggle.info

Code
====

You can find our sources and single-click downloads::

* `Main repository`_ on Github
* Documentation_ for all releases and current development tree
* Download as a tar/zip file the `current trunk`_
* Downloads of all `available releases`_

.. _main repository: http://github.com/binarybottle/mindboggle
.. _Documentation: http://www.mindboggle.info/documentation
.. _current trunk: http://github.com/binarybottle/mindboggle/archives/master
.. _available releases: http://github.com/binarybottle/mindboggle/downloads

License
=======

The mindboggle package, including all examples, code snippets and attached
documentation, is licensed under the terms of the Apache 2.0 license.
Please see the LICENSE file in the mindboggle distribution.

File tree
=========

mindboggle/
    README             # this file
    INSTALL            # installation instructions
    LICENSE            # copyright and license information
    THANKS             # acknowledgments
    setup.py           # installation program
    data/              # atlas and template data
    doc/               # documentation (build with Sphinx)
    mindboggle_cpp/    # C++ code that needs to be compiled
    mindboggle/        # Python software directory
        pipeline.py    # main program (nipype pipeline)
        extract/       # software for extracting features
        label/         # software for labeling image data
        measure/       # software for measuring shapes
        utils/         # utilities
        <misc. files for distributing software>
