#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Mindboggle: package for automated human brain image labeling and morphometry.

"""

import os
import sys
from os.path import join as pjoin

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

# For some commands, use setuptools.
if len(set(('develop', 'bdist_egg', 'bdist_rpm', 'bdist', 'bdist_dumb',
            'bdist_wininst', 'install_egg_info', 'egg_info', 'easy_install',
            )).intersection(sys.argv)) > 0:
    # setup_egg imports setuptools setup, thus monkeypatching distutils.
    import setup_egg

from distutils.core import setup

# Python 2 to 3 build
from nisext.py3builder import build_py
# Commit hash writing, and dependency checking
from nisext.sexts import get_comrec_build, package_check
cmdclass = {'build_py': get_comrec_build('mindboggle', build_py)}

# Get version and release info, which is all stored in mindboggle/info.py
ver_file = os.path.join('mindboggle', 'info.py')
exec(open(ver_file).read())

# Do dependency checking
package_check('numpy', NUMPY_MIN_VERSION)

extra_setuptools_args = {}
if 'setuptools' in sys.modules:
    extra_setuptools_args = dict(
        tests_require=['nose'],
        test_suite='nose.collector',
        zip_safe=False,
        extras_require = dict(
            doc='Sphinx>=0.3',
            test='nose>=0.10.1')
    )

def main(**extra_args):
    setup(name=NAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          url=URL,
          download_url=DOWNLOAD_URL,
          license=LICENSE,
          classifiers=CLASSIFIERS,
          author=AUTHOR,
          author_email=AUTHOR_EMAIL,
          platforms=PLATFORMS,
          version=VERSION,
          requires=REQUIRES,
          provides=PROVIDES,
          packages     = ['mindboggle',
                          'mindboggle.evaluate',
                          'mindboggle.features',
                          'mindboggle.labels',
                          'mindboggle.labels.protocol',
                          'mindboggle.shapes',
                          'mindboggle.tables',
                          'mindboggle.utils'],
          #                'mindboggle.testing',
          #                'mindboggle.tests',
          #                'mindboggle.benchmarks',
          # The package_data spec has no effect for me (on python 2.6) -- even
          # changing to data_files doesn't get this stuff included in the source
          # distribution -- not sure if it has something to do with the magic
          # above, but distutils is surely the worst piece of code in all of
          # python -- duplicating things into MANIFEST.in but this is admittedly
          # only a workaround to get things started -- not a solution

          package_data = {'mindboggle':
                          [pjoin('labels', 'protocol', '*.txt')
                           #pjoin('extract', 'medial_surfaces', '*'),
                           #pjoin('measure', 'surface_measures', '*.cpp')
                          ]},
          #scripts      = [pjoin('bin', 'parrec2nii'),
          #                pjoin('bin', 'nib-ls'),
          #                ],
          cmdclass = cmdclass,
          **extra_args
         )

if __name__ == "__main__":
    main(**extra_setuptools_args)


