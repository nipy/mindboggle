#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Mindboggle : package for automated human brain image labeling and morphometry.

export PYTHONPATH=$PYTHONPATH:/desk/temp/lib/python2.7/site-packages/
python setup.py install --prefix=/desk/temp
python setup.py develop

"""

from distribute_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

setup(
    name = "Mindboggle",
    version = "0.01",
    packages = find_packages(),
    #scripts = ['mindboggle/pipeline.sh'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires = ['docutils>=0.3'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
    },

    # metadata for upload to PyPI
    author = "Arno Klein",
    author_email = "arno@binarybottle.com",
    description = "Automated human brain image automated anatomical labeling and shape analysis",
    long_description = "Mindboggle is a package for automated anatomical labeling and morphometry of human brain images",
    license = "Apache 2.0",
    keywords = "Mindboggle human brain MRI automated labeling parcellation morphometry",
    url = "http://mindboggle.info/",
    download_url = "http://mindboggle.info/",
)


"""
import sys
from glob import glob

# Import build helpers
try:
    from nisext.sexts import package_check, get_comrec_build
except ImportError:
    raise RuntimeError('Need nisext package from nibabel installation'
                       ' - please install nibabel first')

from build_docs import cmdclass, INFO_VARS

# Add custom commit-recording build command
cmdclass['build_py'] = get_comrec_build('mindboggle')

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    # The quiet=True option will silence all of the name setting warnings
    config.get_version('mindboggle/__init__.py') # sets config.version
    config.add_subpackage('mindboggle', 'mindboggle')
    return config

################################################################################
# For some commands, use setuptools

if len(set(('develop', 'bdist_egg', 'bdist_rpm', 'bdist', 'bdist_dumb',
            'bdist_wininst', 'install_egg_info', 'egg_info', 'easy_install',
            )).intersection(sys.argv)) > 0:
    from setup_egg import extra_setuptools_args

# extra_setuptools_args can be defined from the line above, but it can
# also be defined here because setup.py has been exec'd from
# setup_egg.py.
if not 'extra_setuptools_args' in globals():
    extra_setuptools_args = dict()

# Hard and soft dependency checking
package_check('nibabel', INFO_VARS['NIBABEL_MIN_VERSION'])
package_check('numpy', INFO_VARS['NUMPY_MIN_VERSION'])
package_check('nipype', INFO_VARS['NIPYPE_MIN_VERSION'])

################################################################################
# Import the documentation building classes.

try:
    from build_docs import cmdclass
except ImportError:
    "" Pass by the doc build gracefully if sphinx is not installed ""
    print "Sphinx is not installed, docs cannot be built"
    cmdclass = {}


################################################################################

def main(**extra_args):
    from numpy.distutils.core import setup

    setup(name=INFO_VARS['NAME'],
          maintainer=INFO_VARS['MAINTAINER'],
          maintainer_email=INFO_VARS['MAINTAINER_EMAIL'],
          description=INFO_VARS['DESCRIPTION'],
          long_description=INFO_VARS['LONG_DESCRIPTION'],
          url=INFO_VARS['URL'],
          download_url=INFO_VARS['DOWNLOAD_URL'],
          license=INFO_VARS['LICENSE'],
          classifiers=INFO_VARS['CLASSIFIERS'],
          author=INFO_VARS['AUTHOR'],
          author_email=INFO_VARS['AUTHOR_EMAIL'],
          platforms=INFO_VARS['PLATFORMS'],
          version=INFO_VARS['VERSION'],
          requires=INFO_VARS['REQUIRES'],
          configuration = configuration,
          cmdclass = cmdclass,
          scripts = glob('bin/*'),
          **extra_args)



if __name__ == "__main__":
    main(**extra_setuptools_args)
"""
