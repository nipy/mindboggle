import os

#from .info import (LONG_DESCRIPTION as __doc__,
#                   __version__)
#__doc__ += """
#"""

# Set up package information function
from .pkg_info import get_pkg_info as _get_pkg_info
get_info = lambda : _get_pkg_info(os.path.dirname(__file__))

# module imports
#from . import blah as blah
# object imports
#from .blah import blah, blah
