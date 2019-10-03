import os

#from .info import (LONG_DESCRIPTION as __doc__,
#                  __version__)
#__doc__ += """
#"""

# Set up package information function
try:
    from .pkg_info import get_pkg_info as _get_pkg_info
except:
    get_info = lambda: ""
else:
    get_info = lambda : _get_pkg_info(os.path.dirname(__file__))

# module imports
#from . import blah as blah
# object imports
#from .blah import blah, blah

INIT_MSG = "Running {packname} version {version} (latest: {latest})".format
latest = {"version": 'Unknown'}
try:
    from .version import __version__
    import etelemetry
    latest = etelemetry.get_project("nipy/mindboggle")
except Exception as e:
    print("Could not check for version updates: ", e)
finally:
    print(INIT_MSG(packname='mindboggle',
                   version=__version__,
                   latest=latest["version"]))

