#!/usr/bin/env python
"""
Utility functions.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#-----------------------------------------------------------------------------
# Functions for running shell commands:
#-----------------------------------------------------------------------------
def execute(cmd, type='subprocess'):
    """
    Execute command by either subprocess.call or os.system.

    Parameters
    ----------
    cmd : sequence or string
        command with arguments
    type : string
        how to execute {os, subprocess}

    Examples
    --------
    >>> from mindboggle.utils.utils import execute
    >>> cmd = ['ls', '-l', '-a', '.']
    >>> type = 'subprocess'
    >>> execute(cmd, type)
    >>> type = 'os'
    >>> execute(cmd, type)
    >>> cmd = 'ls -l -a .'
    >>> execute(cmd)

    """
    from subprocess import call
    import sys

    if isinstance(cmd, str):
        print(cmd)
    else:
        print(' '.join(cmd))

    # Use subprocess.call:
    if type == 'subprocess':
        try:
            retcode = call(cmd)
            if retcode < 0:
                print >>sys.stderr, "Child terminated by signal", -retcode
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

    # Use os.system:
    elif type == 'os':
        from os import system

        if isinstance(cmd, str):
            pass
        else:
            cmd = ' '.join(cmd)
        try:
            system(cmd)
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

    else:
        sys.exit('Select either "subprocess" or "os" for execution type.')


