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
def subprocess_call(cmd,args):
    from subprocess import call
    from sys import stderr
    print(cmd + ' ' + args)
    try:
        retcode = call(cmd + ' ' + args, shell=True)
        if retcode < 0:
            print >>stderr, "Child terminated by signal", -retcode
    except OSError, e:
        print >>stderr, "Execution failed:", e


def os_system(cmd):
    from os import system
    from sys import stderr
    print(cmd)
    try:
        system(cmd)
    except OSError, e:
        print >>stderr, "Execution failed:", e
