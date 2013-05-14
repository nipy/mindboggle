# Functions for running shell commands:

from subprocess import call
from os import system
from sys import stderr

def subprocess_call(cmd,args):
    print(cmd + ' ' + args)
    try:
        retcode = call(cmd + ' ' + args, shell=True)
        if retcode < 0:
            print >>stderr, "Child terminated by signal", -retcode
    except OSError, e:
        print >>stderr, "Execution failed:", e

def os_system(cmd):
    print(cmd)
    try:
        system(cmd)
    except OSError, e:
        print >>stderr, "Execution failed:", e
