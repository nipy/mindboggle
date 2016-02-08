#!/usr/bin/env python
"""
Utility functions.

Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def execute(cmd, type='os'):
    """
    Execute command by either subprocess.call or os.system.

    Parameters
    ----------
    cmd : sequence (string also permitted if type=='os')
        command with arguments
    type : string
        how to execute {os, subprocess}

    Examples
    --------
    >>> from mindboggle.guts.utilities import execute # doctest: +SKIP
    >>> cmd = ['date', '-r', '0'] # doctest: +SKIP
    >>> type = 'subprocess' # doctest: +SKIP
    >>> execute(cmd, type) # doctest: +SKIP
    Wed Dec 31 19:00:00 EST 1969
    >>> type = 'os' # doctest: +SKIP
    >>> execute(cmd, type) # doctest: +SKIP
    Wed Dec 31 19:00:00 EST 1969
    >>> cmd = 'date -r 0' # doctest: +SKIP
    >>> execute(cmd) # doctest: +SKIP
    Wed Dec 31 19:00:00 EST 1969

    """
    from subprocess import call

    verbose = False
    if verbose:
        if isinstance(cmd, str):
            print(cmd)
        else:
            print(' '.join(cmd))

    # Use subprocess.call:
    if type == 'subprocess':
        try:
            retcode = call(cmd)
            if retcode < 0:
                raise IOError("Child terminated by signal: retcode {0}".
                              format(retcode))
        except OSError, e:
            raise OSError("Execution failed: {0}".format(e))

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
            raise OSError("Execution failed: {0}".format(e))

    else:
        raise IOError('Select either "subprocess" or "os" for execution type.')


def list_strings(string1='', string2='', string3='', string4=''):
    """
    Put strings in a list.

    Parameters
    ----------
    string1 : string
    string2 : string
    string3 : string
    string4 : string

    Returns
    -------
    string_list : list of strings

    Examples
    --------
    >>> from mindboggle.guts.utilities import list_strings
    >>> string1 = 'a b c'
    >>> string2 = 'd e f'
    >>> string3 = ''
    >>> string4 = 'j k l'
    >>> output_file = ''
    >>> string_list = list_strings(string1, string2, string3, string4)
    >>> string_list
    ['a b c', 'd e f', 'j k l']

    """

    string_list = []
    if string1 and isinstance(string1, str):
        string_list.append(string1)
    if string2 and isinstance(string1, str):
        string_list.append(string2)
    if string3 and isinstance(string1, str):
        string_list.append(string3)
    if string4 and isinstance(string1, str):
        string_list.append(string4)

    return string_list


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()