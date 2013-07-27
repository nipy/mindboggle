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


#-----------------------------------------------------------------------------
# Get data through a URL call:
#-----------------------------------------------------------------------------
def retrieve_data(data_file, uri, hashes={}, cache_env='', cache_dir=''):
    """
    Get data file through a URL call and check its hash.

    Steps ::
        If hashes provided:
            1. Check hash table for data file.
            2. Check hash subdirectory within cache directory for data file.
            3. If data file not in cache, download file, compute hash,
               and verify hash.
            4. If hash correct, save file.
        Otherwise, simply download file.

    Parameters
    ----------
    data_file : string
        data file name
    uri : string
        URL for data file
    hashes : dictionary
        file names and md5 hashes
    cache_env : string
        environment variable name for cache path
    cache_dir : string
        in case cache_env is not set, use as cache directory

    Returns
    -------
    data_path : string
        data file name (full path)

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.utils import retrieve_data
    >>> data_file = 'OASIS-TRT-20_DKT31_CMA_jointfusion_labels_in_MNI152.nii.gz'
    >>> uri = 'http://mindboggle.info/data/atlases/jointfusion/'
    >>> hashes = {}
    >>> hashes[data_file] = '082f19c118b428e49fbb56c55756c676'
    >>> cache_env = 'MINDBOGGLE_CACHE'
    >>> cache_dir = os.path.join(os.environ['HOME'], 'mindboggle_temp')
    >>> retrieve_data(data_file, uri, hashes, cache_env, cache_dir)

    """
    import os
    import sys
    import urllib
    import hashlib
    import shutil

    #-------------------------------------------------------------------------
    # If hashes provided, go through steps to check/download file:
    #-------------------------------------------------------------------------
    if hashes:

        if not cache_env:
            cache_env = 'MINDBOGGLE_CACHE'
        if not cache_dir:
            cache_dir = os.path.join(os.environ['HOME'], 'hash_temp')

        #---------------------------------------------------------------------
        # Check hash table for file:
        #---------------------------------------------------------------------
        if data_file not in hashes.keys():
            sys.exit("Image file '{0}' not in hash table.".format(data_file))
        else:
            stored_hash = hashes[data_file]

        #---------------------------------------------------------------------
        # Create missing cache and hash directories:
        #---------------------------------------------------------------------
        if cache_env in os.environ.keys():
            cache_dir = os.environ[cache_env]
        if not os.path.exists(cache_dir):
            print("Create missing cache directory: {0}".format(cache_dir))
            os.mkdir(cache_dir)
        hash_dir = os.path.join(cache_dir, stored_hash)
        if not os.path.exists(hash_dir):
            print("Create missing hash directory: {0}".format(hash_dir))
            os.mkdir(os.path.join(hash_dir))

        #---------------------------------------------------------------------
        # Check hash subdirectory for file:
        #---------------------------------------------------------------------
        data_path = os.path.join(hash_dir, data_file)
        if os.path.exists(data_path):
            return data_path

        #---------------------------------------------------------------------
        # If file not in cache, download file, compute hash, and verify:
        #---------------------------------------------------------------------
        else:
            print("Retrieve file: {0}".format(uri+data_file))

            # Download file as a temporary file:
            temp_file, foo = urllib.urlretrieve(uri+data_file)

            # Compute the file's hash:
            data_hash = hashlib.md5(open(temp_file, 'rb').read()).hexdigest()

            # If hash matches name of the hash directory, save file:
            if os.path.join(cache_dir, data_hash) == hash_dir:
                print("Copy file to cache directory: {0}".format(data_path))
                shutil.copyfile(temp_file, data_path)
                return data_path
            else:
                print("Retrieved file's hash does not matched stored hash.")

    #-------------------------------------------------------------------------
    # If hashes not provided, simply download file:
    #-------------------------------------------------------------------------
    else:
        # Download file as a temporary file:
        data_path, foo = urllib.urlretrieve(uri+data_file)
        print("Retrieved file: {0}".format(data_path))
        return data_path
