#!/usr/bin/env python
"""
Functions for input/output via URIs.

Authors:
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


#-----------------------------------------------------------------------------
# Get data through a URL call:
#-----------------------------------------------------------------------------
def retrieve_data(data_file, url, hashes={}, cache_env='', cache='',
                  return_missing=False):
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
    url : string
        URL for data file
    hashes : dictionary
        file names and md5 hashes (if empty, simply download file from url)
    cache_env : string
        environment variable name for cache path
    cache : string
        in case cache_env is not set, use as cache directory
    return_missing : Boolean
        if data_file not in hash, simply return data_file

    Returns
    -------
    data_path : string
        data file name (full path)

    Examples
    --------
    >>> import os
    >>> from mindboggle.data import hashes_url
    >>> from mindboggle.utils.io_uri import retrieve_data
    >>> hashes, url, cache_env, cache = hashes_url()
    >>> data_file = hashes.keys()[0]
    >>> retrieve_data(data_file, url, hashes, cache_env, cache)

    """
    import os
    import sys
    import shutil

    from mindboggle.utils.io_uri import get_data, get_hash

    #-------------------------------------------------------------------------
    # If hashes provided, go through steps to check/download file:
    #-------------------------------------------------------------------------
    if hashes:

        if not cache_env:
            cache_env = 'MINDBOGGLE_CACHE'
        if not cache:
            cache = os.path.join(os.environ['HOME'], 'hash_temp')

        #---------------------------------------------------------------------
        # Check hash table for file:
        #---------------------------------------------------------------------
        if data_file not in hashes.keys():
            if return_missing:
                data_path = data_file
                print("Retrieved file not in hashes: {0}".format(data_path))
                return data_path
            else:
                sys.exit("Data file '{0}' not in hash table.".
                format(data_file))
        else:
            stored_hash = hashes[data_file]

            #-----------------------------------------------------------------
            # Create missing cache and hash directories:
            #-----------------------------------------------------------------
            if cache_env in os.environ.keys():
                cache = os.environ[cache_env]
            if not os.path.exists(cache):
                print("Create missing cache directory: {0}".format(cache))
                os.mkdir(cache)
            hash_dir = os.path.join(cache, stored_hash)
            if not os.path.exists(hash_dir):
                print("Create missing hash directory: {0}".format(hash_dir))
                os.mkdir(os.path.join(hash_dir))
    
            #-----------------------------------------------------------------
            # Check hash subdirectory for file:
            #-----------------------------------------------------------------
            data_path = os.path.join(hash_dir, data_file)
            if os.path.exists(data_path):
                return data_path
    
            #-----------------------------------------------------------------
            # If file not in cache, download file, compute hash, and verify:
            #-----------------------------------------------------------------
            else:
                print("Retrieve file: {0}".format(url+data_file))
    
                # Download file as a temporary file:
                temp_file = get_data(url+data_file)
    
                # Compute the file's hash:
                data_hash = get_hash(temp_file)
    
                # If hash matches name of the hash directory, save file:
                if os.path.join(cache, data_hash) == hash_dir:
                    print("Copy file to cache: {0}".format(data_path))
                    shutil.copyfile(temp_file, data_path)
                    return data_path
                else:
                    print("Retrieved file's hash does not match stored hash.")

    #-------------------------------------------------------------------------
    # If hashes not provided, simply download file:
    #-------------------------------------------------------------------------
    else:
        # Download file as a temporary file:
        data_path = get_data(url+data_file)
        print("Retrieved file: {0}".format(data_path))
        return data_path


def get_data(url, output_file=''):
    """
    Get data from URL.

    Parameters
    ----------
    url : string
        URL for data file
    output_file : string
        name of output file (full path)

    Returns
    -------
    output_file : string
        name of output file (full path)

    Examples
    --------
    >>> from mindboggle.utils.io_uri import get_data
    >>> from mindboggle.data import hashes_url
    >>> hashes, url = hashes_url()
    >>> data_file = hashes.keys()[0]
    >>> output_file = 'test_output.nii.gz'
    >>> get_data(url, output_file)

    """
    import urllib

    # Download file as a temporary file:
    if output_file:
        output_file, foo = urllib.urlretrieve(url, output_file)
    # Download file to specified output:
    else:
        output_file, foo = urllib.urlretrieve(url)

    return output_file


def get_hash(data_file):
    """
    Get hash of data file.

    Parameters
    ----------
    data_file : string
        data file name

    Returns
    -------
    hash : string
        hash of data file

    Examples
    --------
    >>> import os
    >>> from mindboggle.data import hashes_url
    >>> from mindboggle.utils.io_uri import get_hash
    >>> cache_env = 'MINDBOGGLE_CACHE'
    >>> hashes, url = hashes_url()
    >>> data_file = hashes.keys()[0]
    >>> data_path = os.path.join(os.environ[cache_env], data_file)
    >>> get_hash(data_path)

    """
    import hashlib

    # Compute the file's hash:
    hash = hashlib.md5(open(data_file, 'rb').read()).hexdigest()

    return hash