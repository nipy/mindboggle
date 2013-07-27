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
    >>> from mindboggle.utils.io_uri import retrieve_data
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
    import shutil

    from mindboggle.utils.io_uri import get_data, get_hash

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
            temp_file = get_data(uri+data_file)

            # Compute the file's hash:
            data_hash = get_hash(data_file)

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
        data_path = get_data(uri+data_file)
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
    >>> data_file = 'OASIS-TRT-20_DKT31_CMA_jointfusion_labels_in_MNI152.nii.gz'
    >>> uri = 'http://mindboggle.info/data/atlases/jointfusion/'
    >>> url = uri + data_file
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
    >>> from mindboggle.utils.io_uri import get_hash
    >>> data_file = 'OASIS-TRT-20_DKT31_CMA_jointfusion_labels_in_MNI152.nii.gz'
    >>> get_hash(data_file)

    """
    import hashlib

    # Compute the file's hash:
    hash = hashlib.md5(open(data_file, 'rb').read()).hexdigest()

    return hash