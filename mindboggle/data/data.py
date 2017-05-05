#!/usr/bin/env python
"""
Functions for fetching information related to Mindboggle data.

Authors:
    - Arno Klein, 2017  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2017,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def fetch_file_path(file_name):
    """
    Return the path to files in the data directory in Mindboggle.

    Parameters
    ----------
    file_name : str
        name of the file whose full path is to be returned

    Returns
    -------
    file_path : str
        full path to file

    """
    import os

    file_path = os.path.join(os.path.dirname(__file__), file_name)

    return file_path
