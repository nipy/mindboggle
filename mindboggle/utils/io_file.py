#!/usr/bin/env python
"""
This Python library reads and writes different file types.

Authors:
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein, 2012  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import re

def read_columns(filename, n_columns=1, trail=False):
    """
    Read n-column text file.

    Parameters
    ----------
    filename :  name of text file [string]
    n_columns :  number of columns to extract [integer]
    trail :  combine all remaining columns as a string
            in the final list [Boolean]

    Returns
    -------
    columns :  a list of lists of strings, one list per column of text.

    """
    import re

    Fp = open(filename, 'r')
    lines = Fp.readlines()
    columns = [[] for x in range(n_columns)]
    for line in lines:
        if line:
            row = re.findall(r'\S+', line)
            if len(row) >= n_columns:
                for icolumn in range(n_columns):
                    if trail and icolumn == n_columns - 1:
                        columns[icolumn].append(' '.join(row[icolumn::]))
                    else:
                        columns[icolumn].append(row[icolumn])
            else:
                import os
                os.error('The number of columns in {0} is less than {1}.'.format(
                         filename, n_columns))
    Fp.close()

    return columns

def write_table(labels, columns, column_names, table_file):
    """
    Write table with label column, value columns, and column names.

    Parameters
    ----------
    labels :  list of integers
        label numbers (same length as values)
    columns :  list of lists of floats or integers
        values (each list is a column of values)
    column_names :  list of strings
        names of columns
    table_file : string
        name of output table file

    Returns
    -------
    table_file : string
        name of output table file

    Examples
    --------
    >>> from mindboggle.utils.io_file import write_table
    >>> labels = [0,1,3,5]
    >>> columns = [0.12,0.36,0.75,0.03]
    >>> column_names = ['label', 'volume']
    >>> table_file = 'label_volume_shapes.txt'
    >>> write_table(labels, columns, column_names, table_file)

    """
    import sys

    #-----------------------
    # Check format of inputs
    #-----------------------
    # If the list contains integers or floats, put in a list.
    if isinstance(columns[0], int) or isinstance(columns[0], float):
        columns = [columns]
    # If the list contains all lists, accept format.
    elif all([isinstance(x, list) for x in columns]):
        pass
    else:
        print "io_file.py: Error: columns contains unacceptable elements."
        print "io_file.py: columns type is:", type(columns)
        print "io_file,py: columns length is:", len(columns)
        print "io_file.py: columns[0] type is:", type(columns[0])
        sys.exit()
    # If column_names is a string, create a list containing
    # as many of this string as there are columns.
    if isinstance(column_names, str):
        column_names = [column_names for x in columns]
    elif isinstance(column_names, list):
        if len(column_names) < len(columns):
            column_names = [column_names[0] for x in columns]
        else:
            pass
    else:
        print "Error: column_names is neither a list nor a string"
        sys.exit()

    #----------------------------
    # Open table file for writing
    #----------------------------
    Fp = open(table_file, 'w')

    if column_names:
        Fp.write("\t".join(column_names) + "\n")

    for irow, label in enumerate(labels):
        row = ["{0}\t".format(label)]
        for column in columns:
            row.append(str(column[irow]))
        Fp.write("\t".join(row) + "\n")

    Fp.close()

    return table_file

def write_list(table_file, List, header=""):
    """
    Write a list to a file, each line of which is a list element.
    """

    Fp = open(filename,'w')

    if header:
        Fp.write(header + '\n')

    for Element in List:
        Fp.write(str(Element) + '\n')

    Fp.close()

"""
def write_table_means(filename, column_names, labels, *values):
    ""
    Make a table of mean values per label.

    NOTE:  untested

    Parameters
    ----------
    filename :  output filename (without path)
    column_names :  names of columns [list of strings]
    labels :  list (same length as values)
    *values :  arbitrary number of lists, each containing a value per label

    Returns
    -------
    table_file :  table file

    ""
    import os
    from shapes.measure import mean_value_per_label

    columns = []
    for value_list in values:
        mean_values, label_list = mean_value_per_label(value_list, labels)
        columns.append(mean_values)

    filename = os.path.join(os.getcwd(), filename)
    write_table(label_list, columns, column_names, filename)

    return filename
"""
