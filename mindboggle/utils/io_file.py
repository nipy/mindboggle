#!/usr/bin/python
"""
This Python library reads and writes different file types.

Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
import re

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
    #import re

    Fp = open(filename, 'r')
    lines = Fp.readlines()
    columns = [[] for x in range(n_columns)]
    for line in lines:
        if len(line) > 0:
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

def write_table(labels, columns, column_names, filename):
    """
    Write table with label column, value column, and column names.

    Parameters
    ----------
    labels :  list (same length as values)
    columns :  list of lists of values (each list is a column of values)
    column_names :  names of columns [list of strings]

    """

    Fp = open(filename,'w')

    if len(column_names):
        Fp.write("\t".join(column_names) + "\n")

    for irow, label in enumerate(labels):
        row = ["{0}\t".format(label)]
        for column in columns:
            row.append(str(column[irow]))
        Fp.write("\t".join(row) + "\n")

    Fp.close()

def write_list(filename, List, header=""):
    """
    Write a list to a file, each line of which is a list element.
    """

    Fp = open(filename,'w')

    if len(header):
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
    from measure.measure_functions import mean_value_per_label

    columns = []
    for value_list in values:
        mean_values, label_list = mean_value_per_label(value_list, labels)
        columns.append(mean_values)

    filename = os.path.join(os.getcwd(), filename)
    write_table(label_list, columns, column_names, filename)

    return filename
"""
