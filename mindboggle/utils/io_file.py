#!/usr/bin/env python
"""
This Python library reads and writes different file types.

Authors:
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net
    - Arno Klein, 2012  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def read_columns(filename, n_columns=1, trail=False):
    """
    Read n-column text file. Assumes space(s) as delimiter.

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

def write_columns(columns, column_names, output_table, delimiter=',',
                  quote=True, input_table=''):
    """
    Write table with columns and column names.  Assumes space(s) as delimiter.

    If there is an input table file to append to, assume a 1-line header.

    Parameters
    ----------
    columns :  list of lists of floats or integers
        values (each list is a column of values)
    column_names :  list of strings
        names of columns
    output_table : string
        name of output table file
    delimiter : string
        delimiter between columns, such as ','
    bracket : string
        string bracketing each element, such as '"'
    input_table : string (default is empty string)
        name of table file to which the columns are to be appended

    Returns
    -------
    output_table : string
        name of output table file

    Examples
    --------
    >>> from mindboggle.utils.io_file import write_columns
    >>> labels = ['category one', 'category two', 'category three', 'category four']
    >>> values = [0.12, 0.36, 0.75, 0.03]
    >>> values2 = [32, 87, 53, 23]
    >>> columns = [labels, values]
    >>> column_names = ['label', 'value']
    >>> output_table = 'write_columns.csv'
    >>> delimiter = ','
    >>> quote = True
    >>> input_table = ''
    >>> write_columns(columns, column_names, output_table, delimiter, quote, input_table)
    >>> write_columns(values2, 'value 2', output_table, delimiter,
    >>>               quote, input_table=output_table)

    """
    import os
    import sys
    from mindboggle.utils.io_file import read_columns

    output_table = os.path.join(os.getcwd(), output_table)
    if quote:
        q = '"'
    else:
        q = ''

    #-----------------------
    # Check format of inputs
    #-----------------------
    # If the list contains integers or floats, put in a list:
    if columns:
        if isinstance(columns[0], int) or isinstance(columns[0], float) or \
           isinstance(columns[0], str):
            columns = [columns]
        # If the list contains all lists, accept format:
        elif all([isinstance(x, list) for x in columns]):
            pass
        else:
            print("Error: columns contains unacceptable elements.")
            print("columns type is: {0}".format(type(columns)))
            print("columns length is: {0}".format(len(columns)))
            print("columns[0] type is: {0}".format(type(columns[0])))
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
            print("Error: column_names is neither a list nor a string")
            sys.exit()

        #-----------------------------------
        # Read columns from input table file
        #-----------------------------------
        if input_table:
            input_columns = read_columns(input_table, n_columns=1, trail=True)
            input_names = input_columns[0][0]
            input_columns = input_columns[0][1::]
        #else:
        #    input_names = ''
        #    input_columns = ['' for x in columns[0]]

        #--------------
        # Write to file
        #--------------
        Fp = open(output_table, 'wa')
        if column_names:
            column_names = [q+x+q for x in column_names]
            if input_table:
                Fp.write(delimiter.join([input_names,
                                         delimiter.join(column_names) + "\n"]))
            else:
                Fp.write(delimiter.join(column_names) + "\n")
        #else:
        #    Fp.write(input_names + "\n")

        for irow in range(len(columns[0])):
            if input_table:
                Fp.write(input_columns[irow] + delimiter)
            for icolumn, column in enumerate(columns):
                if icolumn < len(columns)-1:
                    Fp.write('{0}{1}{2}{3}'.format(
                        q, column[irow], q, delimiter))
                else:
                    Fp.write('{0}{1}{2}'.format(q, column[irow], q))
            Fp.write("\n")

        Fp.close()

    else:
        print("NOTE: 'columns' is empty. Nothing written.")

    return output_table

def write_rows(filename, list_of_lines, header=""):
    """
    Write a list to a file, one line per list element.

    Parameters
    ----------
    filename : string
        name of output file
    list_of_lines :  list
        each element is written to file as a line
    header : string (default is empty string)
        header to write at the top of the file

    Returns
    -------
    filename : string
        name of output file

    """

    Fp = open(filename, 'w')

    if header:
        Fp.write(header + '\n')

    for element in list_of_lines:
        Fp.write(str(element) + '\n')

    Fp.close()

    return filename
