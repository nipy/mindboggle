"""
Utilities to query the MindboggleDB

    1. Initialize the database
    2. Read csv files containing fundus stats
"""
from bulbs.rest import Resource
from mbdb.base import *
import csv, sys

def read_stats(file):
    """
    read_stats is used to parse the stats files from a mindboggle analysis

    input: a text file and attempts to guess how it is delimited
    output: a list of tuples where each tuple is a row

    """
    filename = file
    with open(filename, 'rb') as f:
        dialect = csv.Sniffer().sniff(f.readline())
        f.seek(0)
        reader = csv.reader(f, dialect)
        header = "".join(reader.next())
        print "Reading %s ..." % header
        try:
            stats = list()
            for row in reader:
                stats.append(tuple(row))
            return stats
        except csv.Error, e:
            sys.exit('file %s, line %d: %s' % (filename, reader.line_num, e))




# fundus_read_stats


# create a root node
mb = Database(name = 'MindboggleDB')

# create a project node
project = Project(name = 'CUMC12')

# connect mb and project
has_proj = hasProject(mb, project)

# create a subject node
subj1 = Subject(name = 'Subject1')

# connect the project to the subject
has_subj = hasSubject(project, subj1)

# create another subject node
subj2 = Subject(name = 'Subject2')

# connect the project to this subject
has_subj = hasSubject(project, subj2)

# create a basins node for the first subject
basins1 = Basins(name = 'Basins')

# connect the subject to their basins
has_basins1 = hasBasins(subj1, basins1)

# create a basin node
basin1 = Basin(name = 'basin')

# connect the basin to the basins
has_basin1 = hasBasin(basins1, basin1)




