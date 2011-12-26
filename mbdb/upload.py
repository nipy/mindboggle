"""
Utilities for uploading data to MindboggleDB

    1. Initialize the database
    2. Read csv files containing fundus stats
"""
from bulbs.rest import Resource
from mbdb.base import *
import csv, sys

def set_db_url(server=None):
    """
    set_db_url is used to configure the default server for all mbdb.base domain objects

    http://50.18.175.117:8182/graphs
    """
    NodeMB.resource = Resource(db_url=server)
    RelationshipMB.resource = Resource(db_url=server)
    print "db_url set to", server


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
            # stats.append(f.name)
            for row in reader:
                stripped = [float(value) for value in row if len(value) != 0]
                stats.append(tuple(stripped))
            return stats
        except csv.Error, e:
            sys.exit('file %s, line %d: %s' % (filename, reader.line_num, e))


def create_db(db_name):
    db = Database(name=db_name)
    return db


def create_project(project_name, db):
    project = Project(name=project_name)
    ContatinedIn(project, db)
    return project


def create_subject(subject_name, project):
    subject = Subject(name=subject_name)
    ContatinedIn(subject, project)
    return subject


def create_fundus(fundus_name, subject):
    fundus = Fundus(name=fundus_name)
    ContatinedIn(fundus, subject)
    return fundus


def set_fundus_stats(subject, stats):
    for featureVector in range(len(stats)):
        featureVectorName = subject.name + '-' + str(featureVector)
        fundus = create_fundus(featureVectorName, subject)
        fundus.curvature = stats[featureVector][0]
        fundus.convexity = stats[featureVector][1]
        fundus.depth = stats[featureVector][2]
        fundus.thickness = stats[featureVector][3]
        fundus.length = stats[featureVector][4]
        fundus.save()
        print stats[featureVector]




