#!/usr/local/bin/python

"""
The start of a bulb framework for the MindboggleDB using the Neo4j graph database.
"""
from bulbs.model import Model, Node, Relationship, Resource
from bulbs.property import Property, String, Integer

__author__ = "Nolan Nichols"
__copyright__ = ""
__credits__ = ["Noah Lee", "Arno Klein" "Nolan Nichols"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Nolan Nichols"
__email__ = "bnniii@uw.edu"
__status__ = "Prototype"

Model.resource = Resource('http://50.18.175.117:8182/graphs/mindboggle')

class Database(Node):
    element_type = "database"

    name = Property(String, nullable=False)
            
    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class Project(Node):
    element_type = "project"

    name = Property(String, nullable=False)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class Subject(Node):
    element_type = "subject"

    name = Property(String, nullable=False)
    age = Property(Integer)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class Basins(Node):
    element_type = "basins"

    name = Property(String, nullable=False)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class Basin(Node):
    element_type = "basin"

    name = Property(String, nullable=False)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class hasProject(Relationship):
    label = "has_project"

    @property
    def database(self):
        return Database.get(self.outV)

    @property
    def project(self):
        return Project.get(self.inV)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class hasSubject(Relationship):
    label = "has_subject"

    @property
    def project(self):
        return Project.get(self.outV)

    @property
    def subject(self):
        return Subject.get(self.inV)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class hasBasins(Relationship):
    label = "has_basins"

    @property
    def subject(self):
        return Subject.get(self.outV)

    @property
    def basins(self):
        return Basins.get(self.inV)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name

class hasBasin(Relationship):
    label = "has_basin"

    @property
    def basins(self):
        return Basins.get(self.outV)

    @property
    def basin(self):
        return Basin.get(self.inV)

    def __unicode__(self):
        # include code to create relationships and to index the node
        return self.name