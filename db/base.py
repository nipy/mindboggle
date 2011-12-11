"""
Base model of MindboggleDB

Base is a generic set of classes that model the vertices (nodes) and
edges (arcs) in the Mindboggle graph database implementation.

 Domain Objects
    Database
    Project
    Subject
    Basin
    SulcalSurface
    Fundus
    Pit

"""

from bulbs.model import Model, Node, Relationship, Resource
from bulbs.property import Property, String, Integer

class Database(Node):
    """
    Database is the root node of mbdb Base model

    Relationships
        has-a Project, Owner, User, Location
    """
    element_type = "db"

    name = Property(String, nullable=False)

    def __unicode__(self):
        return self.name

class Project(Node):
    """
    Project is the concept of a collection of participants in a study - potentially with a set of
    overlapping metadata attributes

    Relationships
        has-a Subject, labelingProtocol, PI, IRB
    """
    element_type = "project"

    name = Property(String, nullable=False)

    def __unicode__(self):
        return self.name

class Subject(Node):
    """
    Subject is the concept of a participant in a Project with a set of data collected about them

    Relationships
        has-a sulcalBasin, sulcalSurface, Fundus, sulcalPit
        is-a Person
    """
    element_type = "subject"

    name = Property(String, nullable=False)
    age = Property(Integer)

    def __unicode__(self):
        return self.name, self.age

class Basin(Node):
    """
    Basin is anatomical entity or image feature?

    Relationships
        has-a sucalSurface, Fundus, sulcalPit
    """
    element_type = "basin"

    name = Property(String, nullable=False)

    def __unicode__(self):
        return self.name

class hasProject(Relationship):
    """
    hasProject is...
    """
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
    """
    hasSubject is...
    """
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

class hasBasin(Relationship):
    """
    hasBasin is...
    """
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