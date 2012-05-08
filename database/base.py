"""
Base model of MindboggleDB

Base is a generic set of classes that model the vertices (nodes) and
edges (arcs) in the Mindboggle graph database implementation.

 Domain Objects
    Database -
    Project - a collection of subjects
    Person - an individual with a roll
    Subject - a participant in a project -
    Sulcus - a surface (mesh) is a fold of the brain.
    Ribbon - a medial surface (mesh) within a fold, extending from a fundus.
    Fundus - a curve (polyline) runs through the pits (via a minimum spanning tree algorithm).
    Pit - a point (vertex) of maximal depth or curvature within a neighborhood on a brain surface.

    Not Implemented (could follow the XCEDE data model):
        Subject_Group
        Visit
        Study
        Episode
        Acquisition
"""
from bulbs.model import Node, Relationship
from bulbs.property import Property, String, Integer, Float

#  Base Node and Relationship class for MBDB
class NodeMB(Node):
    """
    NodeMB is the root node for all vertices in MBDB
    """
    element_type = "node"

    def __unicode__(self):
        return self.element_type


class RelationshipMB(Relationship):
    """
    RelationshipMB is the root node for all vertices in MBDB
    """
    element_type = "relationship"

    def __unicode__(self):
        return self.element_type

# Vertices
class Database(NodeMB):
    """
    Database is the root node of mbdb domain model
    """
    element_type = "database"
    name = Property(String, nullable=False)

    #def after_initialized(self):
    #    self.create_index(index_keys = self.name)

    #def after_created(self):
    #   self.index.put_unique(self.eid, key = self.element_type, value = self.name)

    def __unicode__(self):
        return self.name


class Project(Database):
    """
    Project is the concept of a collection of participants in a study - potentially with a set of
    overlapping metadata attributes

    Relationships
        contained_in Database
    """
    element_type = "project"
    name = Property(String, nullable=False)

    def after_created(self):
        self.create_index(index_keys=self.name)

    def __unicode__(self):
        return self.name


class Person(NodeMB):
    """
    Project is the concept of a collection of participants in a study - potentially with a set of
    overlapping metadata attributes

    Relationships
        contained_in Project
    """
    element_type = "project"
    name = Property(String, nullable=False)

    def __unicode__(self):
        return self.name


class Subject(Project, Person):
    """
    Subject is the concept of a participant in a Project with a set of data collected about them

    Relationships
        contained_in Project
        is-a Person
    """
    element_type = "subject"

    name = Property(String, nullable=False)
    age = Property(Integer)

    def __unicode__(self):
        return self.name


class Sulcus(Subject):
    """
    Sulcus is anatomical entity with a set of image features

    Relationships
        contained_in Subject
    """
    element_type = "sulcus"

    name = Property(String, nullable=False)

    def __unicode__(self):
        return self.name


class Ribbon(Sulcus):
    """
    Sulcus is anatomical entity with a set of image features

    Relationships
        contained_in Subject
    """
    element_type = "ribbon"

    name = Property(String, nullable=False)

    def __unicode__(self):
        return self.name


class Fundus(Sulcus):
    """
    Sulcus is anatomical entity with a set of image features

    Relationships
        contained_in Subject
    """
    element_type = "fundus"

    name = Property(String, nullable=False)

    curvature = Property(Float)
    convexity = Property(Float)
    depth = Property(Float)
    thickness = Property(Float)
    length = Property(Float)


    def __unicode__(self):
        return self.name


class Pit(Fundus):
    """
    Pit is anatomical entity with a set of image features

    Relationships
        contained_in Fundus
    """
    element_type = "pit"

    name = Property(String, nullable=False)

    def __unicode__(self):
        return self.name

# Relationship types
class ContatinedIn(RelationshipMB):
    """
    ContainedIn is a relationship type

    Usage:
        The first vertex is contained within the second vertex

        Example:
        subject = Subject(name = "Nolan Nichols")
        project = Project(name = "Awesome Project")
        relationship = ContainedIn(Subject,Project)
    """
    label = "contained_in"
    name = "contained_in"

    # incoming node
    @property
    def outVObject(self):
        return self.label

    @outVObject.setter
    def outVObject(self, value):
        outVertex = value.__class__
        outVertex.get(self.outV)

    # outgoing node
    @property
    def inVObject(self):
        return self.label

    @inVObject.setter
    def inVObject(self, value):
        inVertex = value.__class__
        inVertex.get(self.inV)

    def __unicode__(self):
        return self.label


class IsA(Relationship):
    """
    IsA is a relationship type...
    """
    label = "is_a"
    name = "is_a"

    # incoming node
    @property
    def outVObject(self):
        return self.label

    @outVObject.setter
    def outVObject(self, value):
        outVertex = value.__class__
        outVertex.get(self.outV)

    # outgoing node
    @property
    def inVObject(self):
        return self.label

    @inVObject.setter
    def inVObject(self, value):
        inVertex = value.__class__
        inVertex.get(self.inV)

    def __unicode__(self):
        return self.label
