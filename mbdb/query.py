"""Utilities to query the MindboggleDB
"""

### currently all below is using generic bulbs framework 

from bulbs.graph import Graph
from bulbs.model import Model

# connect to mbdb
g = Graph(db_url="http://localhost:8182/graphs/mindboggle")

# create root node
mb = g.vertices.create({'name':'MindboggleDB',
                        'url':"http://localhost:8182/graphs/mindboggle"})

# create a project node
def create_project(project):
    g.vertices.create({'name':project,
                       'label':'Project'})
    g.edges.create(mb,"has_project",proj)

# create a subject node
subj = g.vertices.create({'name':'Subject1',
                          'label':'Subject'})
g.edges.create(proj,"has_subject", subj)

# create a basins node
basins = g.vertices.create({'name':'Basins',
                           'label':'Basins'})
g.edges.create(subj,'has_basins',basins)

# create a basin node
basin = g.vertices.create({'name':'Basin',
                           'label':'Basin1'})
g.edges.create(basins,'has_basin',basin)
