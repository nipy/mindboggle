"""Utilities to query the MindboggleDB
"""
from bulbs.model import Model
from db.models.base import *

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




