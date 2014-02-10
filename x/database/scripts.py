import os

from mbdb.upload import *

# set the project
# TODO need to index all projects with K/V pairs for easy search
#proj = list(Project.get_all())[0]

set_db_url(server="http://50.18.175.117:8182/graphs/mindboggle")

db = create_db('MindBoggleDB')

proj = create_project('MDD', db)

# get a list of the files
dataList = os.listdir('.')

for file in dataList:
    subjectName = file.partition('_')[0]
    subject = create_subject(subjectName, proj)
    stats = read_stats(file)
    set_fundus_stats(subject, stats)