
import os

subject = 'test'
hemisphere = 'lh'
surface = 'pial'
annotname = 'labels.max'

args = ['tksurfer', subject, hemisphere, surface, '-annotation', annotname]
c = ' '.join(args)
os.system(c)
