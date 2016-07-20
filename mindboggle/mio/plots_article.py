import os
import numpy as np
import pandas as pd

from mindboggle.mio.labels import DKTprotocol

# ----------------------------------------------------------------------------
# Setup:
# ----------------------------------------------------------------------------
# Run first step (concatenate tables):
concat_tables = False
dirnum = 1

# ----------------------------------------------------------------------------
# Alternating left, right cortex label numbers:
# ----------------------------------------------------------------------------
dkt = DKTprotocol()
labels_left = dkt.left_cerebrum_cortex_DKT31_numbers
labels_right = dkt.right_cerebrum_cortex_DKT31_numbers
DKT31_names = dkt.DKT31_names
labels = []
label_names = []
for ilabel, label_left in enumerate(labels_left):
    labels.append(label_left)
    labels.append(labels_right[ilabel])
    label_names.append(DKT31_names[ilabel] + ' (left)')
    label_names.append(DKT31_names[ilabel] + ' (right)')

# ----------------------------------------------------------------------------
# Volume overlap files:
# ----------------------------------------------------------------------------
if dirnum == 1:
    directory = 'volume_overlap_mindboggle_vs_ants_filled_labels'
    file_append = '_ants_filled_labels_overlap.csv'
    title = "Volume overlap: Manual labels vs. Mindboggle with ANTs-filled labels"
elif dirnum == 2:
    directory = 'volume_overlap_mindboggle_vs_freesurfer_filled_labels'
    file_append = '_freesurfer_wmparc_filled_labels_overlap.csv'
    title = "Volume overlap: Manual labels vs. Mindboggle with FreeSurfer-filled labels"
elif dirnum == 3:
    directory = 'volume_overlap_mindboggle_vs_freesurfer_filled_labels_no_ants'
    file_append = '_freesurfer_wmparc_filled_labels_in_freesurfer_segmentation_overlap.csv'
    title = "Volume overlap: Manual labels vs. Mindboggle with FreeSurfer-filled labels (no ANTs segmentation)"
subject_list = '/Users/arno/Data/subject_list_Mindboggle101.txt'
path = '/Users/arno/Desktop/tables_mindboggle_vs_manual_label_overlaps/'
html_file = directory + '.html'
jaccard_file = 'subjects_by_labels_jaccards' + file_append
dice_file = 'subjects_by_labels_dices' + file_append

# ----------------------------------------------------------------------------
# Concatenate volume overlap files:
# ----------------------------------------------------------------------------
path = os.path.join(path, directory)
files = os.listdir(path)
subjects = [x.strip(file_append) for x in files]
if concat_tables:
    jaccards = np.zeros((len(files), len(labels)))
    dices = jaccards.copy()
    for ifile, file in enumerate(files):
        file = os.path.join(path, file)
        columns = pd.read_csv(file)
        for ilabel, label in enumerate(labels):
            for irow in range(columns.size):
                if int(columns.iloc[irow][0].split()[0]) == label:
                    row = columns.iloc[irow][0].split()[1::]
                    jaccards[ifile, ilabel] = np.float(row[0])
                    dices[ifile, ilabel] = np.float(row[1])

    df_jaccards = pd.DataFrame(jaccards, index=subjects, columns=labels)
    df_dices = pd.DataFrame(dices, index=subjects, columns=labels)
    df_jaccards.to_csv(jaccard_file)
    df_dices.to_csv(dice_file)

# ----------------------------------------------------------------------------
# Plot heatmap for subjects X labels array:
# ----------------------------------------------------------------------------
from math import pi
from bokeh.models import HoverTool
from bokeh.plotting import ColumnDataSource, figure, show, output_file
from mindboggle.mio.colors import viridis_colormap

fid = open(subject_list, 'r')
subjects = [x.strip() for x in fid.readlines()]

data_file = dice_file
data = pd.read_csv(data_file)
data = data.set_index('Unnamed: 0')
labels = [str(x) for x in labels]

colors = viridis_colormap()
#from matplotlib import cm as cmaps
#import matplotlib.pyplot as plt
#plt.register_cmap(name='viridis', cmap=cmaps.viridis)
#plt.set_cmap(cmaps.viridis)

# Set up the data for plotting. We will need to have values for every
# pair of subject/label names. Map the value to a color.
subjectx = []
labelx = []
label_namex = []
valuex = []
colorx = []
for ilabel, label in enumerate(labels):
    for isubject, subject in enumerate(subjects):
        labelx.append(label)
        label_namex.append(label_names[ilabel])
        subjectx.append(subject)
        valuex.append(data[label][subject])
        rgb = [np.int(255 * x) for x in
               colors[np.int(255 * data.iloc[isubject, ilabel])]]
        hex = "#%02x%02x%02x" % tuple(rgb)
        colorx.append(hex)

output_file(html_file, title="Volume overlap")
source = ColumnDataSource(dict(subject=subjectx, label=labelx,
                               label_name=label_namex,
                               color=colorx, value=valuex))
TOOLS = "hover,save,pan,box_zoom,wheel_zoom"

p = figure(title=title, x_range=label_names, y_range=subjects,
           x_axis_location="above", plot_width=620, plot_height=1010,
           tools=TOOLS)

p.grid.grid_line_color = None
p.axis.axis_line_color = None
p.axis.major_tick_line_color = None
p.axis.major_label_text_font_size = "5pt"
p.axis.major_label_standoff = 0
p.xaxis.major_label_orientation = pi/3
p.rect(x="label_name", y="subject", width=1, height=1, source=source,
       color="color", line_color=None)

p.select_one(HoverTool).tooltips = [
    ('subject', '@subject'),
    ('label', '@label'),
    ('label name', '@label_name'),
    ('value', '@value'),
]

show(p)      # show the plot

