
# Comment out the top line (below the header) with the number of rows,
# so that an appropriately sized numpy array is loaded.

label_file = '/home/arno/Documents/Projects/mindboggle/results/workingdir/Mindboggle_workflow/Atlas_workflow/_hemi_lh_subject_test/Write_label_files/mapflow/_Write_label_files30/lh.insula.label'

import numpy as np

a = np.loadtxt(label_file)
x = a[:,1]
y = a[:,2] 
z = a[:,3]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z, 'b.')
plt.show()
