"""
Run Mindboggle on Mindboggle-101 brains
"""
import os

run_all = True
#run_manual = ' '
run_manual = ' --surface_labels manual '
run_proc = ' --proc 20 '
#run_proc = ' --cluster '

if run_all:
    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20,1,1,2,2,12]
else:
    #    name = 'MMRR-3T7T-2'
    #    numbers = [1]
    name = 'HLN-12'
    numbers = range(1,13)

prefix = 'ants'
ants_dir = '/mnt/nfs-share/brains/Mindboggle101/antsCorticalThickness'
subjects_dir = '/mnt/nfs-share/brains/Mindboggle101/subjects'

for i,name in enumerate(names):
    number = numbers[i]
    for n in range(1,number+1):
        s = 'mindboggle {0}-{1} --all {2} {3} --ants ' \
            '{4}/{0}-{1}/{5}BrainSegmentation.nii.gz ' \
            '--subjects_dir /mnt/nfs-share/brains/Mindboggle101/subjects' \
            .format(name, n, run_manual, run_proc, ants_dir, prefix)
        print(s)
        os.system(s)
