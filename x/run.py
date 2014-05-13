"""
Run Mindboggle on Mindboggle-101 brains
"""
import os

run_all = True

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
ants_dir = '/data/Brains/Mindboggle101/antsCorticalThickness'

for i,name in enumerate(names):
    number = numbers[i]
    for n in range(1,number+1):
        s = 'mindboggle -n 8 {0}-{1} --ants_segments ' \
            '{2}/{0}-{1}/{3}BrainSegmentation.nii.gz ' \
            '--sulci --fundi --thickness --spectra 10 --moments 3 --vertices ' \
            '--surface_labels manual' \
            .format(name, n, ants_dir, prefix)
        print(s)
        os.system(s)
