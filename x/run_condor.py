"""
Run Mindboggle on Mindboggle-101 brains
"""
import os

run_all = 0 #True

if run_all:
    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20, 21, 22, 20, 1, 1, 2, 2, 12]
else:
    names = ['MMRR-3T7T-2']
    numbers = [2]

prefix = 'ants'
ants_dir = '/data-share/Brains/Mindboggle101/antsCorticalThickness'
condor_file = 'mindboggle101.std'
out_dir = '/data-share'

Fp = open(condor_file, 'wa')
Fp.write('Universe   = standard\n'
         'Executable = mindboggle\n'
         'Log        = mindboggle.log\n'
         'Output     = mindboggle.$(Process).out\n'
         'Error      = mindboggle.$(Process).error\n')

for i,name in enumerate(names):
    number = numbers[i]
    for n in range(1,number+1):

        args = '{0}-{1} ' \
               '--ants_segments {2}/{0}-{1}/{3}BrainSegmentation.nii.gz ' \
               '--out {4}/mindboggled ' \
               '--working {4}/mindboggle_working ' \
               '--cache {4}/mindboggle_cache ' \
               '--sulci --fundi --vertices ' \
               '--thickness --spectra 10 --moments 3 ' \
               '--surface_labels manual' \
               .format(name, n, ants_dir, prefix, out_dir, out_dir, out_dir)
        Fp.write('Arguments = {0}\n'.format(args))
        Fp.write('Queue\n')

Fp.close()
