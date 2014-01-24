
import os

names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20']
numbers = [20,21,22,20]
prefix = 'ants'
ants_dir = '/data/Brains/Mindboggle101/antsCorticalThickness' #OASIS-TRT-20_volumes/OASIS-TRT-20-11

for i,name in enumerate(names):
    number = numbers[i]
    for n in range(1,number+1):
        s = 'mindboggle -n 8 {0}-{1} --ants_prefix {2}/{0}_volumes/{0}-{1}/ants --sulci --fundi --vertices --spectra 10 --thickness'.format(name, n, '/data/Brains/Mindboggle101/antsCorticalThickness', 'ants')
        #s = 'mindboggle -n 8 {0}-{1} --ants_prefix {2}/{0}_volumes/{0}-{1}/ants'.format(name, n, ants_dir)
        print(s)
        os.system(s)
