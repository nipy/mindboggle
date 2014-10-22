"""
Copy Mindboggle tables for Mindboggle-101 brains to a separate directory.
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

mb_dir = '/homedir/mindboggled'
out_dir = '/homedir/tables'

s = 'mkdir {0}'.format(out_dir)
print(s); os.system(s)

for i,name in enumerate(names):
    number = numbers[i]
    for n in range(1,number+1):
        s = 'mkdir {0}/{1}-{2}'.format(out_dir, name, n)
        print(s); os.system(s)
        s = 'cp -r {0}/{1}-{2}/tables/* {3}/{1}-{2}/'.format(mb_dir, name, n, out_dir)
        print(s); os.system(s)
