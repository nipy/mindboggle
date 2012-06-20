import libfundi, fileio, libvtk # this line imports my own library
import sys, getopt # this line imports pulic libraries

try:
    opts, args = getopt.getopt(sys.argv[1:],'S:C:T:')
except getopt.GetoptError, err:
    print str(err) # will print something like "option -a not recognized"
    print "\n\t Usage: python extract.py [-S 2nd Surface File] [-C 2nd Curvature File] [-T thickness file] CurvFile SurfFile"
    sys.exit(2)

if len(args)!=2:
    print "\n\t Please give a CurvFile and a SurfFile"
    print "\t Usage: python extract.py [-S 2nd Surface File] [-C 2nd Curvature File] [-T thickness file] CurvFile SurfFile"
    sys.exit()

#print args
#print opts

[CurvFile, SurfFile] = args


print "\tRead in FreeSurfer curvature file:", CurvFile
print "\tRead in FreeSurfer surfer file:", SurfFile

SurfFile2, CurvFile2 = '', ''
ThickFile  = ''

#print opts

libvtk.surf2VTK(SurfFile)

for o,p in opts:
    if o in ['-S']:
        SurfFile2 = p
        print "\t Extraction result be mapped onto the 2nd surface:", SurfFile2
        libvtk.surf2VTK(SurfFile2)
    elif o in ['-C']:
        CurvFile2 = p
        print "\t The surface will be thresholded by the 2nd curvature map:", CurvFile2
    elif o in ['-T']:
        ThickFile = p
        print "\t Thickness file in use:", ThickFile


# new unified extraction function Forrest 2011-10-08 
libfundi.getFundi([CurvFile, SurfFile, SurfFile2, CurvFile2, ThickFile], 'FreeSurfer', '')


# extract basin and pits
#libbasin.getBasin(CurvFile, SurfFile, ToVTK, SurfFile2 = SurfFile2)

# extract fundus curves from pits
#libfundipit.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)

# extract fundus curves from fundus face strips
#libfundifc.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)

# extract fundus curves from vertexes directly as Gang Li's
#libfundivtx.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)