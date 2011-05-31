# All functions for extracting fundi as vertex clouds 

import fileio, libvtk, libbasin, libfundivtx
from numpy import mean, std, abs, matrix, zeros, flatnonzero, sign, array, argmin
import os


def getFundi(CurvFile, SurfFile, ToVTK=True, SurfFile2=''):
    '''String up pits into fundus curves
    '''
    print "Stringing up pits into fundus curves"
    
    Vrtx, Fc = fileio.readSurf(SurfFile)

    NbrLst = libbasin.vrtxNbrLst(len(Vrtx), Fc, SurfFile)   
        
    VrtxCmpntFile = CurvFile + '.cmpnt.vrtx'  # need to run libbasin first to get components
    VrtxCmpnt = fileio.loadCmpnt(VrtxCmpntFile) # need to run libbasin first to get components    
    
 # check whether variables changed 
    
    # the code to lineup pits into fundi    
    
    PitsFile = CurvFile + '.pits'
    Pits = libvtk.loadFundiList(PitsFile)
    PSegs = libfundivtx.lineUp(Pits, NbrLst, VrtxCmpnt, Vrtx, CurvFile) # activated 2011-05-29 19:48

    print "Saving fundus curved obtained via stringing up fundus vertexes into VTK files"

    FPits = CurvFile + '.fundi.sgmt.from.pits'
    fileio.writeFundiSeg(FPits, PSegs)
    
    if ToVTK:
        VTKFile = FPits + "." + SurfFile[-1*SurfFile[::-1].find('.'):] + '.vtk'
        libvtk.seg2VTK(VTKFile, SurfFile, FPits)
        if SurfFile2 != '':
            VTKFile = FPits + "." + SurfFile2[-1*SurfFile2[::-1].find('.'):] + '.vtk'
            libvtk.seg2VTK(VTKFile, SurfFile2, FPits)    
    # End of the code to lineup pits into fundi 
