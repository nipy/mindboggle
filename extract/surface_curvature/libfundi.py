# This file contains all functions to extract fundus curves from per-vertex-value (e.g., curvature) map
# Last updated: 2011-08-09 Forrest Sheng Bao 

import readFreeSurfer as rfs

# level 1 functions: functions to get fundi by different methods ---
def fundiFromPits(CurvFile, SurfFile, ToVTK=True, SurfFile2=''):
    '''String up pits into fundus curves
    '''
    print "Connecting pits into fundus curves"

    Curvature = rfs.readCurv(CurvFile)
    CurvDisp = []
    for x in Curvature:
        if x > 0:
            CurvDisp.append(x)
        else:
            CurvDisp.append(0)

    Vrtx, Fc = rfs.readSurf(SurfFile)

    NbrLst = libbasin.vrtxNbrLst(len(Vrtx), Fc, SurfFile)   
        
    VrtxCmpntFile = CurvFile + '.cmpnt.vrtx'  # need to run libbasin first to get components
    VrtxCmpnt = rfs.loadCmpnt(VrtxCmpntFile) # need to run libbasin first to get components    
    
    # the code to lineup pits into fundi    
    
    PitsFile = CurvFile + '.pits'
    Pits = libvtk.loadFundiList(PitsFile)
    
    PSegs = libfundivtx.lineUp(Pits, NbrLst, VrtxCmpnt, Vrtx, CurvFile, Curvature) # changed 2011-07-21 00:23

    print "Saving fundus curved obtained via stringing up fundus vertexes into VTK files"

    FPits = CurvFile + '.fundi.from.pits'
    fileio.writeFundiSeg(FPits, PSegs)
    
    if ToVTK:
        VTKFile = FPits + "." + SurfFile[-1*SurfFile[::-1].find('.'):] + '.vtk'
        libvtk.seg2VTK(VTKFile, SurfFile, FPits, LUT=[CurvDisp], LUTname=['Curvature'])
        if SurfFile2 != '':
            VTKFile = FPits + "." + SurfFile2[-1*SurfFile2[::-1].find('.'):] + '.vtk'
            libvtk.seg2VTK(VTKFile, SurfFile2, FPits, LUT=[CurvDisp], LUTname=['Curvature'])    
    # End of the code to lineup pits into fundi 


# End of level 1 functions ---

# Level 2 functions: functions that connects special vertexes into curves  ---

# End of level 2 functions ---

# level 3: extracting special vertexes ---
# end of level 3 ---