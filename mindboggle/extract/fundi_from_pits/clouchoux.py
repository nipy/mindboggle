# A prototyping code to extract pits in Clouchoux's approach

def load_curvs(File):
    '''Load mean curvature and Gaussian curvature from a VTK file
    
    PointData
    1. Norm
    2. Depth  (now travel depth)
    3. mean curvature 
    4. Gaussian curvature 
    
    ''' 
    import pyvtk
    VTKReader = pyvtk.VtkData(File)
    Vertexes =  VTKReader.structure.points
    Faces =     VTKReader.structure.polygons
    Depth =     VTKReader.point_data.data[1].scalars
    MCurv = VTKReader.point_data.data[2].scalars
    GCurv = VTKReader.point_data.data[3].scalars
    
    return Vertexes, Faces, Depth, MCurv, GCurv

def clouchoux(MCurv, GCurv):
    '''Judge whether a vertex is a pit in Clouchoux's definition
    
    Parameters
    ===========
    
    MCurv: float
        mean curvature of a vertex
        H in Clouchoux's paper 
    
    GCurv: float
        mean curvature of a vertex 
        K in Clouchoux's paper
      
    Returns 
    ========
    
    True if this is a pit. False, otherwise. 
    
    Notes
    =========
    
    (Since Joachim's code updates all the time, this settings has to be updated frequently)
    
    In Clochoux's paper, the following definitions are used:
        H > 0, K > 0: pit, in Clouchoux's paper
        H < 0, K > 0: peak, in Clouchoux's paper    
    
    If features are computed by ComputePricipalCurvature(), 
    use this settings to get proper pits:
        H > 3, K < 0 (curvatures not normalized)
        H > 0.2, K < 0 (curvatures normalized)
    	
    
    '''
    
    if (MCurv > 3) and (GCurv < 0):
        return True
    else:
        return False
    
def clouchoux_pits(Vertexes, MCurv, GCurv):
    '''Extract pits using Clouchoux's definition
    '''
    Pits = []

    for i in xrange(len(Vertexes)):
        if clouchoux(MCurv[i], GCurv[i]):
            Pits.append(i)
    
    print  len(Pits), "Pits found"
    
    return Pits
    
def write_Clouchoux_Pits(Vertexes, Pits, File, File2=""):
    '''Write Clouchoux-type pits into VTK
    
    Parameters
    ===========
    
    Vertexes: list of 3-tuples of floats
        Coordinates of all points on the mesh
        
    Pits: list of integers
        IDs of vertexes that are pits by Clouchoux's definition
    
    File: string
        The path to write extracted pits
        
    File2: string
        The path to the second surface
    
    '''
    import pyvtk
    pyvtk.VtkData(pyvtk.PolyData(points=Vertexes, vertices=Pits)).tofile(File, 'ascii')
    
    if File2!="":
        VTKReader = pyvtk.VtkData(File2)
        Vertexes =  VTKReader.structure.points
        pyvtk.VtkData(pyvtk.PolyData(points=Vertexes, vertices=Pits)).tofile(File[:-3]+"inflated.vtk", 'ascii')
    
    return 0
    
def pits_extract(File):
    '''Load a VTK file, extract pits and write pits into another VTK file 
    '''
    
    Vertexes, Faces, Depth, MCurv, GCurv = load_curvs(File)
    Pits = clouchoux_pits(Vertexes, MCurv, GCurv)
    write_Clouchoux_Pits(Vertexes, Pits, File[:-3]+"pits.vtk", File[:File.find(".")]+".inflated.vtk")
#    write_Pits(Vertexes, Pits, 'test_pits.vtk')
        
    return 0 

#import sys
#pits_extract(sys.argv[1])
#pits_extract("/lh.clouchoux.vtk")