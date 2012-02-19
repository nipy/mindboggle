# Read in 21 FreeSurfer *.annot files and output one VTK file of all consistent labels for vertexes
 
def read_annot(filepath, orig_ids=False):
    """Read in a Freesurfer annotation from a .annot file.
    From https://github.com/nipy/PySurfer/blob/master/surfer/io.py

    Parameters
    ----------
    filepath : str
    Path to annotation file
    orig_ids : bool
    Whether to return the vertex ids as stored in the annotation
    file or the positional colortable ids
    
    Returns
    -------
    labels : n_vtx numpy array
    Annotation id at each vertex
    ctab : numpy array
    RGBA + label id colortable array
    names : numpy array
    Array of region names as stored in the annot file
    
    """
    import numpy as np
    
    with open(filepath, "rb") as fobj:
        dt = ">i4"
        vnum = np.fromfile(fobj, dt, 1)[0]
        data = np.fromfile(fobj, dt, vnum * 2).reshape(vnum, 2)
        labels = data[:, 1]
        ctab_exists = np.fromfile(fobj, dt, 1)[0]
        if not ctab_exists:
            raise Exception('Color table not found in annotation file')
        n_entries = np.fromfile(fobj, dt, 1)[0]
        if n_entries > 0:
            length = np.fromfile(fobj, dt, 1)[0]
            orig_tab = np.fromfile(fobj, '>c', length)
            orig_tab = orig_tab[:-1]

            names = list()
            ctab = np.zeros((n_entries, 5), np.int)
            for i in xrange(n_entries):
                name_length = np.fromfile(fobj, dt, 1)[0]
                name = np.fromfile(fobj, "|S%d" % name_length, 1)[0]
                names.append(name)
                ctab[i, :4] = np.fromfile(fobj, dt, 4)
                ctab[i, 4] = (ctab[i, 0] + ctab[i, 1] * (2 ** 8) +
                              ctab[i, 2] * (2 ** 16) +
                              ctab[i, 3] * (2 ** 24))
        else:
            ctab_version = -n_entries
            if ctab_version != 2:
                raise Exception('Color table version not supported')
            n_entries = np.fromfile(fobj, dt, 1)[0]
            ctab = np.zeros((n_entries, 5), np.int)
            length = np.fromfile(fobj, dt, 1)[0]
            _ = np.fromfile(fobj, "|S%d" % length, 1)[0] # Orig table path
            entries_to_read = np.fromfile(fobj, dt, 1)[0]
            names = list()
            for i in xrange(entries_to_read):
                _ = np.fromfile(fobj, dt, 1)[0] # Structure
                name_length = np.fromfile(fobj, dt, 1)[0]
                name = np.fromfile(fobj, "|S%d" % name_length, 1)[0]
                names.append(name)
                ctab[i, :4] = np.fromfile(fobj, dt, 4)
                ctab[i, 4] = (ctab[i, 0] + ctab[i, 1] * (2 ** 8) +
                                ctab[i, 2] * (2 ** 16))
        ctab[:, 3] = 255
    if not orig_ids:
        ord = np.argsort(ctab[:, -1])
        labels = ord[np.searchsorted(ctab[ord, -1], labels)]
    return labels, ctab, names 

#def labelMap(Label):
#    '''Given a label as in http://surfer.nmr.mgh.harvard.edu/fswiki/AnnotFiles, return its index
#    
#    '''
#    
#    Map = {1639705:1, 2647065:2, 10511485:3, 6500:4, 3294840:5,\
#           6558940:6, 660700:7, 9231540:8, 14433500:9, 7874740:10,\
#           9180300:11, 9182740:12, 3296035:13, 9211105:14, 4924360:15,\
#           3302560:16, 3988500:17, 3988540:18, 9221340:19, 3302420:20,\
#           1326300:21, 3957880:22, 1316060:23, 14464220:24, 14423100:25,\
#           11832480:26, 9180240:27, 8204875:28, 10542100:29, 9221140:30,\
#           14474380:31, 1351760:32, 6553700:33, 11146310:34, 13145750:35, 2146559:36, 0:0}
#
#    return Map[Label]

def load21(AnnotPath, Subject):
    '''Load 21 annotation files for each of the two hemispheres of a subject
    
    Parameters
    ===========
    
    Subject: string
        The ID an MDD subject
    
    LeftFiles: list of strings
        Each element is an annotation file for the subject's left hemisphere from 21 atlases 

    RightFiles: list of strings
        Each element is an annotation file for the subject's right hemisphere from 21 atlases
        
    AnnotPath: string
        Path where all annotation files are saved
    
    Returns 
    ========
    
    LeftLabels: list of 21 lists of integers
        Each element is a list of labels for all vertexes on Subject's left hemisphere
        
    RightLabels: list of 21 lists of integers
        Each element is a list of labels for all vertexes on Subject's right hemisphere 
    
    Notes
    ======
    
    The labels from read_annot range from 1 to 35 
    
    '''
    
    print "Loading 21 annotations "
    
    from os import listdir
    Allfiles  = listdir(AnnotPath)
    LeftFiles, RightFiles = [], []
    for File in Allfiles:
        if File.find(Subject) > -1:
            if File[0] == 'l':
                LeftFiles.append(File)
            elif File[0] == 'r':        
                RightFiles.append(File)
            else:
                print "unable to match any file names"
    
    LeftLabels, RightLabels = [], []
    for File in LeftFiles:
        Labels, ColorTable, Names = read_annot(AnnotPath+File)
        LeftLabels.append(Labels)

    for File in RightFiles:
        Labels, ColorTable, Names = read_annot(AnnotPath+File)
        RightLabels.append(Labels)
    
    print "21 annotations loaded"
    
    return LeftLabels, RightLabels
    
def vote21(LeftLabels, RightLabels):
    '''For each vertex, let 21 labels of it vote. If they do not all agree, the vertex is unlabeled

    Parameters 
    ========
    LeftLabels: list of 21 lists of integers
        Each element is a list of labels for all vertexes on Subject's left hemisphere
        
    RightLabels: list of 21 lists of integers
        Each element is a list of labels for all vertexes on Subject's right hemisphere
        
    LeftNum: integer
        number of vertexes on the left hemisphere
        
    RightNum: integer
        number of vertexes on the right hemisphere
        
    Returns
    ========
    
    LeftAssign: list of integers
        The ``voting'' result of all vertexes on left hemisphere  
    
    RightAssign: list of integers
        The ``voting'' result of all vertexes on right hemisphere  
        
    LeftConsensus: list of integers
        Number of consensus labels for all vertexes on left hemisphere 
  
    RightConsensus: list of integers
        Number of consensus labels for all vertexes on right hemisphere
  
    '''
    from collections import Counter
    
    print "Begin voting"
    
    LeftNum, RightNum = len(LeftLabels[0]), len(RightLabels[0])
    LeftAssign = [-1 for i in xrange(LeftNum)]   # if no consistent vote, the label is -1, meaning this vertex is unlabeled
    RightAssign = [-1 for i in xrange(RightNum)]
    LeftDiff = [1 for i in xrange(LeftNum)]   
    RightDiff = [1 for i in xrange(RightNum)]
    LeftConsensus = [21 for i in xrange(LeftNum)]   
    RightConsensus = [21 for i in xrange(RightNum)]
    
    for Vrtx in xrange(LeftNum):
        Votes = Counter([LeftLabels[i][Vrtx] for i in xrange(21)])
        if len(Votes) == 1:
            LeftAssign[Vrtx] = Votes.most_common(1)[0][0]
        else: 
            LeftConsensus[Vrtx] = Votes.most_common(1)[0][1]
            LeftDiff[Vrtx] = len(Votes)
            
    for Vrtx in xrange(RightNum):
        Votes = Counter([RightLabels[i][Vrtx] for i in xrange(21)])
        if len(Votes) == 1:
            RightAssign[Vrtx] = Votes.most_common(1)[0][0]
        else:
            RightConsensus[Vrtx] = Votes.most_common(1)[0][1]
            RightDiff[Vrtx] = len(Votes)
            
    print "Voting done"
    return LeftAssign, RightAssign, LeftConsensus, RightConsensus, LeftDiff, RightDiff 

def labeling(SurfPath, Subject, AnnotPath):
    '''Load vtk surfaces from SurfPath, and write assigned labels of Subject into SurfPath as VTK files, 
    according to labels from 21 atlases in AnnotPath
    
    '''
    
    LeftLabels, RightLabels = load21(AnnotPath, Subject)
    LeftAssign, RightAssign, LeftConsensus, RightConsensus, LeftDiff, RightDiff = vote21(LeftLabels, RightLabels)
    
    import pyvtk
    VTKReader = pyvtk.VtkData(SurfPath+"lh.pial.vtk")
    Vertexes =  VTKReader.structure.points
    Faces =     VTKReader.structure.polygons
    pyvtk.VtkData(pyvtk.PolyData(points=Vertexes, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(LeftAssign, name='Assigned_Label'),\
                                  pyvtk.Scalars(LeftDiff, name='Diff_Labels'),\
                                  pyvtk.Scalars(LeftConsensus, name='Common_Labels'))).\
                  tofile(SurfPath+'lh.assign.pial.vtk', 'ascii')
    
    VTKReader = pyvtk.VtkData(SurfPath+"lh.inflated.vtk")
    Vertexes =  VTKReader.structure.points
    pyvtk.VtkData(pyvtk.PolyData(points=Vertexes, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(LeftAssign, name='Assigned_Label'),\
                                  pyvtk.Scalars(LeftDiff, name='Diff_Labels'),\
                                  pyvtk.Scalars(LeftConsensus, name='Common_Labels'))).\
                  tofile(SurfPath+'lh.assign.inflated.vtk', 'ascii')

    VTKReader = pyvtk.VtkData(SurfPath+"rh.pial.vtk")
    Vertexes =  VTKReader.structure.points
    Faces =     VTKReader.structure.polygons
    pyvtk.VtkData(pyvtk.PolyData(points=Vertexes, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(RightAssign, name='Assigned_Label'),\
                                  pyvtk.Scalars(RightDiff, name='Diff_Labels'),\
                                  pyvtk.Scalars(RightConsensus, name='Common_Labels'))).\
                  tofile(SurfPath+'rh.assign.pial.vtk', 'ascii')

    VTKReader = pyvtk.VtkData(SurfPath+"rh.inflated.vtk")
    Vertexes =  VTKReader.structure.points
    pyvtk.VtkData(pyvtk.PolyData(points=Vertexes, polygons=Faces),\
                  pyvtk.PointData(pyvtk.Scalars(RightAssign, name='Assigned_Label'),\
                                  pyvtk.Scalars(RightDiff, name='Diff_Labels'),\
                                  pyvtk.Scalars(RightConsensus, name='Common_Labels'))).\
                  tofile(SurfPath+'rh.assign.inflated.vtk', 'ascii')

# test 
labeling('/forrest/data/MRI/MDD/50332/surf/', '50332', '/forrest/data/MRI/MDD/atlas_to_patients/')

# do some real work
#import os
#for DirIndx, Dir in enumerate(os.listdir('/forrest/data/MRI/MDD')):
#    if len(Dir) == 5 and DirIndx <= 30:
#        print Dir
#        labeling('/forrest/data/MRI/MDD/'+Dir+'/surf/', Dir, '/forrest/data/MRI/MDD/atlas_to_patients/')