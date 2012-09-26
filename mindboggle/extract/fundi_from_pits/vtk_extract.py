def print_help():
    print "\n  Usage: python vtk_extract.py  [OPTIONS] InputVTK\n"
    print "  Options: "
    print "  --thick ThicknessFile"
    print "    the thickness file provided by FreeSurfer\n"
    print "  --second SecondSurfaceFile"
    print "    the 2nd surface to save feature(s) extracted, in FreeSurfer surface format\n"
    print "  --convex ConvexityFile"
    print "    the convexity file provided by FreeSurfer\n"
    print "  --fundi FundiVTK"
    print "    the VTK output file to save fundi, onto the same surface as VTKFile\n"
    print "  --fundi2 Fundi2VTK"
    print "    the VTK output file to save fundi, onto the second surface in SecondSurfaceFile\n"
    print "  --pits PitsVTK"
    print "    the VTK output file to save pits, onto the same surface as VTKFile\n"
    print "  --pits2 Pits2VTK"
    print "    the VTK output file to save pits, onto the second surface in SecondSurfaceFile\n"
    print "  --sulci SulciVTK"
    print "    the VTK output file to save sulci, onto the same surface as VTKFile\n"
    print "  --sulci2 Sulci2VTK"
    print "    the VTK output file to save sulci, onto the second surface in SecondSurfaceFile\n"
    print "  --sulciThld SulciThld"
    print "    SulciThld: the value to threshold the surface to get sulci (default = 0.2)"
    print "  Examples: "
    print "    python vtk_extract.py --thick lh.thickness --second lh.inflated --convex lh.sulc --fundi lh.pial.fundi.vtk --pits lh.pial.pits.vtk --sulci lh.pial.sulci.vtk --fundi2 lh.inflated.fundi.vtk --pits2 lh.inflated.pits.vtk --sulci2 lh.inflated.sulci.vtk --sulciThld 0.15 lh.depth.vtk\n"

def check_opt(opts, args):
    '''Check whether opts and args satisfy constraints 
    '''
    ## Check whether mandatory input is provided 
    if len(args)!=1:
        print "\n [ERROR]Please give at least one input VTK file containing a per-vertex map"
        print "  Since we use getopt module to process options, please provide options ahead of arguments"
        print_help()
        sys.exit()
        
    ## Check whether second surface is provided 
    Need2nd, Given2nd = False, False # A second surface is needed if options to save features into second surface are on. 
    for o,p in opts:
        if p == "":
            print "[ERROR] The option",o, "is missing parameter. Please check usage and provide it."
#            print_help()
            sys.exit()
        if o in ["--fundi2","--pits2","--sulci2"]:
            Need2nd = True # we will need the option --second and its parameter
        if o=="--second":
            Given2nd = True   
#    print Need2nd, Given2nd
    if Need2nd and not Given2nd: # second surface is needed but not provided
        print "[ERROR] Second surface is needed but not provided. Please check usage and provide a second surface using --second option."
        sys.exit()

def process_opt(opts,args):
    '''Give option parameters to variables that represent I/O file names.
    Do this AFTER check_opt.
    '''
    
    InputColor = "\033[32m"
    OutputColor = "\033[36m"
    EndColor = "\033[0m"
    
    ThickFile, ConvexFile, SurfFile2 = '', '', ''
    FundiVTK, PitsVTK, SulciVTK = '','',''
    Fundi2VTK, Pits2VTK, Sulci2VTK = '','',''
    SulciThld = 0.2
    
    [VTKFile] = args
    print "  [Input] VTK mesh:" + InputColor + VTKFile + EndColor
    
    for o,p in opts:
        if o=='--second':
            SurfFile2 = p
            print "  [Input] second surface:" + InputColor + SurfFile2 + EndColor
#            libvtk.surf2VTK(SurfFile2)
        elif o=='--convex':
            ConvexFile = p
            print "  [Input] convexity file:" + InputColor + ConvexFile + EndColor
        elif o=='--thick':
            ThickFile = p
            print "  [Input] thickness file:" + InputColor + ThickFile + EndColor
        elif o=='--fundi':
            FundiVTK = p
            print "  [Output] Fundi in:" + OutputColor + FundiVTK + EndColor
        elif o=='--pits':
            PitsVTK = p
            print "  [Output] Pits in:" + OutputColor + PitsVTK + EndColor
        elif o=='--sulci':
            SulciVTK = p
            print "  [Output] Sulci in:" + OutputColor + SulciVTK + EndColor
        elif o=='--fundi2':
            Fundi2VTK = p
            print "  [Output] Fundi on 2nd surface in:" + OutputColor +  Fundi2VTK + EndColor
        elif o=='--pits2':
            Pits2VTK = p
            print "  [Output] Pits on 2nd surface in:" +  OutputColor + Pits2VTK + EndColor
        elif o=='--sulci2':
            Sulci2VTK = p
            print "  [Output] Sulci on 2nd surface in:" + OutputColor + Sulci2VTK + EndColor
        elif o =='--sulciThld':
            if p == '':
                print "  [ERROR] Please provide a threshold value for option --", o
                print "To check usage, run: python extract.py"
                sys.exit()
            SulciThld = float(p)        
    return VTKFile, ThickFile, ConvexFile, SurfFile2, \
    FundiVTK, PitsVTK, SulciVTK, \
    Fundi2VTK, Pits2VTK, Sulci2VTK, SulciThld

if __name__ == "__main__":
    # this script extracts fundi/pits/sulci from a VTK file from the output from Joachim's code
    
    import libfundi # this line imports my own library
    import sys, getopt # this line imports pulic libraries


    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["thick=","second=","convex=","fundi=","fundi2=","pits=","pits2=","sulci=","sulci2=","sulciThld="])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        print_help()
        sys.exit(2)
    
    check_opt(opts, args)
    
    InputFiles = process_opt(opts, args)
    
    
    libfundi.getFeatures(InputFiles, 'vtk','')