def print_help():
    print "\n  Usage: python extract.py [OPTIONS] SurfFile\n"
    print "  Options (at least one of --curv and --convex are needed): "
    print "  --thick ThicknessFile"
    print "    the thickness file provided by FreeSurfer\n"
    print "  --second SecondSurfaceFile"
    print "    the 2nd surface to save feature(s) extracted, in FreeSurfer surface format\n"
    print "  --conv ConvexityFile"
    print "    the convexity file provided by FreeSurfer\n"
    print "  --curv CurvatureFile"
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
    print "  --use MapToExtractFeatures"
    print "    Choose a map to use for extracting features. Other maps will be used as shape descriptors associated with features only. Values as follows."
    print "    conv: use convexity map (usually of suffix .sulc in FreeSurfer) to extract features (default)"
    print "    curv: use curvatiure map (usually of suffix .curv in FreeSurfer) to extract features"
    print "  --sulciThld SulciThld"
    print "    SulciThld: the value to threshold the surface to get sulci (default = 0.2)"
    print "  Examples: "
    print "    python extract.py --use conv --thick lh.thickness --second lh.inflated --conv lh.sulc --curv lh.curv --fundi lh.pial.fundi.vtk --pits lh.pial.pits.vtk --sulci lh.pial.sulci.vtk --fundi2 lh.inflated.fundi.vtk --pits2 lh.inflated.pits.vtk --sulci2 lh.inflated.sulci.vtk --sulciThld 0.15 lh.pial \n"

def check_opt(opts, args):
    '''Check whether opts and args satisfy constraints 
    '''
    ## Check whether mandatory input is provided 
    if len(args)!=1:
        print "\n [ERROR] At least 1 FreeSurfer surface is needed. "
        print "  Since we use getopt module to process options, please provide options ahead of arguments"
        print "To check usage, run: python extract.py"
        print_help()
        sys.exit()
        
    ## Check whether second surface is provided 
    Need2nd, Given2nd = False, False # A second surface is needed if options to save features into second surface are on. 
    # check whether at least one of --curv and --convex is provided; check wether curv or convex is given according to --use flag parameter
    GivenCurv, GivenConvex = False, False
    
    Use = 'conv'
    
    for o,p in opts:
        if p == "":
            print "[ERROR] The option",o, "is missing parameter. Please check usage and provide it."
            print "To check usage, run: python extract.py"
            sys.exit()
        elif o in ["--fundi2","--pits2","--sulci2"]:
            Need2nd = True # we will need the option --second and its parameter
        elif o=="--second":
            Given2nd = True   
        elif o=="--curv": 
            GivenCurv = True
        elif o=="--conv":
            GivenConvex = True
        elif o=="--use":
            Use=p
#    print Need2nd, Given2nd
    if Need2nd and not Given2nd: # second surface is needed but not provided
        print "[ERROR] Second surface is needed but not provided. Please check usage and provide a second surface using --second option."
        print "To check usage, run: python extract.py"
        sys.exit()
    if not GivenCurv and not GivenConvex:
        print "[ERROR] At least one map in FreeSurfer curvature format needed after --curv or --convex is needed."
        print "To check usage, run: python extract.py"
    if Use=="conv" and not GivenConvex:
        print "[ERROR] Convexity map is specified or by default set for feature extraction but not given."
        print "To check usage, run: python extract.py"
        sys.exit()
    if Use=="curv" and not GivenCurv:
        print "[ERROR] Curvature map is specified for feature extraction but not given."
        print "To check usage, run: python extract.py"
        sys.exit()

def process_opt(opts, args):
    """Give option parameters to variables that represent I/O file names.
    Do this AFTER check_opt.
    """

    InputColor = "\033[32m"
    OutputColor = "\033[36m"
    EndColor = "\033[0m"
    
    ThickFile, ConvFile, CurvFile, SurfFile2 = '', '', '', ''
    FundiVTK, PitsVTK, SulciVTK = '','',''
    Fundi2VTK, Pits2VTK, Sulci2VTK = '','',''
    Use = 'conv'
    SulciThld = 0.2

    [SurfFile] = args
    print "  [Input] FreeSurfer surface file:" + InputColor + SurfFile  + EndColor 
    
    for o,p in opts:
        if o == '--second' :
            SurfFile2 = p
            print "  [Input] 2nd surface:" + InputColor + SurfFile2 + EndColor
        elif o == '--conv':
            ConvFile = p
            print "  [Input] Convexity File:" + InputColor + ConvFile + EndColor
        elif o == '--thick':
            ThickFile = p
            print "  [Input] Thickness File:" + InputColor + ThickFile + EndColor
        elif o == '--curv':
            CurvFile = p
            print "  [Input] Curvature File:" + InputColor + CurvFile + EndColor
        elif o == '--fundi':
            FundiVTK = p
            print "  [Output] Fundi in:" + OutputColor + FundiVTK + EndColor
        elif o == '--pits':
            PitsVTK = p
            print "  [Output] Pits in:", OutputColor, PitsVTK, EndColor
        elif o == '--sulci':
            SulciVTK = p 
            print "  [Output] Sulci in:", OutputColor, SulciVTK, EndColor
        elif o == '--fundi2':
            Fundi2VTK = p
            print "  [Output] Fundi on 2nd surface in:", OutputColor, Fundi2VTK, EndColor
        elif o == '--pits2':
            Pits2VTK = p
            print "  [Output] Pits on 2nd surface in:", OutputColor, Pits2VTK, EndColor
        elif o == '--sulci2':
            Sulci2VTK = p 
            print "  [Output] Sulci on 2nd surface in:", OutputColor, Sulci2VTK, EndColor
        elif o == '--use':
            Use = p
            if Use == 'curv':
                print "  [Method] Extract features on curvature map"
            elif Use == 'conv':
                print "  [Method] Extract features on convexity map"
            else:
                print "  [Error] --use parameter",p, "is unrecognized. Check usage."
                print "To check usage, run: python extract.py" 
        elif o =='--sulciThld':
            if p == '':
                print "  [ERROR] Please provide a threshold value for option --", o
                print "To check usage, run: python extract.py"
                sys.exit()
            SulciThld = float(p)
        else:
            print "[ERROR] Option", o, "is unrecognized. Please check usage."
            sys.exit()

    return SurfFile, SurfFile2, ThickFile, CurvFile, ConvFile,\
        FundiVTK, PitsVTK, SulciVTK, Fundi2VTK, Pits2VTK, Sulci2VTK, Use, SulciThld

if __name__ == "__main__":
    import sys, getopt # this line imports pulic libraries
    import libfundi

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["thick=","second=","conv=","fundi=","fundi2=", "pits=","pits2=","sulci=","sulci2=", "curv=", "use=", "sulciThld="])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
    #    print_help()
        print "To check usage, run: python extract.py"
        sys.exit(2)
    
    check_opt(opts, args)
    InputFiles = process_opt(opts, args)
    
    # new unified extraction function Forrest 2011-10-08 
    libfundi.getFeatures(InputFiles, 'FreeSurfer', '')