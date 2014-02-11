#!/usr/bin/python
"""
Extracting features from VTK input files. 

Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

For algorithmic details, please check:     
    Forrest S. Bao, et al., Automated extraction of nested sulcus features from human brain MRI data, 
    IEEE EMBC 2012, San Diego, CA

Dependencies:
    python-vtk: vtk's official Python binding
    numpy

Last updated on 2012-09-30: 
    1. Disable secondary surface input and the ability to output features onto a secondary surface.
    2. If users want to visualize features on a secondary surface, check maps_to_second_surface.py in utils folder. 
    
"""

import sys, getopt # this line imports pulic libraries

def print_help():
    print "\n  Usage: python vtk_extract.py  [OPTIONS] InputVTK\n"
    print "  InputVTK"
    print "    The VTK file that contains a map on which features will be extracted (by default, a depth map)\n"
    print "  Options: "
    print "  --thick ThicknessFile"
    print "    the thickness file provided by FreeSurfer"
#    print "  --second SecondSurfaceFile"
#    print "    the 2nd surface to save feature(s) extracted, in FreeSurfer surface format"
    print "  --convex ConvexityFile"
    print "    the convexity file provided by FreeSurfer"
    print "  --fundi FundiVTK"
    print "    the VTK output file to save fundi, onto the same surface as VTKFile"
#    print "  --fundi2 Fundi2VTK"
#    print "    the VTK output file to save fundi, onto the second surface in SecondSurfaceFile"
    print "  --pits PitsVTK"
    print "    the VTK output file to save pits, onto the same surface as VTKFile"
#    print "  --pits2 Pits2VTK"
#    print "    the VTK output file to save pits, onto the second surface in SecondSurfaceFile"
    print "  --sulci SulciVTK"
    print "    the VTK output file to save sulci, onto the same surface as VTKFile"
#    print "  --sulci2 Sulci2VTK"
#    print "    the VTK output file to save sulci, onto the second surface in SecondSurfaceFile"
    print "  --sulciThld SulciThld"
    print "    SulciThld: the value to threshold the surface to get sulci (default = 0.2)"
    print "  --clouchoux [followed by no argument]"
    print "    whether extract pits using Clouchoux's definition"
    print "  --meancurv MeanCurvVTK"
    print "    the VTK input that contains a mean curvature map"
    print "  --gausscurv GaussCurvVTK"
    print "    the VTK input that contains a Gaussian curvature map"
    print "  Examples: "
    print "    python vtk_extract.py --thick lh.thickness --convex lh.sulc\
 --fundi lh.pial.fundi.vtk --pits lh.pial.pits.vtk --sulci lh.pial.sulci.vtk\
 --sulciThld 0.15 --clouchoux --meancurv lh.mean.curv.vtk --gausscurv lh.gauss.curv.vtk lh.depth.vtk\n"

def check_opt(opts, args):
    '''Check whether opts and args satisfy constraints 
    '''
    
    for o, p in opts:
        if o in ('-h', "--help"):
            print_help()
            sys.exit()
    
    ## Check whether mandatory input is provided 
    if len(args)!=1:
        print "  [ERROR]Please give at least one input VTK file containing a per-vertex map"
        print "  Since we use getopt module to process options, please provide options ahead of arguments\n"
        print "  Please run [python vtk_extract.py -h] or [python vtk_extract.py --help] to get the usage."
        sys.exit()
        
    ## Check whether second surface is provided 
#    Need2nd, Given2nd = False, False # A second surface is needed if options to save features into second surface are on.
    
    Clouchoux, MeanCurvVTK, GaussCurvVTK = False, False, False
    for o,p in opts:
        if p == "" and not o in ['--clouchoux', '--help', '-h']:
            print "[ERROR] The option",o, "is missing argument."
            print "   Please run [python vtk_extract.py -h] or [python vtk_extract.py --help] to get the usage."
            sys.exit()
#        elif o in ["--fundi2","--pits2","--sulci2"]:
#            Need2nd = True # we will need the option --second and its parameter
#        elif o=="--second":
#            Given2nd = True
        elif o == "--clouchoux":
            Clouchoux = True
        elif o == "--meancurv":
            MeanCurvVTK = True
        elif o == "gausscurv":
            GaussCurvvTK = True 

#    if Need2nd and not Given2nd: # second surface is needed but not provided
#        print "[ERROR] Second surface is needed but not provided. Please check usage and provide a second surface using --second option."
#        sys.exit()
        
    if Clouchoux and not (MeanCurvVTK or GaussCurvVTK):
        print "[ERROR] Clouchoux type pits will be extracted but mean curvature map and/or Gaussian curvature map is not given."
        sys.exit()
       

def process_opt(opts,args):
    '''Give option parameters to variables that represent I/O file names.
    Do this AFTER check_opt.
    '''
    
    InputColor = "\033[32m"
    OutputColor = "\033[36m"
    EndColor = "\033[0m"
    
    ThickFile, ConvexFile = '', ''
    FundiVTK, PitsVTK, SulciVTK = '','',''
#    Fundi2VTK, Pits2VTK, Sulci2VTK = '','',''
    SulciThld = 0.2
    Clouchoux = False
    MeanCurvVTK, GaussCurvVTK = "", ""
    
    [VTKFile] = args
    print "  [Input] VTK mesh:" + InputColor + VTKFile + EndColor
    
    for o,p in opts:
#        if o=='--second':
#            SurfFile2 = p
#            print "  [Input] second surface:" + InputColor + SurfFile2 + EndColor
        if o=='--convex':
            ConvexFile = p
            print "  [Input] convexity file:" + InputColor + ConvexFile + EndColor
        elif o=='--thick':
            ThickFile = p
            print "  [Input] thickness file:" + InputColor + ThickFile + EndColor
        elif o == '--clouchoux':
                print "  [Settings] Will extract pits using Clouchoux's definition "
                Clouchoux = True
        elif o == "--meancurv":
            MeanCurvVTK = p
            print "  [Input] mean curvature file:" + InputColor + MeanCurvVTK + EndColor
        elif o == "--gausscurv":
            GaussCurvVTK = p
            print "  [Input] Gaussian curvature file:" + InputColor + GaussCurvVTK + EndColor
        elif o =='--sulciThld':
            if p == '':
                print "  [ERROR] Please provide a threshold value for option --", o
                print "To check usage, run: python extract.py"
                sys.exit()
            SulciThld = float(p)
            print "  [Settings] Threshold the surface using value:", SulciThld
        elif o=='--fundi':
            FundiVTK = p
            print "  [Output] Fundi in:" + OutputColor + FundiVTK + EndColor
        elif o=='--pits':
            PitsVTK = p
            print "  [Output] Pits in:" + OutputColor + PitsVTK + EndColor
        elif o=='--sulci':
            SulciVTK = p
            print "  [Output] Sulci in:" + OutputColor + SulciVTK + EndColor
#        elif o=='--fundi2':
#            Fundi2VTK = p
#            print "  [Output] Fundi on 2nd surface in:" + OutputColor +  Fundi2VTK + EndColor
#        elif o=='--pits2':
#            Pits2VTK = p
#            print "  [Output] Pits on 2nd surface in:" +  OutputColor + Pits2VTK + EndColor
#        elif o=='--sulci2':
#            Sulci2VTK = p
#            print "  [Output] Sulci on 2nd surface in:" + OutputColor + Sulci2VTK + EndColor


 


    return VTKFile, ThickFile, ConvexFile, MeanCurvVTK, GaussCurvVTK,\
        FundiVTK, PitsVTK, SulciVTK,\
        SulciThld, Clouchoux

if __name__ == "__main__":
    # this script extracts fundi/pits/sulci from a VTK file from the output from Joachim's code
    
    import libfundi # this line imports my own library

    try:
        opts, args = getopt.getopt(sys.argv[1:],"h",
                                   ["clouchoux", "help","thick=", "convex=","fundi=","pits=","sulci=","sulciThld=","meancurv=", "gausscurv="])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        print "   Please run [python vtk_extract.py -h] or [python vtk_extract.py --help] to get the usage."
        print "If you use options like --second, --pits2, --fundi2, --sulci2,",
        "please note that we have disabled mapping features onto secondary surface.",
        "Please check maps_to_second_surface.py under utils directory to do so."
        sys.exit(2)
    check_opt(opts, args)
    
    InputArguments = process_opt(opts, args)
    
    libfundi.getFeatures(InputArguments, 'vtk','')