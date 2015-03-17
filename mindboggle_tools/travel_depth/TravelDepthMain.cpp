#include "../TravelDepth.h"
#include "../MeshAnalyser.h"

#include <iostream>
#include <fstream>
#include <stdio.h>


void print_help()
{
    printf(
    "Usage: TravelDepthMain [Options] InputVTKMesh MeanCurvatureOutput\n"
    "Options: \n"
    "   -n: Normalize the output values beteween 0 and 1\n"
    "   -w outputWrapperMesh: export the wrapper mesh\n"
    "Example: TravelDepthMain -n -w wrapper.vtk lh.pial.vtk  lh.travel_depth.vtk \n"
        );
}

int main(int argc, char** argv)
{

    if(argc < 3)
    {
        print_help();
        return -1;
    }

    time_t start= time(NULL);

    bool normalized = false;
    bool wrapperExported = false;
    int wrapperNameId;

    /* Processing options and arguments */
    for (int i=1;i<argc;i++) // we may need to use getopt or getlongopt later - Forrest, 2012/05/29
    {
        if (argv[i][0] != '-')
            continue; // no more options
        switch (argv[i][1])
        {
            case 'n': //normalize the values between 0 and 1.
                normalized = true;
                cout<<"Travel depth will be normalized."<<endl;
                break;
            case 'w': //export the wrapper mesh
                wrapperExported = true;
                i++;
                wrapperNameId = i;
                cout<<"Wrapper will be exported."<<endl;
                break;

            default:
                cout<<"[ERROR] Unrecognized argument flag or option. Check usage. ";
                print_help();
        } // end of switch
    } // end of looping thru arguments

    MeshAnalyser* depthComputer = new MeshAnalyser(argv[argc-2]);
    depthComputer->ComputeTravelDepthFromClosed(normalized);
    depthComputer->WriteIntoFile(argv[argc-1],(char*)"depth");
    
    if(wrapperExported) {
            depthComputer->WriteIntoFile(argv[wrapperNameId],(char*)"closed");
    }
    

    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



