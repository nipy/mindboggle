#include "../MeshAnalyser.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

void print_help()
{
    printf(
    "Usage: TravelDepthMain [Options] InputVTKMesh DepthOutput\n"
    "Options: \n"
    "   -d Depth: choose which depth to compute \n"
    "       0 -- travel depth (default) \n"
    "       1 -- Euclidean depth \n"
    "Example: TravelDepthMain -d 1 lh.pial.vtk lh.euc_depth.vtk \n"
        );
}


int main(int argc, char** argv)
{
    time_t start= time(NULL);

    if (argc<2)
    {
         print_help();
         return -1;
    }

    /* Defautl values for options */
    int Depth=0; 
    bool EucDepth=false; 

    /* Select values for options*/
    for (int i=1; i<argc; i++)
    {
        if (argv[i][0] != '-') continue; // no more options
        switch (argv[i][1])
        {
            case 'd' :  // select curvature computation method 
                Depth = atoi(argv[i+1]);
                break;
            
            default: //wrong option
                cout<<"unrecognized option "<<argv[i]<<endl;
                print_help();
                return -2;

        }
    }

    MeshAnalyser* ma = new MeshAnalyser(argv[argc-2]);

    /* Select depth to compute */
    switch (Depth)
    {
        case 0: 
            ma->ComputeTravelDepthFromClosed(false);
            ma->WriteIntoFile(argv[argc-1], (char*)"depth");
        case 1:
            ma->ComputeTravelDepthFromClosed(false);
            ma->WriteIntoFile(argv[argc-1], (char*)"depth");
    }

    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



