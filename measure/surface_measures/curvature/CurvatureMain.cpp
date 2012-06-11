#include "../MeshAnalyser.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

void print_help()
{
    printf(
    "Usage: CurvatureMain [Options] InputVTKMesh MeanCurvatureOutput\n"
    "Options: \n"
    "   -m Method: set the method used to compute curvature(s) (default 0)\n"
    "       0 -- Use ComputePrincipalCurvatures() function to compute both mean and Gaussian curvatures based on the relative direction of the normal vectors in a small neighborhood\n"
    "       1 -- Use ComputeBothCurvatures() function to compute both mean and Gaussian curvatures based on the local ratios betwen a filtered surface and the orginal surface area\n"
    "       2 -- Use ComputeCurvature() function to compute the mean curvature based on the direction of the displacement vectors during a laplacian filtering\n"
    "   -g GaussianCurvVTK: set the output file to save Gaussian curvature (no default)\n"
    "   -x MaxCurvVTK: set the output file to save maximum curvature (no default)\n"
    "   -i MinCurvVTK: set the output file to save minimum curvature (no default)\n"
    "Example: CurvatureMain -i lh.min_curv.vtk -m 0 -x  lh.max_curv.vtk -g lh.gaussian_curv.vtk lh.pial.vtk lh.mean_curv.vtk \n"
        );
}

int main(int argc, char** argv)
{
    cout<<endl;
    time_t start= time(NULL);

    if(argc < 3)
    {
        print_help();
        return -1;
    }
    
    /* Default values for options */
	int Method=0;  // using ComputePrincipalCurvatures()
    bool Gaussian=false, Max=false, Min=false;  // whether output Gaussian, max and min curvatures
    char * GaussianVTK, * MaxVTK, * MinVTK; // files to save Gaussian, max and min curvatures  

    /* Select values for options */
    for (int i=1;i<argc;i++) // we may need to use getopt or getlongopt later - Forrest, 2012/05/29
    {
        if (argv[i][0] != '-') continue; // no more options
        switch (argv[i][1])
        {
            case 'm' :  // select curvature computation method 
                Method = atoi(argv[i+1]);
                cout<<"Using method "<<Method<<" to compute curvature(s)..."<<endl;
                break;

            case 'g' :  // whether output Gaussian curvature 
                if (Method==2) 
                {
                    cout<<"[ERROR]: Method 2 does NOT compute Gaussian curvature."<<endl;
                    return -2;
                }
                Gaussian = true;
                GaussianVTK = argv[i+1];
                break;

            case 'x' : // whether output max curvature 
                Max = true;
                MaxVTK = argv[i+1];
                break;

            case 'i' : // whether output max curvature 
                Min = true;
                MinVTK = argv[i+1];
                break;

        }
    }

    MeshAnalyser* ma = new MeshAnalyser(argv[argc-2]); // the second last input is the inputVTK
    switch (Method)
    {
        case 2:
            ma->ComputeCurvature(0.7); 
        case 1: 
            ma->ComputeBothCurvatures(0.7); 
        default: // =0  
            ma->ComputePrincipalCurvatures();
    }

    ma->WriteIntoFile(argv[argc-1], (char*)"curv");  // the very last input is the VTK for mean curv 

    if(Gaussian)
    {
        cout<<"Saving Gaussian curvature into file "<<GaussianVTK<<endl;
        ma->WriteIntoFile(GaussianVTK, (char*)"gCurv");
    }

    if(Max)
    {
        cout<<"Saving maximum curvature into file "<<MaxVTK<<endl;
        ma->WriteIntoFile(MaxVTK, (char*)"curv1");
    }
    if(Min)
    {
        cout<<"Saving minimum curvature into file "<<MinVTK<<endl;
        ma->WriteIntoFile(MinVTK, (char*)"curv2");
    }

    cout<<"Elapsed time: "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



