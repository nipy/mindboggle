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
    "   -n Neighborhood: set the neighborhood size for computing curvature (default 0.7)\n"
    "   -g GaussianCurvVTK: set the output file to save Gaussian curvature (no default)\n"
    "   -x MaxCurvVTK: set the output file to save maximum curvature (no default)\n"
    "   -i MinCurvVTK: set the output file to save minimum curvature (no default)\n"
    "   -d DirectionVTK: set the output file to save minimal curvature's direction (no default)\n"
    "Example: CurvatureMain -m 0 -n 0.75 -i lh.min_curv.vtk -x lh.max_curv.vtk -g lh.gaussian_curv.vtk -d lh.min_dir.vtk lh.pial.vtk  lh.mean_curv.vtk \n"
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
    bool Gaussian=false, Max=false, Min=false, MinDir=false;  // whether output Gaussian, max and min curvatures, and directions of minimal curvature
    char * GaussianVTK, * MaxVTK, * MinVTK, * MinDirVTK; // files to save Gaussian, max and min curvatures, and directions of minimal curvature
    double NeighborhoodSize = 0.7; // default neighborhood size for computing curvature

    /*Check whether we have mandatory I/Os*/
    int Mandatory = 0; // number of mandatory I/O (input mesh and mean curvature VTK) available 
    for (int i=1;i<argc;i++) // we may need to use getopt or getlongopt later - Forrest, 2012/05/29
    {
        if (argv[i][0] != '-') 
        { 
            Mandatory++;
        }
        else
            i++; // skip the value after a flag
    }
    if (Mandatory < 2)
        {
            cout<<"[ERROR] Not sufficient mandatory I/Os. Check usage please. "<<endl;
            print_help();
            return -4;
        }

    /* Processing options and arguments */
    for (int i=1;i<argc;i++) // we may need to use getopt or getlongopt later - Forrest, 2012/05/29
    {
        if (argv[i][0] != '-') 
            continue; // no more options
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

            case 'i' : // whether output min curvature 
                Min = true;
                MinVTK = argv[i+1];
                break;

            case 'n': //set neighborhood size for computing curvature
                if (Method==0)
                        cout <<"[Warning]: Method 0 does not need neighborhood size though you provided it."<<endl;
                NeighborhoodSize = atof(argv[i+1]);
                break;

            case 'd': // whether output directions of min curvature
                if ((Method==2) or (Method==1))
                {
                    cout<<"[ERROR]: Method 2 or 1 does NOT compute directions of minimal curvature. Remove flag -d and its argument please."<<endl;
                    return -3;
                }
                MinDir = true;
                MinDirVTK = argv[i+1];
                break;
        
            default: 
                cout<<"[ERROR] Unrecognized argument flag or option. Check usage. ";
                print_help();
        } // end of switch
    } // end of looping thru arguments

    /*Compute curvatures*/
    MeshAnalyser* ma = new MeshAnalyser(argv[argc-2]); // the second last input is the inputVTK
    switch (Method)
    {
        case 2:
            ma->ComputeCurvature(NeighborhoodSize);
        case 1: 
            ma->ComputeBothCurvatures(NeighborhoodSize);
        default: // =0  
            vtkDoubleArray* minDirections= ma->ComputePrincipalCurvatures();
            if(MinDir)
            {
                cout<<"Saving directions of minimal curvature into file "<<MinDirVTK<<endl;
                double dir[3];

                ofstream myfile(MinDirVTK);
                myfile.clear();


                for(int i = 0; i<minDirections->GetNumberOfTuples() ; i++)
                {
                    minDirections->GetTuple(i,dir);
                    myfile<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<endl;
                }
                myfile.close();
            }
    }

    /*Write results into VTK files*/
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



