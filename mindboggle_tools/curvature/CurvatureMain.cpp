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
    "       1 -- Use ComputeBothCurvatures() function to compute both mean and Gaussian curvatures based on the local ratios between a filtered surface and the original surface area\n"
    "       2 -- Use ComputeCurvature() function to compute the mean curvature based on the direction of the displacement vectors during a Laplacian filtering\n"
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
<<<<<<< HEAD

=======
    
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    /* Default values for options */
    int Method=0;  // using ComputePrincipalCurvatures()
    bool Gaussian=false, Max=false, Min=false, MinDir=false;  // whether output Gaussian, max and min curvatures, and directions of minimal curvature
    char * GaussianVTK, * MaxVTK, * MinVTK, * MinDirVTK; // files to save Gaussian, max and min curvatures, and directions of minimal curvature
    double NeighborhoodSize = 0.7; // default neighborhood size for computing curvature

    /*Check whether we have mandatory I/Os*/
<<<<<<< HEAD
    int Mandatory = 0; // number of mandatory I/O (input mesh and mean curvature VTK) available
    for (int i=1;i<argc;i++) // we may need to use getopt or getlongopt later - Forrest, 2012/05/29
    {
        if (argv[i][0] != '-')
        {
=======
    int Mandatory = 0; // number of mandatory I/O (input mesh and mean curvature VTK) available 
    for (int i=1;i<argc;i++) // we may need to use getopt or getlongopt later - Forrest, 2012/05/29
    {
        if (argv[i][0] != '-') 
        { 
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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
<<<<<<< HEAD
        if (argv[i][0] != '-')
            continue; // no more options
        switch (argv[i][1])
        {
            case 'm' :  // select curvature computation method
=======
        if (argv[i][0] != '-') 
            continue; // no more options
        switch (argv[i][1])
        {
            case 'm' :  // select curvature computation method 
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
                Method = atoi(argv[i+1]);
                cout<<"Using method "<<Method<<" to compute curvature(s)..."<<endl;
                break;

<<<<<<< HEAD
            case 'g' :  // whether output Gaussian curvature
                if (Method==2)
=======
            case 'g' :  // whether output Gaussian curvature 
                if (Method==2) 
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
                {
                    cout<<"[ERROR]: Method 2 does NOT compute Gaussian curvature."<<endl;
                    return -2;
                }
                Gaussian = true;
                GaussianVTK = argv[i+1];
                break;

<<<<<<< HEAD
            case 'x' : // whether output max curvature
=======
            case 'x' : // whether output max curvature 
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
                Max = true;
                MaxVTK = argv[i+1];
                break;

<<<<<<< HEAD
            case 'i' : // whether output min curvature
=======
            case 'i' : // whether output min curvature 
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
                Min = true;
                MinVTK = argv[i+1];
                break;

            case 'n': //set neighborhood size for computing curvature
<<<<<<< HEAD
//                if (Method==0)
//                        cout <<"[Warning]: Method 0 does not need neighborhood size though you provided it."<<endl;
=======
                if (Method==0)
                        cout <<"[Warning]: Method 0 does not need neighborhood size though you provided it."<<endl;
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
                NeighborhoodSize = atof(argv[i+1]);
                break;

            case 'd': // whether output directions of min curvature
<<<<<<< HEAD
                if ((Method==2) || (Method==1))
=======
                if ((Method==2) or (Method==1))
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
                {
                    cout<<"[ERROR]: Method 2 or 1 does NOT compute directions of minimal curvature. Remove flag -d and its argument please."<<endl;
                    return -3;
                }
                MinDir = true;
                MinDirVTK = argv[i+1];
                break;
<<<<<<< HEAD

            default:
=======
        
            default: 
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
                cout<<"[ERROR] Unrecognized argument flag or option. Check usage. ";
                print_help();
        } // end of switch
    } // end of looping thru arguments

    /*Compute curvatures*/
    MeshAnalyser* ma = new MeshAnalyser(argv[argc-2]); // the second last input is the inputVTK
    switch (Method)
    {
        case 2:
<<<<<<< HEAD
            ma->ComputeCurvature(NeighborhoodSize, 20);
            break;
        case 1:
            ma->ComputeBothCurvatures(NeighborhoodSize);
            break;
        default: // =0
            vtkDoubleArray* minDirections= ma->ComputePrincipalCurvatures(NeighborhoodSize);
=======
            ma->ComputeCurvature(NeighborhoodSize);
        case 1: 
            ma->ComputeBothCurvatures(NeighborhoodSize);
        default: // =0  
            vtkDoubleArray* minDirections= ma->ComputePrincipalCurvatures();
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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
<<<<<<< HEAD
            break;
    }
    /*Write results into VTK files*/
    ma->WriteIntoFile(argv[argc-1], (char*)"curv");  // the very last input is the VTK for mean curv
=======
    }

    /*Write results into VTK files*/
    ma->WriteIntoFile(argv[argc-1], (char*)"curv");  // the very last input is the VTK for mean curv 

>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    if(Gaussian)
    {
        cout<<"Saving Gaussian curvature into file "<<GaussianVTK<<endl;
        ma->WriteIntoFile(GaussianVTK, (char*)"gCurv");
    }
<<<<<<< HEAD
=======

>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    if(Max)
    {
        cout<<"Saving maximum curvature into file "<<MaxVTK<<endl;
        ma->WriteIntoFile(MaxVTK, (char*)"curv1");
    }
<<<<<<< HEAD

=======
   
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    if(Min)
    {
        cout<<"Saving minimum curvature into file "<<MinVTK<<endl;
        ma->WriteIntoFile(MinVTK, (char*)"curv2");
    }
<<<<<<< HEAD
=======

>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    cout<<"Elapsed time: "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



