#include "../MeshAnalyser.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    if(argc < 3)
    {
        cout<<"Usage: CurvatureMain input MeanCurvatureOutput [GaussianCurvatureOutput] [MaximalCurvatureOutput] [MinimalCurvatureOutput]"<<endl;
    }

    MeshAnalyser* ma = new MeshAnalyser(argv[1]);
    vtkDoubleArray* minDirections = ma->ComputePrincipalCurvatures();

    ma->WriteIntoFile(argv[2], (char*)"curv");

    if(argc > 3)
    {
        ma->WriteIntoFile(argv[3], (char*)"gCurv");
    }
    if(argc > 4)
    {
        ma->WriteIntoFile(argv[4], (char*)"curv1");
    }
    if(argc > 5)
    {
        ma->WriteIntoFile(argv[5], (char*)"curv2");
    }

    double dir[3];

    ofstream myfile(argv[6]);
    myfile.clear();


    for(int i = 0; i<minDirections->GetNumberOfTuples() ; i++)
    {
        minDirections->GetTuple(i,dir);
        myfile<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<endl;
    }

    myfile.close();
    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



