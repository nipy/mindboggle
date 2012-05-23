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
    ma->ComputePrincipalCurvatures();

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


    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



