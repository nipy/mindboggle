function [Adjustment NewMaxima Edge]= TaylorExtremum (InputVolume)

[M N R P] = size(InputVolume);
Edge = 0;
Adjustment = [0 0 0];

%%%%%%%%%%%%%%%%%%%%%% To have a accurate Hessian matrix value The
%%%%%%%%%%%%%%%%%%%%%% input volume size needs to be at least one extra
%%%%%%%%%%%%%%%%%%%%%% from each corner, i. e. for a 3x3x3 filter we
%%%%%%%%%%%%%%%%%%%%%% need the input volume to be in the size of 5x5x5

if (P ==1 && R~=1 && N~=1 && M~=1) 
    
    if (M~=3 || N~=3 || R~=3)
        fprintf('\nthe input matrix is not the right size\n');
        return
    end

    D_X = (InputVolume(3,2,2) - InputVolume(1,2,2))/2.0;
    D_Y = (InputVolume(2,3,2) - InputVolume(2,1,2))/2.0;
    D_Z = (InputVolume(2,2,3) - InputVolume(2,2,1))/2.0;

    JacobianVector = [ D_X  D_Y  D_Z ];
    
    %%%%%%%%%%%%%%%%%%%%%% Compute Hessian Matrix
    D_XX = InputVolume(3,2,2) - 2*InputVolume(2,2,2) + InputVolume(1,2,2);
    D_XY = 0.25*(InputVolume(3,3,2) + InputVolume(1,1,2) - InputVolume(3,1,2) - InputVolume(1,3,2));
    D_XZ = 0.25*(InputVolume(3,2,3) + InputVolume(1,2,1) - InputVolume(3,2,1) - InputVolume(1,2,3));
    
    D_YX = D_XY;
    D_YY = InputVolume(2,3,2) - 2*InputVolume(2,2,2) + InputVolume(2,1,2);
    D_YZ = 0.25*(InputVolume(2,3,3) + InputVolume(2,1,1) - InputVolume(2,3,1) - InputVolume(2,1,3));
  
    D_ZX = D_XZ;
    D_ZY = D_YZ;
    D_ZZ = InputVolume(2,2,3) - 2*InputVolume(2,2,2) + InputVolume(2,2,1);

    HessianMatrix = [ D_XX  D_XY  D_XZ
                      D_YX  D_YY  D_YZ
                      D_ZX  D_ZY  D_ZZ ];

     
    Adjustment = -(HessianMatrix^(-1))*JacobianVector';
   
    NewMaxima = InputVolume(2,2,2) + .5*JacobianVector*Adjustment;


    %%%%%%%%%%%%%%%%%%%%%% The Edge Detector
    Trace_H = trace(HessianMatrix(1:2,1:2));
    Deter_H = det(HessianMatrix(1:2,1:2));
    Score = Trace_H^2/Deter_H;
    
    if (Score<12.1 && Score>0)
        Edge = 0;
    else
        Edge = 1;
    end
    
elseif (P >1 && R~=1 && N~=1 && M~=1) 
    
    
    if (M~=3 || N~=3 || R~=3 || P~=3)
        fprintf('the input matrix is not the right size');
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Compute Jacobinan Vector
    D_X = (InputVolume(3,2,2,2) - InputVolume(1,2,2,2))/2.0;
    D_Y = (InputVolume(2,3,2,2) - InputVolume(2,1,2,2))/2.0;
    D_Z = (InputVolume(2,2,3,2) - InputVolume(2,2,1,2))/2.0;
    D_W = (InputVolume(2,2,2,3) - InputVolume(2,2,2,1))/2.0;

    JacobianVector = [ D_X  D_Y  D_Z  D_W];
    
   
    %%%%%%%%%%%%%%%%%%%%%% Compute Hessian Matrix
    D_XX = InputVolume(3,2,2,2) - 2*InputVolume(2,2,2,2) + InputVolume(1,2,2,2);
    D_XY = 0.25*(InputVolume(3,3,2,2) + InputVolume(1,1,2,2) - InputVolume(3,1,2,2) - InputVolume(1,3,2,2));
    D_XZ = 0.25*(InputVolume(3,2,3,2) + InputVolume(1,2,1,2) - InputVolume(3,2,1,2) - InputVolume(1,2,3,2));
    D_XW = 0.25*(InputVolume(3,2,2,3) + InputVolume(1,2,2,1) - InputVolume(3,2,2,1) - InputVolume(1,2,2,3));
    
    D_YX = D_XY;
    D_YY = InputVolume(2,3,2,2) - 2*InputVolume(2,2,2,2) + InputVolume(2,1,2,2);
    D_YZ = 0.25*(InputVolume(2,3,3,2) + InputVolume(2,1,1,2) - InputVolume(2,3,1,2) - InputVolume(2,1,3,2));
    D_YW = 0.25*(InputVolume(2,3,2,3) + InputVolume(2,1,2,1) - InputVolume(2,3,2,1) - InputVolume(2,1,2,3));
    
    D_ZX = D_XZ;
    D_ZY = D_YZ;
    D_ZZ = InputVolume(2,2,3,2) - 2*InputVolume(2,2,2,2) + InputVolume(2,2,1,2);
    D_ZW = 0.25*(InputVolume(2,2,3,3) + InputVolume(2,2,1,1) - InputVolume(2,2,3,1) - InputVolume(2,2,1,3));
     
    D_WX = D_XW;
    D_WY = D_YW;
    D_WZ = D_ZW;
    D_WW = InputVolume(2,2,2,3) - 2*InputVolume(2,2,2,2) + InputVolume(2,2,2,1);

    HessianMatrix = [ D_XX  D_XY  D_XZ  D_XW
                      D_YX  D_YY  D_YZ  D_YW
                      D_ZX  D_ZY  D_ZZ  D_ZW
                      D_WX  D_WY  D_WZ  D_WW  ];

    
    Adjustment = -(HessianMatrix^(-1))*JacobianVector';
   
    NewMaxima = InputVolume(2,2,2,2) + .5*JacobianVector*Adjustment;

 end
