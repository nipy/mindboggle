function [Adjustment NewMaxima Edge]= TaylorMaxima (InputVolume)

[M N R P] = size(InputVolume);

%%%%%%%%%%%%%%%%%%%%%% To have a accurate Hessian matrix value The
%%%%%%%%%%%%%%%%%%%%%% input volume size needs to be at least one extra
%%%%%%%%%%%%%%%%%%%%%% from each corner, i. e. for a 3x3x3 filter we
%%%%%%%%%%%%%%%%%%%%%% need the input volume to be in the size of 5x5x5

if (P ==1 && R~=1 && N~=1 && M~=1) 
    
    if (M~=5 || N~=5 || R~=3)
        fprintf('the input matrix is not the right size');
        return
    end
    EnlargedVolume = zeros(5,5,5);
    EnlargedVolume(:,:,2:4) = InputVolume;
    EnlargedVolume(:,:,1) = InputVolume(:,:,1);
    EnlargedVolume(:,:,5) = InputVolume(:,:,3);

    [D_Y D_X D_Z] =gradient(EnlargedVolume);
    JacobianVector = [ D_X(3,3,3)  D_Y(3,3,3)  D_Z(3,3,3) ];
    
    %%%%%%%%%%%%%%%%%%%%%% Compute Hessian Matrix
    [D_XY D_XX D_XZ] = gradient(D_X);
    [D_YY D_YX D_YZ] = gradient(D_Y);
    [D_ZY D_ZX D_ZZ] = gradient(D_Z);
  
  
    HessianMatrix = [ D_XX(3,3,3)  D_XY(3,3,3)  D_XZ(3,3,3)
                      D_YX(3,3,3)  D_YY(3,3,3)  D_YZ(3,3,3)
                      D_ZX(3,3,3)  D_ZY(3,3,3)  D_ZZ(3,3,3) ];


    Trace_H = trace(HessianMatrix(1:2,1:2));
    Deter_H = det(HessianMatrix(1:2,1:2));
    
    if (Trace_H^2/Deter_H) > 12.1
        Edge = 0;
    else
        Edge = 1;
    end
    
    Adjustment = -(HessianMatrix^(-1))*JacobianVector';
    
    if norm(Adjustment) > .86
        Adjustment = [0 0 0]';
    end
    
    NewMaxima = EnlargedVolume(3,3,3) + .5*JacobianVector*Adjustment;



elseif (P >1 && R~=1 && N~=1 && M~=1) 
    
    
    if (M~=5 || N~=5 || R~=5 || P~=3)
        fprintf('the input matrix is not the right size');
        return
    end
    EnlargedVolume = zeros(5,5,5,5);
    EnlargedVolume(:,:,:,2:4) = InputVolume;
    EnlargedVolume(:,:,:,1) = InputVolume(:,:,:,1);
    EnlargedVolume(:,:,:,5) = InputVolume(:,:,:,3);
    
    %%%%%%%%%%%%%%%%%%%%%% Compute Jacobinan Vector
    [D_Y D_X D_Z D_W] =gradient(EnlargedVolume);
    JacobianVector = [ D_X(3,3,3,3)  D_Y(3,3,3,3)  D_Z(3,3,3,3) D_W(3,3,3,3) ];
    
    %%%%%%%%%%%%%%%%%%%%%% Compute Hessian Matrix
    [D_XY D_XX D_XZ D_XW] = gradient(D_X);
    [D_YY D_YX D_YZ D_YW] = gradient(D_Y);
    [D_ZY D_ZX D_ZZ D_ZW] = gradient(D_Z);
    [D_WY D_WX D_WZ D_WW] = gradient(D_W);
  
  
    HessianMatrix = [ D_XX(3,3,3,3)  D_XY(3,3,3,3)  D_XZ(3,3,3,3)  D_XW(3,3,3,3)
                      D_YX(3,3,3,3)  D_YY(3,3,3,3)  D_YZ(3,3,3,3)  D_YW(3,3,3,3)
                      D_ZX(3,3,3,3)  D_ZY(3,3,3,3)  D_ZZ(3,3,3,3)  D_ZW(3,3,3,3)
                      D_WX(3,3,3,3)  D_WY(3,3,3,3)  D_WZ(3,3,3,3)  D_WW(3,3,3,3)];

    Trace_H = trace(HessianMatrix(1:3,1:3));
    Deter_H = det(HessianMatrix(1:3,1:3));
    
    if (Trace_H^2/Deter_H) > 12.1
        Edge = 0;
    else
        Edge = 1;
    end
    
    Adjustment = -(HessianMatrix^(-1))*JacobianVector';
  
    if norm(Adjustment) > .86
        Adjustment = [0 0 0 0]';
    end

    NewMaxima = EnlargedVolume(3,3,3) + .5*JacobianVector*Adjustment;

end
