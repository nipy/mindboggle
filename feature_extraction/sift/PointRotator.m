function RotatedPoint = PointRotator(InputPoint, Degree, RotationCenter)
 
[M N R] = size(InputPoint);



if N==2
    RotationCenter = [(RotationCenter(2)) (RotationCenter(1))];

    for i = 1:M
        TransformedPoint = InputPoint(i,1:N) - RotationCenter;
        RotationMatrix = [cosd(Degree) -sind(Degree); sind(Degree) cosd(Degree)];
        RotatedTransformedPoint = RotationMatrix * TransformedPoint';
        RotatedPoint(i,1:N) = RotatedTransformedPoint' + RotationCenter;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%% This is coresponds to VolumeRotator Order where x
%%%%%%%%%%%%%%%%%%%%%%%%% to y to z order is conserved.
if N==3

TetaX = Degree(1);
TetaY = Degree(2);
TetaZ = Degree(3);
RotatedPoint = InputPoint;

% Rx = [ 1               0              0
%        0               cosd(TetaX)   -sind(TetaX)
%        0               sind(TetaX)    cosd(TetaX) ];
%       
% Ry = [ cosd(TetaY)     0              sind(TetaY)
%        0               1              0
%       -sind(TetaY)     0              cosd(TetaY) ];
%  
%         
% Rz = [ cosd(TetaZ)    -sind(TetaZ)    0
%        sind(TetaZ)     cosd(TetaZ)    0
%        0               0              1           ];
% 
%     RotationMatrix = Rx*Ry*Rz;
%     %RotationCenter = [(RotationCenter(1)) (RotationCenter(3)) (RotationCenter(2))];
% 
%     for i = 1:M
%         TransformedPoint = InputPoint(i,1:N) - RotationCenter;
%         RotatedTransformedPoint = RotationMatrix * TransformedPoint';
%         RotatedPoint(i,1:N) = RotatedTransformedPoint' + RotationCenter;
%     end
% end 

    if (TetaZ)
%rotation along z axes
    RotationCenterZ = [(RotationCenter(1)) (RotationCenter(2))];
    Rz = [ cosd(TetaZ)    -sind(TetaZ)  
           sind(TetaZ)     cosd(TetaZ)];  
       
    for i = 1:M
        TransformedPoint = InputPoint(i,[1,2]) - RotationCenterZ;
        RotatedTransformedPoint = Rz * TransformedPoint';
        RotatedPoint(i,[1,2]) = RotatedTransformedPoint' + RotationCenterZ;
    end
    elseif(TetaY) 
%rotation along y axes
    RotationCenterY = [(RotationCenter(1)) (RotationCenter(3))];
    Ry = [ cosd(TetaY)    -sind(TetaY)  
           sind(TetaY)     cosd(TetaY)];  
       
    for i = 1:M
        TransformedPoint = RotatedPoint(i,[1,3]) - RotationCenterY;
        RotatedTransformedPoint = Ry * TransformedPoint';
        RotatedPoint(i,[1,3]) = RotatedTransformedPoint' + RotationCenterY;
    end
    elseif(TetaX)
%rotation along x axes
    RotationCenterX = [(RotationCenter(2)) (RotationCenter(3))];
    Rx = [ cosd(TetaX)    -sind(TetaX)  
           sind(TetaX)     cosd(TetaX)];  
       
    for i = 1:M
        TransformedPoint = RotatedPoint(i,[2,3]) - RotationCenterX;
        RotatedTransformedPoint = Rx * TransformedPoint';
        RotatedPoint(i,[2,3]) = RotatedTransformedPoint' + RotationCenterX;
    end
    end
    
end  
%   
% 
% [M N R] = size(Norm_X);
% PointsImage = zeros(M,N,R);
% [m n r] = size(InputPoint);
% for i = 1:m
%     PointsImage(round(InputPoint(i,1)), round(InputPoint(i,2)), round(InputPoint(i,3))) = 1;
% end
% 
% NewVol = VolumeRotator(PointsImage, -Degree);
% 
% 
% 





