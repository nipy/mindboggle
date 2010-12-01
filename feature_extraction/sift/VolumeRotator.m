function RotatedVolume = VolumeRotator(InputVolume, Degree)


%%%%%%%%%%%%%%%%%%%%%%% For Test Purposes
% Degree = [0 0 10];
% InputVolume = Norm_X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Getting the image size and Rotation Components
[M N R] = size(InputVolume);

if R ==1
    RotatedVolume = imrotate(InputVolume, Degree, 'bicubic', 'crop');
    
    
elseif R>3
    
TetaX = Degree(1);
TetaY = Degree(2);
TetaZ = Degree(3);



%%%%%%%%%%%%%%%%%%%%%%% Start Rotation from x toward z access
for i =1:M
    X_rotated(i,:,:) = imrotate(squeeze(InputVolume(i,:,:)),TetaX,'bicubic','crop');
end
for i =1:N
    XY_rotated(:,i,:) = imrotate(squeeze(X_rotated(:,i,:)),TetaY,'bicubic','crop');
end

for i =1:R
    RotatedVolume(:,:,i) = imrotate(XY_rotated(:,:,i),TetaZ,'bicubic','crop');
end

end

%%%%%%%%%%%%%%%%%%%%%%% For Test purposes 
% figure(1)
% imshow(squeeze(InputVolume(:,:,20)),[min(inputVolume(:)) max(InputVolume(:))])
% figure(2)
% imshow(squeeze(RotatedVolume(:,:,20)), [min(RotatedVolume(:)) max(RotatedVolume(:))])
% 
% figure(1)
% imshow(squeeze(InputVolume(:,50,:)),[min(inputVolume(:)) max(InputVolume(:))])
% figure(2)
% imshow(squeeze(RotatedVolume(:,50,:)), [min(RotatedVolume(:)) max(RotatedVolume(:))])
% 
% figure(1)
% imshow(squeeze(InputVolume(50,:,:)),[min(inputVolume(:)) max(InputVolume(:))])
% figure(2)
% imshow(squeeze(RotatedVolume(50,:,:)), [min(RotatedVolume(:)) max(RotatedVolume(:))])