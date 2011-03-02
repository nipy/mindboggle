clear
clc

addpath /data/export/home/ray/Images/MNI_Transfered/
addpath /data/export/home/ray/NIFTI_20100106/
addpath /data/export/home/ray/SIFT/mFiles/
addpath /data/export/home/ray/Images
addpath /data/export/home/ray/SURF
addpath /data/export/home/ray/SIFT/Vedaldi


% image1 = load_nii('/data/export/home/ray/Images/ColumbiaScans/brain1_1.nii');
% X = image1.img;
% Slice_X = X(:, :, 70);

X = imread('/data/export/home/ray/Images/NonMedicals/Mollas.jpg'); 

%Slice_X = rgb2gray(X);
Slice_X = X;
        
[M N R] = size(Slice_X);
BackGround = zeros(M+N);
fprintf(['\n Volume size is: ' int2str(M) 'x' int2str(N) 'x' int2str(R) '\n']);
        
NormalizationType = 3;
NgL=256;
Norm_X = ImageNormalizer(Slice_X, NormalizationType, NgL);
StartPoint1 = round(N/2);
StartPoint2 = round(M/2);

BackGround((StartPoint1+1):(M+StartPoint1),(StartPoint2+1):(N+StartPoint2)) = Norm_X;
[M N R] = size(BackGround);

clear image1
clear Slice_X


  SamplePerScale = 3 ; InitialSigma = 2.0159; Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5; Ther = 0.03; SearchArea=3; EdgeDetect=1;
%  %%%%%%%%%%%%%%%%%%%% Keypoints of the Original Image 
  [KeyPoints_SIFT  Ignord_KeyPoints] = SIFT_KeyPoints(BackGround, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
        
        
        
% Set this option to true if you want to see more information
  Options.verbose=false; 
% Get the Key Points
  Ipts=OpenSurf(BackGround,Options);
  SURF_KeyPoints = [Ipts.x ;Ipts.y ;Ipts.scale]';
  
  %%%% Vedaldi Keypoints 
  [frames,descriptors,gss,dogss] = sift(BackGround,'verbosity', 0);
  Vedaldi_Keypoints = (frames(1:3,:))';

%   
%   for i = 3:3
   Rotation = 90;
  Rotated_X = imrotate(BackGround, Rotation, 'bicubic', 'crop');
%   NormalizationType = 3;
%   NgL=256;
%   Rotated_X = ImageNormalizer(Rotated_X, NormalizationType, NgL);
% 
%   
%   Ipts=OpenSurf(Rotated_X,Options);
%   SURF_Rotated_KeyPoints = [Ipts.x ;Ipts.y ;Ipts.scale]';
%   Rotated_Back_Point = PointRotator(SURF_Rotated_KeyPoints(:, 1:2), -Rotation, [M/2 N/2]);
%   
   [Rotated_KeyPoints_SIFT Ignord_KeyPoints] = SIFT_KeyPoints(Rotated_X, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
   Rotated_Back_SIFT_Point = PointRotator(Rotated_KeyPoints_SIFT(:, 1:2), Rotation, [(M+1)/2 (N+1)/2]);
% 
%   [frames,descriptors,gss,dogss] = sift(Rotated_X,'verbosity', 0);
%   Rotated_Vedaldi_Keypoints = (frames(1:3,:))';
%   Rotated_Back_Vedaldi_Point = PointRotator(Rotated_Vedaldi_Keypoints(:, 1:2), -Rotation, [M/2 N/2]);
% 
%  InitialSigma = 2;
%  [MatchPersentage(i+1) MatchingPoints] = KeyPointsMatching(SURF_KeyPoints(:,1:3), [Rotated_Back_Point [Ipts.scale]'] , InitialSigma);
%  NumberOfSURFKeyPoints(i+1) = length(MatchingPoints);
%  
  [SIFT_MatchPersentage(1) MatchingPoints] = KeyPointsMatching(KeyPoints_SIFT(:,1:3), [Rotated_Back_SIFT_Point Rotated_KeyPoints_SIFT(:,3)] , InitialSigma);
%   NumberOfSIFTKeyPoints(i+1) = length(MatchingPoints);
% 
%    [Vedaldi_MatchPersentage(i+1) MatchingPoints] = KeyPointsMatching( Vedaldi_Keypoints(:,1:3), [Rotated_Back_Vedaldi_Point Rotated_Vedaldi_Keypoints(:,3)] , InitialSigma);
%   NumberOfVedaldiKeyPoints(i+1) = length(MatchingPoints);
% 
%   end
% 
%   x = 0:30:360;
%   plot(x,MatchPersentage,'-r+')
%   hold on
%   plot(x,SIFT_MatchPersentage, '--bo')
%   hold on
%   plot(x,Vedaldi_MatchPersentage, '-.g*')
%   legend ('SURF', 'SIFT', 'Vedaldi')
%   xlabel ('Degree')
%   ylabel('Keypoint Matching Percentage')
%   grid on
%   
%   figure(2)
%   plot(x,NumberOfSURFKeyPoints, '-r+')
%   hold on
%   plot(x,NumberOfSIFTKeyPoints, '--bo')
%   hold on
%   plot(x,NumberOfVedaldiKeyPoints,'-.g*')
%   legend ('SURF', 'SIFT', 'Vedaldi')
%   grid on
%   xlabel('Degree')
%   ylabel('Number of Keypoints')
%   
   figure(3)
   imshow((imrotate((BackGround),0)),[min(BackGround(:)) max(BackGround(:))]);
   hold on
   plot(KeyPoints_SIFT(:,1), KeyPoints_SIFT(:,2), 'bo', 'MarkerSize', 5)
   hold on
   plot(Rotated_Back_SIFT_Point(:,1), Rotated_Back_SIFT_Point(:,2), 'r+', 'MarkerSize', 5)

   
   figure(4)
   imshow((imrotate((BackGround),0)),[min(BackGround(:)) max(BackGround(:))]);
   hold on
   plot(Rotated_KeyPoints_SIFT(:,1), Rotated_KeyPoints_SIFT(:,2), 'bo', 'MarkerSize', 5)
   hold on
   plot(Rotated_Back_SIFT_Point(:,1), Rotated_Back_SIFT_Point(:,2), 'r+', 'MarkerSize', 5)
   hold on
   
   figure(4)
   imshow((imrotate((BackGround),0)),[min(BackGround(:)) max(BackGround(:))]);
   hold on
   plot(Vedaldi_Keypoints(:,1), Vedaldi_Keypoints(:,2), 'r+', 'MarkerSize', 5)
   hold on
   plot(Rotated_Back_Vedaldi_Point(:,1), Rotated_Back_Vedaldi_Point(:,2), 'co', 'MarkerSize', 5)
   hold on 
   
   plot(SURF_KeyPoints(:,1), SURF_KeyPoints(:,2), 'c*', 'MarkerSize', 5)
   
   
%    
% x = (1:10)';
% y = (ones(1,10))'
% Square = [10*y, x; x, 10*y; x, y; y, x]
% RotatedShape = PointRotator(Square, 90, [5.5 5.5])
% plot(Square(:,1),Square(:,2))
% hold on
% plot(RotatedShape(:,1),RotatedShape(:,2), 'r')
% xlim([-1 15]) 
% ylim([-1 15]) 
% grid on
% 


I = round(200*rand(11,11));
Rotated_I = imrotate(I,90,'bicubic', 'crop');
figure(1); mesh(double(I));
figure(2); mesh(double(Rotated_I));

