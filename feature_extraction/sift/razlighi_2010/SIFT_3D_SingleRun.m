
clear all
clc

addpath /data/export/home/ray/Images/ColumbiaScans/
addpath /data/export/home/ray/NIFTI_20100106

addpath /Users/ray/Dropbox/Research/SIFT/mFiles/
addpath /Users/ray/Dropbox/Research/Images/NonMedicals/
addpath /Users/ray/Dropbox/Research/Images/
addpath /Users/ray/Dropbox/Research/SIFT/Vedaldi
addpath /Users/ray/Dropbox/Research/NIFTI_20100106

        
        image1 = load_untouch_nii('brain2_1_Reg2_NMI152.nii');
        X = image1.img;
        
        %X = X(64:192, 31:92, 64:192);
        
        [M N R] = size(X);
        fprintf(['\n Volume size is: ' int2str(M) 'x' int2str(N) 'x' int2str(R) '\n']);
        
        NormalizationType = 3;
        NgL=256;

        Norm_X = ImageNormalizer(X, NormalizationType, NgL);
        [M N R] = size(Norm_X);
        clear image1
        clear X

 
  SamplePerScale = 3 ; InitialSigma = 2.015874; Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5; Ther = 0.001; SearchArea=3; EdgeDetect = false;
%  %%%%%%%%%%%%%%%%%%%% Keypoints of the Original Image 
[KeyPoints_SIFT1  Ignord_KeyPoints1] = SIFT_KeyPoints_V2(Norm_X, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);


Rotation = [45 0 0];
Norm_X2 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
[KeyPoints_SIFT2  Ignord_KeyPoints2] = SIFT_KeyPoints_V2(Norm_X2, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
RotatedPoint = PointRotator(KeyPoints_SIFT2(:,1:3), [-Rotation(1,1) Rotation(1,2) -Rotation(1,3)], [(M)/2 (N-1)/2 (R-1)/2]);
RotatedPoint45 = [RotatedPoint KeyPoints_SIFT2(:,4)];
[MatchPersentage45 MatchingPoints45] = KeyPointsMatching(KeyPoints_SIFT1, RotatedPoint45, 2);


Rotation = [-45 0 0];
Norm_X3 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
[KeyPoints_SIFT3  Ignord_KeyPoints3] = SIFT_KeyPoints_V2(Norm_X3, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
RotatedPoint = PointRotator(KeyPoints_SIFT3(:,1:3), [-Rotation(1,1) -Rotation(1,2) -Rotation(1,3)], [(M)/2 (N-1)/2 (R-1)/2]);
RotatedPoint_45 = [RotatedPoint KeyPoints_SIFT3(:,4)];
[MatchPersentage_45 MatchingPoints_45] = KeyPointsMatching(KeyPoints_SIFT1, RotatedPoint_45, 2);

[MatchPersentage1 MatchingPoints1] = KeyPointsMatching(MatchingPoints45, MatchingPoints_45, 2);



        image1 = load_untouch_nii('brain2_2_Reg2_NMI152.nii');
        X = image1.img;
        
        %X = X(64:192, 31:92, 64:192);
        
        [M N R] = size(X);
        fprintf(['\n Volume size is: ' int2str(M) 'x' int2str(N) 'x' int2str(R) '\n']);
        
        NormalizationType = 3;
        NgL=256;

        Norm_X = ImageNormalizer(X, NormalizationType, NgL);
        [M N R] = size(Norm_X);
        clear image1
        clear X

 
  SamplePerScale = 3 ; InitialSigma = 2.015874; Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5; Ther = 0.001; SearchArea=3; EdgeDetect = false;
%  %%%%%%%%%%%%%%%%%%%% Keypoints of the Original Image 
[KeyPoints_SIFT1  Ignord_KeyPoints1] = SIFT_KeyPoints_V2(Norm_X, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);


Rotation = [45 0 0];
Norm_X2 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
[KeyPoints_SIFT2  Ignord_KeyPoints2] = SIFT_KeyPoints_V2(Norm_X2, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
RotatedPoint = PointRotator(KeyPoints_SIFT2(:,1:3), [-Rotation(1,1) Rotation(1,2) -Rotation(1,3)], [(M)/2 (N-1)/2 (R-1)/2]);
RotatedPoint45 = [RotatedPoint KeyPoints_SIFT2(:,4)];
[MatchPersentage45 MatchingPoints45] = KeyPointsMatching(KeyPoints_SIFT1, RotatedPoint45, 2);


Rotation = [-45 0 0];
Norm_X3 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
[KeyPoints_SIFT3  Ignord_KeyPoints3] = SIFT_KeyPoints_V2(Norm_X3, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
RotatedPoint = PointRotator(KeyPoints_SIFT3(:,1:3), [-Rotation(1,1) -Rotation(1,2) -Rotation(1,3)], [(M)/2 (N-1)/2 (R-1)/2]);
RotatedPoint_45 = [RotatedPoint KeyPoints_SIFT3(:,4)];
[MatchPersentage_45 MatchingPoints_45] = KeyPointsMatching(KeyPoints_SIFT1, RotatedPoint_45, 2);

[MatchPersentage2 MatchingPoints2] = KeyPointsMatching(MatchingPoints45, MatchingPoints_45, 2);


[MatchPersentage MatchingPoints] = KeyPointsMatching(MatchingPoints1, MatchingPoints2, 2);







% figure(11)
% scatter3(KeyPoints_SIFT1(: ,1), KeyPoints_SIFT1(: ,2), KeyPoints_SIFT1(: ,3), '+b')
% hold on
% scatter3(RotatedPoint20(: ,1), RotatedPoint20(: ,2), RotatedPoint20(: ,3), '+r')

ShowSliceKeyPoints(Norm_X, MatchingPoints45, MatchingPoints_45, 3, 90)

%  %%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by 20 degree
%  Rotation = (round(20*rand(1,3))) - 10;
%  Rotation = [20 20 20];
%  Norm_X2 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
%  KeyPoints2 = SIFT_KeyPoints(Norm_X2, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
%  RotatedPoint = PointRotator(KeyPoints2(:,1:3), [-Rotation(1,1) Rotation(1,2) -Rotation(1,3)], [M/2 N/2 R/2]);
%  RotatedPoint20 = [RotatedPoint KeyPoints2(:,4)];
% 
%  %%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by -20 degree
%  Rotation = (round(20*rand(1,3))) - 10;
%  Rotation = [-20 -20 -20];
%  Norm_X3 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
%  KeyPoints3 = SIFT_KeyPoints(Norm_X3, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
%  RotatedPoint = PointRotator(KeyPoints3(:,1:3), [-Rotation(1,1) Rotation(1,2) -Rotation(1,3)], [M/2 N/2 R/2]);
%  RotatedPoint_20 = [RotatedPoint KeyPoints3(:,4)];
% 
%  
%  
%  InitialSigma = 2;
%  
% %%%%%%%%%%%%%%%%%%%% Computation of the Matches  
% [MatchPersentage1 MatchingPoints1] = KeyPointsMatching(KeyPoints1, RotatedPoint20, InitialSigma);
% 
% [MatchPersentage2 MatchingPoints2] = KeyPointsMatching(KeyPoints1, RotatedPoint_20, InitialSigma);
%  
% [MatchPersentage3 MatchingPoints3] = KeyPointsMatching(MatchingPoints1, MatchingPoints2, InitialSigma);
%  
%  
% [MatchPersentage4 MatchingPoints4 NearestDistance] = KeyPointsMatching(MatchingPoints3, RotatedPoint5, InitialSigma);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%% Visualization of the matches
% % figure(1)
% 
% imshow((imrotate(fliplr(Norm_X(:,:,120)),90)),[min(Norm_X(:)) max(Norm_X(:))]);
% hold on
% 
% scatter3(KeyPoints1(:,1), KeyPoints1(:,2), KeyPoints1(:,3), 'o', 'b')
% hold on
% scatter3(RotatedPoint_10(:,1), RotatedPoint_10(:,2), RotatedPoint_10(:,3), '*', 'r')
% hold on
% 
% 
% 
% 
% 
% subimage(1,3,1); ShowSliceKeyPoints(Norm_X, KeyPoints1, RotatedPoint_10, 1, 60);
% subimage(1,3,2); ShowSliceKeyPoints(Norm_X, KeyPoints1, RotatedPoint_10, 2, 30);
% subimage(1,3,3); ShowSliceKeyPoints(Norm_X, KeyPoints1, RotatedPoint_10, 3, 60);
% 

% 
% 
% % 
% % [frames,descriptors,gss,dogss]=sift(Norm_X2,'Verbosity', 1);
% % Vedaldi = frames(1:2, :)';
% % RotatedPoint = PointRotator(Vedaldi, -90, [128 128]);
% % FlipedRotatedPoint = PointFliper(RotatedPoint, 2, [128 128]);
% % scatter(FlipedRotatedPoint(:,1),FlipedRotatedPoint(:,2),'.','g')
% %  
% 
% %%%%%%%%%%%%%%%%%%%% Keypoints of the Original Image 
% [frames,descriptors,gss,dogss]=sift(Norm_X,'Verbosity', 1);
% KeyPoints1 = frames(1:3, :)';
% 
% %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by 10 degree
%  Norm_X2 = imrotate(Norm_X,10,'bilinear','crop');
%  [frames,descriptors,gss,dogss]=sift(Norm_X2,'Verbosity', 1);
% 
%  KeyPoints2 = frames(1:3, :)';
%  RotatedPoint = PointRotator(KeyPoints2(:,1:2), -10, [128 128]);
%  RotatedPoint10 = [RotatedPoint KeyPoints2(:,3)];
%  
%  
%  %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by -10 degree
%  Norm_X3 = imrotate(Norm_X,-10,'bilinear','crop');
%  [frames,descriptors,gss,dogss]=sift(Norm_X3,'Verbosity', 1);
% 
%  KeyPoints3 = frames(1:3, :)';
%  RotatedPoint = PointRotator(KeyPoints3(:,1:2), 10, [128 128]);
%  RotatedPoint_10 = [RotatedPoint KeyPoints3(:,3)];
% 
%  
%  %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by 5 degree
%  Norm_X2 = imrotate(Norm_X,20,'bilinear','crop');
%  [frames,descriptors,gss,dogss]=sift(Norm_X2,'Verbosity', 1);
%  KeyPoints2 = frames(1:3, :)';
%  RotatedPoint = PointRotator(KeyPoints2(:,1:2), -20, [128 128]);
%  RotatedPoint5 = [RotatedPoint KeyPoints2(:,3)];
%  
% 
%  InitialSigma = 5;
%  
% %%%%%%%%%%%%%%%%%%%%%% Computation of the Matches  
% [MatchPersentage1 MatchingPoints1] = KeyPointsMatching(KeyPoints1, RotatedPoint10, InitialSigma);
% 
% [MatchPersentage2 MatchingPoints2] = KeyPointsMatching(KeyPoints1, RotatedPoint_10, InitialSigma);
%  
% [MatchPersentage3 MatchingPoints3] = KeyPointsMatching(MatchingPoints1, MatchingPoints2, InitialSigma);
%  
%  
% [MatchPersentage4 MatchingPoints4 NearestDistance] = KeyPointsMatching(MatchingPoints3, RotatedPoint5, InitialSigma);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%% Visualization of the matches
% figure(1)
% imshow((imrotate(fliplr(Norm_X),90)),[min(Norm_X(:)) max(Norm_X(:))]);
% hold on
% 
% scatter(RotatedPoint5(:,1),RotatedPoint5(:,2),'o','b')
% scatter(MatchingPoints3(:,1),MatchingPoints3(:,2),'+','r')
% 


% 
%  


 

 % Rotated_X = imrotate(X,27,'bilinear','crop');
%  
% RotatedImage_KeyPoints = SIFT_KeyPoints(Rotated_X);
% 
% KeyPointsImage = zeros(size(X));
% 
% for i = 1:size(Original_KeyPoints)
%    KeyPointsImage(round(Original_KeyPoints(i,2)),round(Original_KeyPoints(i,3)),round(Original_KeyPoints(i,4))) = 100;
% end
% 
% 
% Rotated_KeyPointsImage = imrotate(KeyPointsImage,27,'crop');
% 
% 
% Index1 = 0;
% Index2 = 0;
% Index3 = 0;
% 
% for i = 1:67
% 
%    [Temp1, Temp2, Temp4] = find(Rotated_KeyPointsImage(:,:,i));
%    Temp3 = ones(size(Temp2))*i;
%     
%    Index1 = [Index1 ; Temp1];
%    Index2 = [Index2 ; Temp2];
%    Index3 = [Index3 ; Temp3];
% 
% end
% 
% Rotated_KeyPoints_3 = [Index1 Index2 Index3];
% 
% 
% 
% 
% scatter3(RotatedImage_KeyPoints(: ,2), RotatedImage_KeyPoints(: ,3), RotatedImage_KeyPoints(: ,4), '+')
%  hold on
% 
% 
% scatter3(Index1, Index2, Index3, 'o')
%  
% figure(4)
% subplot (2, 3, 1); imshow(uint8(255*Norm_X(:,:,70)/max(max(max(Norm_X)))));
% subplot (2, 3, 2); imshow(uint8(255*DoG1(:,:,70)/max(max(max(DoG1)))));
% subplot (2, 3, 3); imshow(uint8(255*DoG2(:,:,70)/max(max(max(DoG1)))));
% subplot (2, 3, 4); imshow(uint8(255*DoG3(:,:,70)/max(max(max(DoG1)))));
% subplot (2, 3, 5); imshow(uint8(255*DoG4(:,:,70)/max(max(max(DoG1)))));
% subplot (2, 3, 6); imshow(uint8(255*DoG5(:,:,70)/max(max(max(DoG1)))));
% 
% subplot (2, 3, 1); imshow(uint8(255*double(X(:,:,70))/max(max(max(Norm_X)))));
% subplot (2, 3, 2); imshow(uint8(255*Norm_X(:,:,70)/max(max(max(Norm_X)))));
% subplot (2, 3, 3); imshow(uint8(255*FilteredImage2(:,:,70)/max(max(max(FilteredImage2)))));
% subplot (2, 3, 4); imshow(uint8(255*FilteredImage3(:,:,70)/max(max(max(FilteredImage3)))));
% subplot (2, 3, 5); imshow(uint8(255*FilteredImage4(:,:,70)/max(max(max(FilteredImage4)))));
% subplot (2, 3, 6); imshow(uint8(255*FilteredImage5(:,:,70)/max(max(max(FilteredImage5)))));
% Key
% 
% 
% 
% figure(2)
% imshow(Norm_X(90:110,90:110,50))
