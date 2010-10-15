
clear all
clc
format short

addpath /data/export/home/ray/Images/MNI_Transfered/
addpath /data/export/home/ray/NIFTI_20100106/
addpath /data/export/home/ray/SIFT/mFiles/
addpath /data/export/home/ray/Images/ColumbiaScans/
addpath /data/export/home/ray/SURF
addpath /data/export/home/ray/SIFT/Vedaldi

    image1 = load_untouch_nii('brain3_1.nii');
    X1 = image1.img;

    NormalizationType = 3;
    NgL=256;
    Slice_X1 = double(squeeze(X1(:,:,90)));
    Norm_X1 = ImageNormalizer(Slice_X1, NormalizationType, NgL);
    [M N R] = size(Norm_X1);
    clear X1
    clear image1

SamplePerScale = 3 ; InitialSigma = 2.0159; Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5; Ther = 0.03; SearchArea=3; EdgeDetect=0;

[KeyPoints1 IgnordKeyPoints1]= SIFT_KeyPoints(Norm_X1, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);

Options.verbose=false; 
Ipts=OpenSurf(Norm_X1,Options);
SURF_KeyPoints1 = [Ipts.x ;Ipts.y ;Ipts.scale]';

[frames,Descriptors1,gss,dogss] = sift(Norm_X1,'verbosity', 0);
Vedaldi_Keypoints1 = (frames(1:3,:))';
 
    image1 = load_untouch_nii('brain2_1_TransformedTo3_1.nii');
    X1 = image1.img;

    NormalizationType = 3;
    NgL=256;
    Slice_X2 = double(squeeze(X1(:,:,90)));

    Norm_X2 = ImageNormalizer(Slice_X2, NormalizationType, NgL);
    [M N R] = size(Norm_X2);
    clear X1
    clear image1


[KeyPoints2 IgnordKeyPoints2]= SIFT_KeyPoints(Norm_X2, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
Options.verbose=false; 
Ipts=OpenSurf(Norm_X2,Options);
SURF_KeyPoints2 = [Ipts.x ;Ipts.y ;Ipts.scale]';

[frames,Descriptors2,gss,dogss] = sift(Norm_X2,'verbosity', 0);
Vedaldi_Keypoints2 = (frames(1:3,:))';

DesrciptorMatches = siftmatch(Descriptors1, Descriptors2, 1.5)

for i = 1:length(DesrciptorMatches)
    point1(i,:) = Vedaldi_Keypoints1(DesrciptorMatches(1,i),1:2);
    point2(i,:) = Vedaldi_Keypoints2(DesrciptorMatches(2,i),1:2);
end

 figure(1);
 subplot(1,2,1);
 imshow((imrotate((Norm_X1),0)),[min(Norm_X1(:)) max(Norm_X1(:))]);
 hold on
 plot(point1(:,1), point1(:,2), 'r+', 'MarkerSize', 5)
 hold on
 plot(Vedaldi_Keypoints1(:,1), Vedaldi_Keypoints1(:,2), 'co', 'MarkerSize', 5)
 title('Subject2_1 / 403 Keypoints(cyans) / 43 match for ther=1.5') 


subplot(1,2,2);
 imshow((imrotate((Norm_X2),0)),[min(Norm_X2(:)) max(Norm_X2(:))]);
 hold on
 plot(point1(:,1), point1(:,2), 'r+', 'MarkerSize', 5)
 hold on
 plot(Vedaldi_Keypoints2(:,1), Vedaldi_Keypoints2(:,2), 'co', 'MarkerSize', 5)
 title('Subject3_1 / 455 Keypoints(cyans) / 43 match for ther=1.5') 



%  
%     image1 = load_nii('Subject1_UMN1_SkullStripped_Transformed_MNI152_T1_1mm.nii');
%     X1 = image1.img;
% 
%     NormalizationType = 3;
%     NgL=256;
%     Slice_X1 = double(squeeze(X1(:,:,90)));
% 
%     Norm_X1 = ImageNormalizer(Slice_X1, NormalizationType, NgL);
%     [M N R] = size(Norm_X1);
%     clear X1
%     clear image1
% 
% 
% [KeyPoints3 IgnordKeyPoints3]= SIFT_KeyPoints(Norm_X1, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
% Options.verbose=false; 
% Ipts=OpenSurf(Norm_X1,Options);
% SURF_KeyPoints3 = [Ipts.x ;Ipts.y ;Ipts.scale]';
% [frames,descriptors,gss,dogss] = sift(Norm_X1,'verbosity', 0);
% Vedaldi_Keypoints3 = (frames(1:3,:))';
% 
%   
%     image1 = load_nii('Subject1_UMN2_SkullStripped_Transformed_MNI152_T1_1mm.nii');
%     X1 = image1.img;
% 
%     NormalizationType = 3;
%     NgL=256;
%     Slice_X1 = double(squeeze(X1(:,:,90)));
% 
%     Norm_X1 = ImageNormalizer(Slice_X1, NormalizationType, NgL);
%     [M N R] = size(Norm_X1);
%     clear X1
%     clear image1
% 
% 
% [KeyPoints4 IgnordKeyPoints4]= SIFT_KeyPoints(Norm_X1, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
% Options.verbose=false; 
% Ipts=OpenSurf(Norm_X1,Options);
% SURF_KeyPoints4 = [Ipts.x ;Ipts.y ;Ipts.scale]';
% [frames,descriptors,gss,dogss] = sift(Norm_X1,'verbosity', 0);
% Vedaldi_Keypoints4 = (frames(1:3,:))';
% 
[MatchPersentage1 MatchingPoints1] = KeyPointsMatching(KeyPoints1, KeyPoints2, 2);
% [MatchPersentage2 MatchingPoints2] = KeyPointsMatching(KeyPoints1, KeyPoints3, 2);
% [MatchPersentage3 MatchingPoints3] = KeyPointsMatching(KeyPoints1, KeyPoints4, 2);
% [MatchPersentage4 MatchingPoints4] = KeyPointsMatching(KeyPoints2, KeyPoints3, 2);
% [MatchPersentage5 MatchingPoints5] = KeyPointsMatching(KeyPoints2, KeyPoints4, 2);
% [MatchPersentage6 MatchingPoints6] = KeyPointsMatching(KeyPoints3, KeyPoints4, 2);
% 

[SMatchPersentage1 SMatchingPoints1] = KeyPointsMatching(SURF_KeyPoints1, SURF_KeyPoints2, 2);
% [SMatchPersentage2 SMatchingPoints2] = KeyPointsMatching(SURF_KeyPoints1, SURF_KeyPoints3, 2);
% [SMatchPersentage3 SMatchingPoints3] = KeyPointsMatching(SURF_KeyPoints1, SURF_KeyPoints4, 2);
% [SMatchPersentage4 SMatchingPoints4] = KeyPointsMatching(SURF_KeyPoints2, SURF_KeyPoints3, 2);
% [SMatchPersentage5 SMatchingPoints5] = KeyPointsMatching(SURF_KeyPoints2, SURF_KeyPoints4, 2);
% [SMatchPersentage6 SMatchingPoints6] = KeyPointsMatching(SURF_KeyPoints3, SURF_KeyPoints4, 2);

[VMatchPersentage1 VMatchingPoints1] = KeyPointsMatching(Vedaldi_Keypoints1, Vedaldi_Keypoints2, 2);
% [VMatchPersentage2 VMatchingPoints2] = KeyPointsMatching(Vedaldi_Keypoints1, Vedaldi_Keypoints3, 2);
% [VMatchPersentage3 VMatchingPoints3] = KeyPointsMatching(Vedaldi_Keypoints1, Vedaldi_Keypoints4, 2);
% [VMatchPersentage4 VMatchingPoints4] = KeyPointsMatching(Vedaldi_Keypoints2, Vedaldi_Keypoints3, 2);
% [VMatchPersentage5 VMatchingPoints5] = KeyPointsMatching(Vedaldi_Keypoints2, Vedaldi_Keypoints4, 2);
% [VMatchPersentage6 VMatchingPoints6] = KeyPointsMatching(Vedaldi_Keypoints3, Vedaldi_Keypoints4, 2);
% 
 
%  
%   Rotated_X = imrotate(Norm_X1, 20, 'bilinear', 'crop');
%   Rotated_Norm_X = imrotate(Rotated_X, -20, 'bilinear', 'crop');
%   
%   NormalizationType = 3;
%   NgL=256;
%   Rotated_Norm_X = ImageNormalizer(Rotated_Norm_X, NormalizationType, NgL);
%  
%   [KeyPoints2 IgnordKeyPoints2]= SIFT_KeyPoints(Rotated_Norm_X, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
% 
%   
figure(1)
imshow(Norm_X1,[min(Norm_X1(:)) max(Norm_X1(:))]);
 hold on
 plot(KeyPoints1(:,1), KeyPoints1(:,2), 'c+', 'MarkerSize', 5)
 hold on 
 plot(KeyPoints2(:,1), KeyPoints2(:,2), 'ro', 'MarkerSize', 5)
%  hold on 
%  plot(KeyPoints3(:,1), KeyPoints3(:,2), 'm*', 'MarkerSize', 5)
%  hold on 
%  plot(KeyPoints4(:,1), KeyPoints4(:,2), 'gs', 'MarkerSize', 5)
% 
 figure(2)
 imshow((imrotate((Norm_X1),0)),[min(Norm_X1(:)) max(Norm_X1(:))]);
 hold on
 plot(SURF_KeyPoints1(:,1), SURF_KeyPoints1(:,2), 'c+', 'MarkerSize', 5)
 hold on 
 plot(SURF_KeyPoints2(:,1), SURF_KeyPoints2(:,2), 'ro', 'MarkerSize', 5)
 hold on 
%  plot(SURF_KeyPoints3(:,1), SURF_KeyPoints3(:,2), 'm*', 'MarkerSize', 5)
%  hold on 
%  plot(SURF_KeyPoints4(:,1), SURF_KeyPoints4(:,2), 'gs', 'MarkerSize', 5)
% 
 
  figure(3)
 imshow((imrotate((Norm_X1),0)),[min(Norm_X1(:)) max(Norm_X1(:))]);
 hold on
 plot(Vedaldi_Keypoints1(:,1), Vedaldi_Keypoints1(:,2), 'c+', 'MarkerSize', 5)
 hold on 
 plot(VMatchingPoints1(:,1), VMatchingPoints1(:,2), 'ro', 'MarkerSize', 5)
 hold on 
  figure(4)
 imshow((imrotate((Norm_X2),0)),[min(Norm_X2(:)) max(Norm_X2(:))]);
 hold on
 plot(Vedaldi_Keypoints2(:,1), Vedaldi_Keypoints2(:,2), 'c+', 'MarkerSize', 5)
 hold on 
 plot(VMatchingPoints1(:,1), VMatchingPoints1(:,2), 'ro', 'MarkerSize', 5)
 hold on 
%  plot(Vedaldi_Keypoints3(:,1), Vedaldi_Keypoints3(:,2), 'm*', 'MarkerSize', 5)
%  hold on 
%  plot(Vedaldi_Keypoints4(:,1), Vedaldi_Keypoints4(:,2), 'gs', 'MarkerSize', 5)
% 
%   
%  InitialSigma = 2;
% [MatchPersentage MatchingPoints] = KeyPointsMatching(KeyPoints1, KeyPoints2, InitialSigma); 
% 


%  Hold on
%  plot(KeyPoints2(:,1), KeyPoints2(:,2), 'co', 'MarkerSize', 5)
%  
% %    figure(2)
% %    imshow((imrotate(fliplr(Norm_X2),90)),[min(Norm_X2(:)) max(Norm_X2(:))]);
% %    hold on
% plot(KeyPoints2(:,1), KeyPoints2(:,2), 'co', 'MarkerSize', 5)
% % 
%  Rotation = (round(20*rand)) - 10;
%  Rotation = 10;
%  Norm_X3 = VolumeRotator(Norm_X1,Rotation);
%  KeyPoints3 = SIFT_KeyPoints(Norm_X3, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea );
%  RotatedPoint = PointRotator(KeyPoints3(:,1:2), -Rotation , [M/2 N/2]);
%  RotatedPoint_10 = [RotatedPoint KeyPoints3(:,3)];
% 
%  [MatchPersentage2 MatchingPoints2] = KeyPointsMatching(KeyPoints1, RotatedPoint_10, InitialSigma);
%  
%  hold on
%  plot(RotatedPoint_10(:,1), RotatedPoint_10(:,2), 'co', 'MarkerSize', 5)


% 
% 
% Norm_X = ImageNormalizer(Slide_X, NormalizationType, NgL);
% [M N R] = size(Norm_X);
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%% Keypoints of the Original Image 
% SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.01;
% KeyPoints1 = SIFT_KeyPoints_1(Norm_X, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther);
% 
% 
% %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by 10 degree
%  Norm_X2 = VolumeRotator(Norm_X,15);
%  SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.01;
%  KeyPoints2 = SIFT_KeyPoints_1(Norm_X2, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther);
%  RotatedPoint = PointRotator(KeyPoints2(:,1:2), -15, [M/2 N/2]);
%  RotatedPoint10 = [RotatedPoint KeyPoints2(:,3)];
%  
%  
%  %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by -10 degree
%  Norm_X3 = VolumeRotator(Norm_X,-15);
%  SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.01;
%  KeyPoints3 = SIFT_KeyPoints_1(Norm_X3, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther );
%  RotatedPoint = PointRotator(KeyPoints3(:,1:2), 15 , [M/2 N/2]);
%  RotatedPoint_10 = [RotatedPoint KeyPoints3(:,3)];
% 
%  
%  %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by 5 degree
%  Norm_X4 = VolumeRotator(Norm_X,40);
%  SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.01;
%  KeyPoints4 = SIFT_KeyPoints_1(Norm_X4, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther );
%  RotatedPoint = PointRotator(KeyPoints4(:,1:2), -40, [M/2 N/2]);
%  RotatedPoint5 = [RotatedPoint KeyPoints4(:,3)];
%  
% 
%  InitialSigma = 2;
%  
% %%%%%%%%%%%%%%%%%%%%% Computation of the Matches  
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
% plot(KeyPoints1(:,1), KeyPoints1(:,2), 'r+', 'MarkerSize', 5)
% 
% figure(2)
% imshow((imrotate(fliplr(Norm_X2),90)),[min(Norm_X(:)) max(Norm_X(:))]);
% hold on
% plot(KeyPoints2(:,1), KeyPoints2(:,2), 'r+', 'MarkerSize', 5)
% 
% figure(3)
% imshow((imrotate(fliplr(Norm_X3),90)),[min(Norm_X(:)) max(Norm_X(:))]);
% hold on
% plot(KeyPoints3(:,1), KeyPoints3(:,2), 'r+', 'MarkerSize', 5)
% 
% figure(4)
% imshow((imrotate(fliplr(Norm_X4),90)),[min(Norm_X(:)) max(Norm_X(:))]);
% hold on
% plot(KeyPoints4(:,1), KeyPoints4(:,2), 'r+', 'MarkerSize', 5)
% 
% figure(5)
% imshow((imrotate(fliplr(Norm_X),90)),[min(Norm_X(:)) max(Norm_X(:))]);
% hold on
% plot(MatchingPoints4(:,1), MatchingPoints4(:,2), 'r+', 'MarkerSize', 5)
% 
