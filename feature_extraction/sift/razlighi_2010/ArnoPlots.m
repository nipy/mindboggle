
clear all
clc




% rand('seed',0)
% %data = rand(10,10,10);
image1 = load_nii('/data/export/home/ray/SIFT/SIFT_Images/anat3T_subject1.nii');
%image1 = load_nii('/Users/ray/Research/SIFT/T1.img');
X = image1.img;

NormalizationType = 3;
NgL=256;

Slide_X = X(:,:,70);
Norm_X = ImageNormalizer(Slide_X, NormalizationType, NgL);
[M N R] = size(Norm_X);
%imwrite(Norm_X,'/Users/ray/Research/SIFT/TestImage.bmp','bmp'); 

clear image1
clear X



%%%%%%%%%%%%%%%%%%%% Keypoints of the Original Image 
FilterSize = 7; SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.0053;
KeyPoints1 = SIFT_KeyPoints(Norm_X, FilterSize, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther);


%%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by 10 degree
 Norm_X2 = VolumeRotator(Norm_X,15);
 FilterSize = 7; SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.0053;
 KeyPoints2 = SIFT_KeyPoints(Norm_X2, FilterSize, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther);
 RotatedPoint = PointRotator(KeyPoints2(:,1:2), -15, [M/2 N/2]);
 RotatedPoint10 = [RotatedPoint KeyPoints2(:,3)];
 
 
 %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by -10 degree
 Norm_X3 = VolumeRotator(Norm_X,-15);
 FilterSize = 7; SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.0053;
 KeyPoints3 = SIFT_KeyPoints(Norm_X3, FilterSize, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther );
 RotatedPoint = PointRotator(KeyPoints3(:,1:2), 15 , [M/2 N/2]);
 RotatedPoint_10 = [RotatedPoint KeyPoints3(:,3)];

 
 %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by 5 degree
 Norm_X4 = VolumeRotator(Norm_X,40);
 FilterSize = 7; SamplePerScale = 3; InitialSigma = 1.6*2^(1/3); Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;
 KeyPoints4 = SIFT_KeyPoints(Norm_X4, FilterSize, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther );
 RotatedPoint = PointRotator(KeyPoints4(:,1:2), -40, [M/2 N/2]);
 RotatedPoint5 = [RotatedPoint KeyPoints4(:,3)];
 

 InitialSigma = 2;
 
%%%%%%%%%%%%%%%%%%%%%% Computation of the Matches  
[MatchPersentage1 MatchingPoints1] = KeyPointsMatching(KeyPoints1, RotatedPoint10, InitialSigma);

[MatchPersentage2 MatchingPoints2] = KeyPointsMatching(KeyPoints1, RotatedPoint_10, InitialSigma);
 
[MatchPersentage3 MatchingPoints3] = KeyPointsMatching(MatchingPoints1, MatchingPoints2, InitialSigma);
 
 
[MatchPersentage4 MatchingPoints4 NearestDistance] = KeyPointsMatching(MatchingPoints3, RotatedPoint5, InitialSigma);


%%%%%%%%%%%%%%%%%%%%%% Visualization of the matches
figure(1)
imshow((imrotate(fliplr(Norm_X),90)),[min(Norm_X(:)) max(Norm_X(:))]);
hold on
plot(KeyPoints1(:,1), KeyPoints1(:,2), 'r+', 'MarkerSize', 10)

figure(2)
imshow((imrotate(fliplr(Norm_X2),90)),[min(Norm_X(:)) max(Norm_X(:))]);
hold on
plot(KeyPoints2(:,1), KeyPoints2(:,2), 'r+', 'MarkerSize', 10)

figure(3)
imshow((imrotate(fliplr(Norm_X3),90)),[min(Norm_X(:)) max(Norm_X(:))]);
hold on
plot(KeyPoints3(:,1), KeyPoints3(:,2), 'r+', 'MarkerSize', 10)

figure(4)
imshow((imrotate(fliplr(Norm_X4),90)),[min(Norm_X(:)) max(Norm_X(:))]);
hold on
plot(KeyPoints4(:,1), KeyPoints4(:,2), 'r+', 'MarkerSize', 10)

figure(5)
imshow((imrotate(fliplr(Norm_X),90)),[min(Norm_X(:)) max(Norm_X(:))]);
hold on
plot(MatchingPoints4(:,1), MatchingPoints4(:,2), 'r+', 'MarkerSize', 10)

