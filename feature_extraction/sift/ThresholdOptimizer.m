%function [Threshold] = ThresholdOptimizer (InputImage)

clear all
clc
addpath /data/export/home/ray/NIFTI_20100106/
addpath /data/export/home/ray/SIFT/mFiles/

    image1 = load_untouch_nii('brain3_1.nii');
    X1 = image1.img;

    NormalizationType = 3;
    NgL=256;
    Slice_X1 = double(squeeze(X1(:,:,80)));
    Norm_X1 = ImageNormalizer(Slice_X1, NormalizationType, NgL);
    [M N R] = size(Norm_X1);
    clear X1
    clear image1
InputImage = Norm_X1;



SamplePerScale = 3 ; InitialSigma = 2.0159; Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5; Ther = 0.001; SearchArea=3; EdgeDetect=0;
 Rotated_X = imrotate(InputImage, 10, 'bicubic', 'crop');
 NormalizationType = 3;
 NgL=256;
 Rotated_X = ImageNormalizer(Rotated_X, NormalizationType, NgL);


for i = 1:10
  
   Ther = i/100;
  
  KeyPoints = SIFT_KeyPoints(InputImage, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
   
  [frames,descriptors,gss,dogss] = sift(InputImage,'threshold', Ther, 'verbosity', 0);
  Vedaldi_Keypoints = (frames(1:3,:))';

  
  Rotated_KeyPoints_SIFT = SIFT_KeyPoints(Rotated_X, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
  Rotated_Back_SIFT_Point = PointRotator(Rotated_KeyPoints_SIFT(:, 1:2), -10, [M/2 N/2]);
  
  [frames,descriptors,gss,dogss] = sift(Rotated_X,'threshold', Ther, 'verbosity', 0);
  Rotated_Vedaldi_Keypoints = (frames(1:3,:))';
  Rotated_Back_Vedaldi_Point = PointRotator(Rotated_Vedaldi_Keypoints(:, 1:2), 10, [M/2 N/2]);

  NumberOfKeypoints_SIFT(i) = length(KeyPoints) ;
  NumberOfKeypoints_Vedal(i) = length(Vedaldi_Keypoints);
  
  if NumberOfKeypoints_SIFT(i) >10
      [MatchPersentageSIFT(i) MatchingPoints] = KeyPointsMatching(KeyPoints(:,1:3), [Rotated_Back_SIFT_Point Rotated_KeyPoints_SIFT(:,3)] , InitialSigma);
       [MatchPersentageVedal(i) MatchingPoints] = KeyPointsMatching(Vedaldi_Keypoints(:,1:3), [Rotated_Back_Vedaldi_Point Rotated_Vedaldi_Keypoints(:,3)] , InitialSigma);
  end
end

 figure(1);
 plot( NumberOfKeypoints_SIFT, 'b');
 hold on
 plot(NumberOfKeypoints_Vedal, 'r');
 figure(2)
 plot(MatchPersentageSIFT, 'b');
 hold on
 plot(MatchPersentageVedal,'r');
 
 Threshold = Ther;
 
 
