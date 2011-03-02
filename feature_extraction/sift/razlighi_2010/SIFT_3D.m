
clear all
clc

%%%%%%%% Please Run this code on different Matlab 
%%%%%%%% For each run change the InitialSigma value to one of these 
%%%%%%%% values {1.2, 1.4, 1.6, 1.8} and save the result workspace
%%%%%%%% with the name of the following line
InitialSigma=1.6;SamplePerScale=3;Nominal_Sigma=.5;StartingOctave=-1;LastOctave=6;Ther=0.00;SearchArea=3;EdgeDetect=0;



%%%%%%%% Filename Setting for one round
Subjects = ['Subject1'; 'Subject3'; 'Subject5']
Locations = ['Iowa1'; 'Iowa2'; 'UMN1 '; 'UMN2 ']

for i = 3
    for j=1:4
        
        if (j==3 || j==4) 
            temp = Locations(j,1:4);
        else
            temp = Locations(j,:);
        end
        
        FileName = ['/data/export/home/ray/SIFT/SIFT_Images/SkullSTripped/' Subjects(i,:) '_'  temp '_SkullStripped.nii'];
        
        image1 = load_nii(FileName);
        
        X = image1.img;
        
        fprintf('\n A new volume is loaded from %s\n', FileName );
        
        [M N R] = size(X);
        fprintf(['\n Volume size is: \n' int2str(M) 'x' int2str(N) 'x' int2str(R) '\n']);
        
        NormalizationType = 3;
        NgL=256;

        Norm_X = ImageNormalizer(X, NormalizationType, NgL);
        [M N R] = size(Norm_X);
        clear image1
        clear X


 %%%%%%%%%%%%%%%%%%%% Keypoints of the Original Image 
 fprintf('\nExtract Keypoints with parameeter:\nInitialSigma= %G\nSamplePerScale= %G\nNominal_Sigma= %G\nStartingOctave= %G\nLastOctave=%G\nTher= %G\nSearchArea= %G\nEdgeDetect= %G\n', InitialSigma,SamplePerScale,Nominal_Sigma,StartingOctave,LastOctave,Ther,SearchArea,EdgeDetect);
 KeyPoints1 = SIFT_KeyPoints(Norm_X, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);


 %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by random degree in
 %%%%%%%%%%%%%%%%%%%% the range of [-10 10] and along all three axes
 Rotation = (round(40*rand(1,3)));
 Norm_X2 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
 KeyPoints2 = SIFT_KeyPoints(Norm_X2, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
 RotatedPoint1 = PointRotator(KeyPoints2(:,1:3), [-Rotation(1,1) Rotation(1,2) -Rotation(1,3)], [M/2 N/2 R/2]);
 RotatedPointPlus = [RotatedPoint1 KeyPoints2(:,4)];

 %%%%%%%%%%%%%%%%%%%% KeyPoints of the rotated image by random degree in
 %%%%%%%%%%%%%%%%%%%% the range of [-10 10] and along all three axes
 Rotation = -Rotation;
 Norm_X3 = VolumeRotator(Norm_X,[Rotation(1,1) Rotation(1,2) Rotation(1,3)]);
 KeyPoints3 = SIFT_KeyPoints(Norm_X3, SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect);
 RotatedPoint2 = PointRotator(KeyPoints3(:,1:3), [-Rotation(1,1) Rotation(1,2) -Rotation(1,3)], [M/2 N/2 R/2]);
 RotatedPointMinus = [RotatedPoint2 KeyPoints3(:,4)];

 
 
 InitialSigma = 2;
 
% %%%%%%%%%%%%%%%%%%%%%% Computation of the Matches  
[MatchPersentage1 MatchingPoints1] = KeyPointsMatching(KeyPoints1, RotatedPointPlus, InitialSigma);

[MatchPersentage2 MatchingPoints2] = KeyPointsMatching(KeyPoints2, RotatedPointMinus, InitialSigma);

[MatchPersentage(i,j) MatchingPoints(i,j).keys] = KeyPointsMatching(MatchingPoints1, MatchingPoints2, InitialSigma);


    end
end

 


