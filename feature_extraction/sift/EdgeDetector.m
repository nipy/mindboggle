clear
clc

addpath 'C:\Users\Ray\Documents\My Dropbox\Research\SIFT\mFiles'
addpath 'C:\Users\Ray\Documents\My Dropbox\Research\Images\NonMedicals'
addpath 'C:\Users\Ray\Documents\My Dropbox\Research\SIFT\Vedaldi'
addpath 'C:\Users\Ray\Documents\My Dropbox\Research\NIFTI_20100106'

addpath /Users/ray/Dropbox/Research/SIFT/mFiles/
addpath /Users/ray/Dropbox/Research/Images/NonMedicals/
addpath /Users/ray/Dropbox/Research/Images/ColumbiaScans
addpath /Users/ray/Dropbox/Research/SIFT/Vedaldi
addpath /Users/ray/Dropbox/Research/NIFTI_20100106


 image1 = load_untouch_nii('brain2_1.nii');
 X = image1.img;
 Slice_X = double(X(:, :, 80));
 %SmoothedImage = imfilter(Slice_X, GaussianKernel(1.3, 1), 'replicate');
 %Slice_X = SmoothedImage;
 Slice_X = doublesize(Slice_X);
 [M N R] = size(Slice_X);
 
 Der_X = (circshift(Slice_X,[-1 0 0]) - circshift(Slice_X,[1 0 0]))/2;
 Der_Y = (circshift(Slice_X,[0 -1 0]) - circshift(Slice_X,[0 1 0]))/2;
 Der_Z = (circshift(Slice_X,[0 0 -1]) - circshift(Slice_X,[0 0 1]))/2;
 
 Der_XX = (circshift(Slice_X,[-1 0 0]) -2*Slice_X + circshift(Slice_X,[1 0 0]));
 Der_YY = (circshift(Slice_X,[0 -1 0]) -2*Slice_X + circshift(Slice_X,[0 1 0]));
 Der_ZZ = (circshift(Slice_X,[0 0 -1]) -2*Slice_X + circshift(Slice_X,[0 0 1]));
 
 Gradient = sqrt (Der_X.^2 + Der_Y.^2 + Der_Z.^2);
 
 Laplacian = sqrt (Der_XX.^2 + Der_YY.^2 + Der_ZZ.^2);
 
 DifferentialGeometric = (Der_X.^2).*Der_XX + 2.*Der_X.*Der_Y.*Der_XY + Der_Y.^2.*Der_YY;
 
 
 BorderImage_Gra = uint8(boolean(floor(Gradient/150)));
 BorderImage_Lap = uint8(boolean(floor(Laplacian/50)));
 
 
  
 figure(1)
 subplot(3,3,1)
    imshow(Slice_X,[min(Slice_X(:)) max(Slice_X(:))]);
 
subplot(3,3,2)
    imshow(Gradient,[min(Gradient(:)) max(Gradient(:))]);

 for i=1:7
 BorderImage_Gra = uint8(boolean(floor(Gradient/(i*15))));
 
subplot(3,3,i+2)
    imshow(BorderImage_Gra,[min(BorderImage_Gra(:)) max(BorderImage_Gra(:))]);
 end

 
 

 
 
figure(2) 
subplot(3,3,1)
    imshow(Slice_X,[min(Slice_X(:)) max(Slice_X(:))]);
 
subplot(3,3,2)
    imshow(Laplacian,[min(Laplacian(:)) max(Laplacian(:))]);

 for i=1:7
 BorderImage_Lap = uint8(boolean(floor(Laplacian/(i*8))));
 
subplot(3,3,i+2)
    imshow(BorderImage_Lap,[min(BorderImage_Lap(:)) max(BorderImage_Lap(:))]);
 end
 
    
    