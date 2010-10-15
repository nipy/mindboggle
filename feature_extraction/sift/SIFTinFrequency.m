clear all
clc
addpath /data/export/home/ray/NIFTI_20100106/
addpath /data/export/home/ray/SIFT/mFiles/


image1 = load_nii('/data/export/home/ray/Images/Standard/MNI152_T1_1mm_brain.nii');
image2 = load_nii('/data/export/home/ray/Images/Standard/MNI152_T1_1mm_brain_mask.nii');

X = image1.img;
Mask_X = int16(image2.img);

[m n r] = size(X);
% NormalizationType = 2;
% NgL=256;
% 
% Norm_X = squeeze(X(:,70,:));
% Norm_X = ImageNormalizer(Norm_X, NormalizationType, NgL);
% figure(1); imshow(Norm_X);

MeanValue = sum(sum(sum(X .* Mask_X)))/sum(sum(sum(Mask_X)))
UnBiased_X = (X - MeanValue * Mask_X);    
        

Norm_X = squeeze(UnBiased_X(:,120,:));

figure(1); imshow(Norm_X);

F1 = fftn(double(Norm_X));
F = fftshift(F1);

Abs_F = abs(F);
figure(2); imshow(Abs_F, [min(Abs_F(:)) max(Abs_F(:))])

OneD_X = Norm_X(93,:);
F1 = fftn(double(OneD_X));
F = fftshift(F1);
figure (3)
plot(OneD_X);
hold on
grid on
plot(abs(F), 'r')

NoBorder_X = OneD_X(26:160);
F1 = fftn(double(NoBorder_X));
F = fftshift(F1);
figure (4)
plot(NoBorder_X);
hold on
grid on
plot(abs(F), 'r')




% mesh(Abs_F)
% figure(2)
%imshow(Abs_F,[min(Abs_F(:)) (max(Abs_F(:))])
% ImageSize = size(F1);
% 
% figure(3);
%TransferFunction = ImageFilterTransferFunction('Type','Ideal LPF', 'Bandwith', 20, 'Filter Order', 2, 'Image Size', ImageSize);

%imshow(TransferFunction)
% plot(TransferFunction(97,:))
% hold on
%  TransferFunction_20 = ImageFilterTransferFunction('Type','Butterworth LPF', 'Bandwith', 20, 'Filter Order', 10, 'Image Size', ImageSize);
%  HPF_TransferFunction_20 = 1 - TransferFunction_20;
%  
%  TransferFunction_40 = ImageFilterTransferFunction('Type','Butterworth LPF', 'Bandwith', 60, 'Filter Order', 10, 'Image Size', ImageSize);
%  HPF_TransferFunction_40 = 1 - TransferFunction_40;
% 
% BPF = HPF_TransferFunction_20 - HPF_TransferFunction_40;
% 
%  plot(HPF_TransferFunction_20(129,:))
%  hold on
%   plot(HPF_TransferFunction_40(129,:),'g')
% hold on
%   plot(BPF(129,:),'r')
% 
%  grid on
% 
 
 % TransferFunction = ImageFilterTransferFunction('Type','Gaussian LPF', 'Bandwith', 20, 'Filter Order', 2, 'Image Size', ImageSize);
% plot(TransferFunction(97,:))
% hold on

% 
% for i=1:16
%     
%     
%     if i~=1
%         Bandwith(i) = Sigma(i-1) - Sigma(i);
%         CenterFrequency(i) = Sigma(i) + Bandwith(i)/2;
%         TransferFunction = ImageFilterTransferFunction('Type','Butterworth BPF', 'Bandwith', Bandwith(i), 'Filter Order', 20, 'Image Size', ImageSize, 'Center Frequency', CenterFrequency(i));
%         plot(TransferFunction(128,:))
%         hold on
% 
%     end
% end

%xlim =[128 256];

clear all

figure(5)
TestImage = zeros(512,512);
TestImage(10:20,10:20) = 1; 
F = fftshift(fftn(TestImage));
mesh(abs(F))
figure (6)
mesh(angle(F))
