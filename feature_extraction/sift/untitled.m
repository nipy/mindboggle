clear all
clc

Input_Image1 = imread('/data/export/home/ray/Images/NonMedicals/Mollas.jpg'); 

[m n r] = size(Input_Image1)
for i =1:r
    Input_Image(:,:,i) = flipud(Input_Image1(:,:,i));
end
[m n r] = size(Input_Image)




figure(1); imshow(Input_Image1)
figure(2); mesh(double(Input_Image1))
figure(3); mesh(double(Input_Image))
xlabel('X');
ylabel('Y');

