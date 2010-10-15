clear all
clc

x = zeros(1000,1000);
x(:, 1:2:1000) = 256;


imshow(x, [0 256]);