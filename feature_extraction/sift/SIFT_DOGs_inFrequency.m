clear all
clc

iter = .5
for k=1:1
Sigma = 1 + .2*3
S = 3;

figure(k)

u = 0:.001: iter;
v = 0:.001: iter;

for i=1:(4*S+1)
    Sig = Sigma*2^((i-5)/S)       
    
    
   
    F(i,:) = exp(-2*pi^2*Sig^2*u.^2 );
    
    
    
    if i~=1
        d(i,:) = F(i-1,:) - F(i,:);
        
        if (i>1 && i<(S+2))
            plot(u, d(i,:), 'b:', 'MarkerSize', 30 )
            hold on
        elseif (i>(S+1) && i<(2*S+2))
            plot(u, d(i,:), 'b.-', 'MarkerSize', 5 )
            hold on
        elseif (i>(2*S+1) && i<(3*S+2))
            plot(u, d(i,:), 'b--', 'MarkerSize', 10 )
            hold on
        elseif (i>(3*S+1) && i<(4*S+2))
            plot(u, d(i,:), 'r-', 'MarkerSize', 10 )
            hold on
        end

    end
    %s(i,:) = ((2*pi*Sigma^2)^-.5)*exp(-(x.^2 )/(2*Sigma^2));
end
grid on 
xlabel('Frequency','fontsize',15);
ylabel('Magnitude','fontsize',15);

%set(gca, 'FontSize', 15);


end
%legend('OCtave -1','OCtave -1','OCtave -1','Octave 0','Octave 0','Octave 0', 'OCtave 1','OCtave 1','OCtave 1', 'Octave 2', 'Octave 2', 'Octave 2','fontsize',20)

% plot(d(i,:),'r')
% hold on
% plot(s,'b')
% 
% 


% 
% 
% 
% % rand('seed',0)
% % %data = rand(10,10,10);
% image1 = load_nii('/data/export/home/ray/SIFT/SIFT_Images/anat3T_subject1.nii');
% %image1 = load_nii('/Users/ray/Research/SIFT/T1.img');
% X = image1.img;
% 
% NormalizationType = 3;
% NgL=256;
% 
% Slide_X = X(:,:,70);
% Norm_X = ImageNormalizer(Slide_X, NormalizationType, NgL);
% [M N R] = size(Norm_X);
% %imwrite(Norm_X,'/Users/ray/Research/SIFT/TestImage.bmp','bmp'); 
% 
% clear image1
% clear X
% 
% 
% 
