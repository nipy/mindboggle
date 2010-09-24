clear all
clc
format short

imread('C:\Users\Ray\Documents\My Dropbox\Research\Images\NonMedicals\Mollas.jpg');





addpath /data/export/home/ray/Images/MNI_Transfered/
addpath /data/export/home/ray/NIFTI_20100106/
addpath /data/export/home/ray/SIFT/mFiles/
addpath /data/export/home/ray/Images
addpath /data/export/home/ray/SURF
addpath /data/export/home/ray/SIFT/Vedaldi


%    N = 1000*GaussianKernel(1.8, 4);
%    N = N(4:10,4:10,4:10);

%N=round(100*rand(5,5,5));

DeltaX = -0.0450;
DeltaY = 0.0120;
DeltaZ = -0.60;

x=(-3+DeltaX):(3+DeltaX);
y=(-3+DeltaY):(3+DeltaY);
z=(-3+DeltaZ):(3+DeltaZ);


for i=1:7
    for j=1:7
        for k=1:7
                N(i,j,k) = exp( -( x(i)^2 + y(j)^2 + z(k)^2 ) );
        end
    end
end

% DeltaX = -0.8450;
% DeltaY = 0.0120;
% DeltaZ = -0.050;
% DeltaW = 0.82;
% 
% x=(-3+DeltaX):(3+DeltaX);
% y=(-3+DeltaY):(3+DeltaY);
% z=(-3+DeltaZ):(3+DeltaZ);
% w=(-3+DeltaW):(3+DeltaW);
% 
% 
% for i=1:7
%     for j=1:7
%         for k=1:7
%             for l=1:7
%                 N(i,j,k,l) = exp( -( x(i)^2 + y(j)^2 + z(k)^2 + w(l)^2 ) );
%             end
%         end
%     end
% end

% FirstDerivative_X = ImageDerivatives(N, 'x');
% FirstDerivative_Y = ImageDerivatives(N, 'y');
% FirstDerivative_Z = ImageDerivatives(N, 'z');
% 
% SecondDerivative_XX = ImageDerivatives(FirstDerivative_X, 'x');
% SecondDerivative_XY = ImageDerivatives(FirstDerivative_X, 'y');
% SecondDerivative_XZ = ImageDerivatives(FirstDerivative_X, 'z');
% SecondDerivative_YY = ImageDerivatives(FirstDerivative_Y, 'y');
% SecondDerivative_YX = ImageDerivatives(FirstDerivative_Y, 'x');
% SecondDerivative_YZ = ImageDerivatives(FirstDerivative_Y, 'z');
% SecondDerivative_ZZ = ImageDerivatives(FirstDerivative_Z, 'z');
% SecondDerivative_ZX = ImageDerivatives(FirstDerivative_Z, 'x');
% SecondDerivative_ZY = ImageDerivatives(FirstDerivative_Z, 'y');
% 
% J = [ FirstDerivative_X(4,4,4); FirstDerivative_Y(4,4,4); FirstDerivative_Z(4,4,4)]
% 
% H = [ SecondDerivative_XX(4,4,4) SecondDerivative_XY(4,4,4) SecondDerivative_XZ(4,4,4)
%       SecondDerivative_YX(4,4,4) SecondDerivative_YY(4,4,4) SecondDerivative_YZ(4,4,4)
%       SecondDerivative_ZX(4,4,4) SecondDerivative_ZY(4,4,4) SecondDerivative_ZZ(4,4,4) ]
% 
% Extremum = - (H^(-1))*J
% NewMaxima = N(4,4,4) + .5*J'*Extremum


%InputVolume = N(3:5,3:5,3:5,3:5);
InputVolume = N(3:5,3:5,3:5);

[Adjustment NewMaxima Edge]= TaylorExtremum (InputVolume)
%InputVolume = N(4,4,4,4)
InputVolume = N(4,4,4)











% x=-3:3;
% y=-3:3;
% z=-3:3;
% 
% for i=1:7
%     for j=1:7
%         for k=1:7
%             N(i,j,k) = (x(i)^2 + y(j)^2 + z(k)^2 +1)^(1);
%         end
%     end
% end
% 
% 
% 


% 
% N = ones(5,5,5);
% N(:,1,:) = ones(5,5)*5;
% N(:,2,:) = ones(5,5)*4;
% N(:,3,:) = ones(5,5)*3;
% N(:,4,:) = ones(5,5)*2;
% N(:,5,:) = ones(5,5)*1;
% 
% N
% 
% B = N(2:6,2:6,2:6);
% /home/ray
% B(4,3,3) = 8.5;
% 
% Adjustment = TaylorMaxima (B)
% 
% 
% M= [  1     5     4     7     1
%      6     2     1     4     8
%      2     5     2     8     0
%      4     6     7     1     6
%      1     6     1     4     7 ];
%  M=M+20;
%  
%  Kernel /home/ray= [-1 0 1
%            -2 0 2
%            -1 0 1];
% D_M = zeros(5,5);       
% for i=2:4
%     for j=2:4
%         D_M(i,j) = sum(sum(M(i-1:i+1,j-1:j+1).*Kernel));
%     end
% end
%         
% D_M
% DD_M = sum(sum(D_M(2:4,2:4).*Kernel))      
%         

% FilterSize = 7; SamplePerScale = 3; InitialSigma = 1.6; Nominal_Sigma = .5; StartingOctave = -1; LastOctave = 5;Ther = 0.0053;
% 
% o =StartingOctave:1:LastOctave
% s = (-1:1:SamplePerScale+1) 
% 
% for i=1:(LastOctave-StartingOctave) 
%     for j = 1:(SamplePerScale+3)
% 
%         Sigma(i,j) = InitialSigma * 2^(o(i)+s(j)/SamplePerScale);
%     end
% end

% Sigma =1;        
% 
% 
% 
% 
% x=-3:.1:3;
% y=-3:.1:3;
% 
% for i=1:61
%     for j=1:61
%             N(i,j) = ((2*pi*(Sigma^2))^(-.5))*exp(-(x(i)^2 + y(j)^2)/(2*(Sigma^2)));
%     end
% end
% mesh(N)
% 
% H = hess(N);
% 
% 

% Mu = [0 0];
% Sigma = [1 .3; .3 1];
% x1 = -3:.2:3; x2 = -3:.2:3;
% [X1,X2] = meshgrid(x1,x2);
% F = mvnpdf([X1(:) X2(:)],Mu,Sigma);
% F = reshape(F,length(x2),length(x1));
% surf(x1,x2,F);
% caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
% axis([-3 3 -3 3 0 .4])
% xlabel('x1'); ylabel('x2'); zlabel('Probability Density');




%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Plane equation in 2D

% x1 = -3:.2:3; x2 = -3:.2:3;
% [X1,X2] = meshgrid(x1,x2);
% 
% H = [1 0; 0 1];
% J = [1 ;1];
% 
% %F = (X'*H)*X;
% 
% F = [X1(:) X2(:)]*J;
% 
% F = reshape(F,length(x2),length(x1));
% surf(x1, x2, F);
% xlabel('x1'); ylabel('x2'); zlabel('F');

% 
% Signal = [
% 
%    180
%     90
%     84
%    133
%    154
%    181
%    170
%    100
%    114
%    184
% ];
% plot(Signal)
% hold on 
% 
% Dsignal = -(Signal(1:end-1) - Signal(2:end));
% plot(Dsignal,'g')
% hold on
% 
% DDsignal = Signal(1:end-2) - 2*Signal(2:end-1) + Signal(3:end)
% %DDsignal(2:end+1) = DDsignal(1:end)
% 
% hold on
% plot(DDsignal, 'r')
% grid on
% 
% NewSignal = zeros(100,1);
% for i = 1:8
%     NewSignal(i*10) = Signal(i);
%     for j=1:4
%         NewSignal(i*10+j) = Signal(i) + Dsignal(i)*(.1*j) + .5*DDsignal(i)*((.1*j)^2);
%     end
%     
% end
% figure(2)
% plot(NewSignal)
% grid on
% 
% figure (3)
% ReverseSignal=flipud(Signal);
% plot(ReverseSignal)
% hold on 
% DReverseSignal = -(ReverseSignal(1:end-1) - ReverseSignal(2:end));
% plot(DReverseSignal,'g')
% hold on
% 
% DDReverseSignal = ReverseSignal(1:end-2) - 2*ReverseSignal(2:end-1) + ReverseSignal(3:end)
% DDReverseSignal(2:end+1) = DDReverseSignal(1:end)
% 
% hold on
% plot(DDReverseSignal, 'r')
% grid on
% 
% NewReverseSignal = zeros(100,1);
% for i = 1:8
%     NewReverseSignal(i*10) = ReverseSignal(i);
%     for j=1:4
%         NewReverseSignal(i*10+j) = ReverseSignal(i) + DReverseSignal(i)*(.1*j) + .5*DDReverseSignal(i)*((.1*j)^2);
%     end
%     
% end
% figure(4)
% SecondHalf = zeros(100,1);
% temp = (flipud(NewReverseSignal));
% SecondHalf(10:end) = temp(1:91);
% plot(SecondHalf);
% hold on
% plot(NewSignal)
% grid on
% 
% figure(5)
% Result = (SecondHalf + NewSignal);
% for i =25:5:89
%     Result(i) = (Result(i-1)+Result(i+1))/2;
% end
% plot(Result(20:90))
% 
% hold on
% 
% 
% ReSignal = zeros(100,1);
% ReSignal(10:10:99) = Signal(2:end);
% for i = 1:9
%     for j=1:9
%         ReSignal(i*10+j) = (1-j*.1)*ReSignal(i*10) + (j*.1)*ReSignal((i+1)*10);
%     end
% end
% plot(ReSignal(10:90))
% grid on
