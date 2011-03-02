function DJ1 = DoubleSize(I)

%  image1 = load_untouch_nii('brain2_1.nii');
%  I = double(image1.img);
%  Slice_X = X(:, :, 70);
%  I = double(Slice_X);
%  
 [M,N,R]=size(I) ;
% 
if R>1
%     for k=1:R
%         J(:,:,k) = imresize(I(:,:,k), [2*M 2*N],'bilinear');
%     end
%     for i=1:2*M
%         DJ(i,:,:) = imresize(squeeze(J(i,:,:)), [2*N 2*R], 'bilinear');
%     end
% else
%     for k=1:R
%         DJ(:,:,k) = imresize_old(I(:,:,k),2,'bilinear');
%     end
% end    
    
    DJ1 = zeros(2*M,2*N,2*R);
    DJ1(1:2:end,1:2:end,1:2:end)=I ;
    DJ1(2:2:end-1,2:2:end-1,2:2:end-1)= (I(1:end-1,1:end-1,1:end-1) + I(2:end,1:end-1,1:end-1) + ...
		I(1:end-1,2:end,1:end-1) + I(2:end,2:end,1:end-1) + I(1:end-1,1:end-1,2:end) + I(2:end,1:end-1,2:end) +...
        I(1:end-1,2:end,2:end) + I(2:end,2:end,2:end))/8;
    
    DJ1(2:2:end-1,2:2:end-1,1:2:end) = (I(1:end-1,1:end-1,1:end) + I(2:end,1:end-1,1:end) + ...
		I(1:end-1,2:end,1:end) + I(2:end,2:end,1:end))/4 ;
	DJ1(2:2:end-1,1:2:end,1:2:end) = 0.5*I(1:end-1,:,1:end) + 0.5*I(2:end,:,1:end) ;
	DJ1(1:2:end,2:2:end-1,1:2:end) = 0.5*I(:,1:end-1,1:end) + 0.5*I(:,2:end,1:end) ;
   
    DJ1(1:2:end,2:2:end-1,2:2:end-1) = (I(1:end,1:end-1,1:end-1) + I(1:end,2:end,1:end-1) + ...
		I(1:end,1:end-1,2:end) + I(1:end,2:end,2:end))/4 ;
	DJ1(1:2:end,2:2:end-1,1:2:end) = 0.5*I(1:end,1:end-1,:) + 0.5*I(1:end,2:end,:) ;
	DJ1(1:2:end,1:2:end,2:2:end-1) = 0.5*I(1:end,:,1:end-1) + 0.5*I(1:end,:,2:end) ;
    
    DJ1(2:2:end-1,1:2:end,2:2:end-1) = (I(1:end-1,1:end,1:end-1) + I(2:end,1:end,1:end-1) + ...
		I(1:end-1,1:end,2:end) + I(2:end,1:end,2:end))/4 ;
	DJ1(2:2:end-1,1:2:end,1:2:end) = 0.5*I(1:end-1,1:end,:) + 0.5*I(2:end,1:end,:) ;
	DJ1(1:2:end,1:2:end,2:2:end-1) = 0.5*I(:,1:end,1:end-1) + 0.5*I(:,1:end,2:end) ;
 
else    
    

	DJ1 = zeros(2*M,2*N) ;
	DJ1(1:2:end,1:2:end) = I ;
	DJ1(2:2:end-1,2:2:end-1) = 0.25*I(1:end-1,1:end-1) + 0.25*I(2:end,1:end-1) + ...
		0.25*I(1:end-1,2:end) + 0.25*I(2:end,2:end) ;
	DJ1(2:2:end-1,1:2:end) = 0.5*I(1:end-1,:) + 0.5*I(2:end,:) ;
	DJ1(1:2:end,2:2:end-1) = 0.5*I(:,1:end-1) + 0.5*I(:,2:end) ;
    
end    
%    sum(sum(abs(double(DJ) - double(DJ1))))
% imshow(DJ1(:,:,161), [min(DJ1(:)) max(DJ1(:))])
% imshow(squeeze(DJ1(254,:,:)), [min(DJ1(:)) max(DJ1(:))])

