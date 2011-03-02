function [KeyPoint IgnordKeyPoint] = SIFT_KeyPoints(Input_Image1,SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect)

%%%%%%%%%%%%%%%%%%%%%%% Introduction
% This function obtains the keypoints of a 2D or 3D image for given
% parameters. 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization
Count = 0;
IgnordKeyPoint = 0;
KeyPoint = [0 0 0 0];
IgnoredCount = 0;
KeypointTemp = 0;
KeypointDiscardedTemp = 0;
Adjustment = zeros(1,4);

%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer the image to the mathematical
%%%%%%%%%%%%%%%%%%%%%%%%%% coordinate system

[m n r] = size(Input_Image1);
for i =1:r
    Input_Image(:,:,i) = flipud(imrotate(Input_Image1(:,:,i),90));
end
clear Input_Image1;
[m n r] = size(Input_Image);



%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the Gaussian Sigma for every scale space
%%%%%%%%%%%%%%%%%%%%%%%%%% in the process
[Sigma dSigma] = GaussianSigma(InitialSigma, Nominal_Sigma, SamplePerScale, StartingOctave, LastOctave);


% Double the size of the image for octave -1
if StartingOctave == -1

	Norm_X = DoubleSize(Input_Image);
else

	Norm_X = double(Input_Image);
end

%Starting the actual loop for finding keypoints
for Octave = 1:LastOctave
    
    %Building Gaussian Kernel for every scale 
    for i = 1:(SamplePerScale+3)
        h(i).kernel = GaussianKernel(dSigma(Octave,i), r);   
    end
       
    %Obtain the blurred images 
    for s = 1:(SamplePerScale+3)
        if (Octave==1 && s==1)
            FilteredImage(s).img= imfilter(Norm_X, h(s).kernel);
        elseif (s==1 && Octave~=1)
            FilteredImage(s).img= NextOctaveInput(1:2:end, 1:2:end, 1:2:end);	
        else
            FilteredImage(s).img = imfilter(FilteredImage(s-1).img, h(s).kernel);
        end
    end
    NextOctaveInput = FilteredImage(SamplePerScale+1).img;

    
    
    % Compute difference of Gaussians
    for s = 1:(SamplePerScale+2)
        DoG(s).img = FilteredImage(s+1).img - FilteredImage(s).img;
    end
     
   
    [m n r] = size(DoG(1).img);
    iter = (round((SearchArea-1)/2));

    % Search for local extrimum 
    if r>3	
    for i=(iter+2):m-(iter+1)
         for j=(iter+2):n-(iter+1)
             for k=(iter+2):r-(iter+1)
                 for s=1:SamplePerScale
                     
                     Sub_DoG(:,:,:,1) = DoG(s).img(i-iter:i+iter,j-iter:j+iter,k-iter:k+iter);
                     Sub_DoG(:,:,:,2) = DoG(s+1).img(i-iter:i+iter,j-iter:j+iter,k-iter:k+iter);
                     Sub_DoG(:,:,:,3) = DoG(s+2).img(i-iter:i+iter,j-iter:j+iter,k-iter:k+iter);
                                              
                     if (abs(DoG(s+1).img(i,j,k)) > Ther)                   
                        if (DoG(s+1).img(i,j,k) == max(Sub_DoG(:)) || DoG(s+1).img(i,j,k) == min(Sub_DoG(:)))
                            
                            index1 = i* 2^(Octave-1+StartingOctave);
                            index2 = j* 2^(Octave-1+StartingOctave);
                            index3 = k* 2^(Octave-1+StartingOctave);
                            
                            Scale = InitialSigma * 2^(Octave-1+StartingOctave+((s-1)/SamplePerScale));
                            
                            temp = zeros(5,5,5,3);
                            temp(:,:,:,1) = DoG(s).img(i-2:i+2,j-2:j+2,k-2:k+2);
                            temp(:,:,:,2) = DoG(s+1).img(i-2:i+2,j-2:j+2,k-2:k+2);
                            temp(:,:,:,3) = DoG(s+2).img(i-2:i+2,j-2:j+2,k-2:k+2);
                            
                            [Adjustment NewMaxima Edge]= TaylorMaxima (temp);
                            
                            if EdgeDetect
                                if Edge
                                    Count = Count +1;
                                    KeyPoint (Count, 1:4) = [index1+Adjustment(1) index2+Adjustment(2) index3+Adjustment(3) Scale+Adjustment(4)];
                                else
                                    %fprintf('\n this keypoint is ignored for being along the edge: (%G, %G, %G)',index1, index2, index3);
                                    IgnoredCount = IgnoredCount+1;
                                    IgnordKeyPoint (IgnoredCount, 1:4) = [index1+Adjustment(1) index2+Adjustment(2) index3+Adjustment(3) Scale+Adjustment(4)];
                                end
                            else
                                Count = Count +1;
                                KeyPoint (Count, 1:4) = [index1+Adjustment(1) index2+Adjustment(2) index3+Adjustment(3) Scale+Adjustment(4)];
                            
                            end
                            
                        end
                     end
                 end
             end
         end
    end
    
    
    else
    for i=(iter+2):m-(iter+1)
         for j=(iter+2):n-(iter+1)
              
             for s=1:SamplePerScale
                                          
                     Sub_DoG(:,:,1) = DoG(s).img(i-iter:i+iter, j-iter:j+iter);
                     Sub_DoG(:,:,2) = DoG(s+1).img(i-iter:i+iter, j-iter:j+iter);
                     Sub_DoG(:,:,3) = DoG(s+2).img(i-iter:i+iter, j-iter:j+iter);
                                              
                     if (abs(DoG(s+1).img(i,j)) > Ther)                   
                        if (DoG(s+1).img(i,j) == max(Sub_DoG(:)) || DoG(s+1).img(i,j) == min(Sub_DoG(:)))
                            
                            index1 = i* 2^(Octave-1+StartingOctave);
                            index2 = j* 2^(Octave-1+StartingOctave);
                            
                            Scale = InitialSigma * 2^(Octave-1+StartingOctave+((s-1)/SamplePerScale));
                            temp = zeros(5,5,3);
                            temp(:,:,1) = DoG(s).img(i-2:i+2,j-2:j+2);
                            temp(:,:,2) = DoG(s+1).img(i-2:i+2,j-2:j+2);
                            temp(:,:,3) = DoG(s+2).img(i-2:i+2,j-2:j+2);
                            
                            [Adjustment NewMaxima Edge]= TaylorMaxima (temp);
                            
                            if EdgeDetect 
                                if Edge
                                    Count = Count +1;
                                    KeyPoint (Count, 1:3) = [index1+Adjustment(1) index2+Adjustment(2) Scale+Adjustment(3)];
                                else
                                    IgnoredCount = IgnoredCount+1;
                                    IgnordKeyPoint (IgnoredCount, 1:3) = [index1+Adjustment(1) index2+Adjustment(2) Scale+Adjustment(3)];
                                    %fprintf('\n this keypoint is ignored for being along the edge: (%G, %G)',index1, index2);
                                end
                            else
                                Count = Count +1;
                                KeyPoint (Count, 1:3) = [index1+Adjustment(1) index2+Adjustment(2) Scale+Adjustment(3)];    
                            end
                        end
                     end
             end
         end
    end
    end
fprintf('\n%G KeyPoints founded in Octave: %G', (Count-KeypointTemp+IgnoredCount-KeypointDiscardedTemp), Octave-1+StartingOctave);
fprintf('\n%G KeyPoints discarded in Octave: %G', IgnoredCount-KeypointDiscardedTemp, Octave-1+StartingOctave);
fprintf('\n%G Valid KeyPoints Octave: %G \n', Count-KeypointTemp, Octave-1+StartingOctave);

KeypointTemp = Count;
KeypointDiscardedTemp = IgnoredCount;
end
fprintf('\nTotal Keypoints are: %G', IgnoredCount+Count);
fprintf('\nTotal discarded edge Keypoints are: %G', IgnoredCount);
fprintf('\nFinal valid Keypoints are: %G\n\n', Count);


