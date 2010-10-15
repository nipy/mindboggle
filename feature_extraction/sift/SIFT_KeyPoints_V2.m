function [KeyPoints FilteredImage] = SIFT_KeyPoints_V2(Input_Image,SamplePerScale, InitialSigma, Nominal_Sigma, StartingOctave, LastOctave, Ther, SearchArea, EdgeDetect)

%%%%%%%%%%%%%%%%%%%%%%% Introduction
% This function obtains the keypoints of a 2D or 3D image for given
% parameters. 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization
TotalCount = 0;
RobustCount = 0;
Octave_Count = 0;
RobustOctaveCount = 0;
IgnordKeyPoint = 0;
Edge = 0;
KeyPoint = [0 0 0 0];
IgnoredCount = 0;
KeypointTemp = 0;
KeypointDiscardedTemp = 0;
Adjustment = zeros(1,4);

%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer the image to the mathematical
%%%%%%%%%%%%%%%%%%%%%%%%%% coordinate system


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
            FilteredImage(s).img= imfilter(Norm_X, h(s).kernel, 'replicate');
        elseif (s==1 && Octave~=1)
            FilteredImage(s).img= NextOctaveInput(1:2:end, 1:2:end, 1:2:end);	
        else
            FilteredImage(s).img = imfilter(FilteredImage(s-1).img, h(s).kernel, 'replicate');
        end
    end
    NextOctaveInput = FilteredImage(SamplePerScale+1).img;

    [m n r] = size(NextOctaveInput);
    
    % Compute difference of Gaussians
    for s = 1:(SamplePerScale+2)
        if r>3
            Diff_of_Gauss(:,:,:,s) = FilteredImage(s+1).img - FilteredImage(s).img;
        else
            Diff_of_Gauss(:,:,s) = FilteredImage(s+1).img - FilteredImage(s).img;
        end
    end
    
   
    
    iter = (round((SearchArea-1)/2));

    % Search for local extrimum 
if r>3	
    
    for i=(iter+2):m-(iter+1)
         for j=(iter+2):n-(iter+1)
             for k=(iter+2):r-(iter+1)
                 for s=1:SamplePerScale
                     if abs(Diff_of_Gauss(i,j,k,s+1)) > (0.8*Ther)
                        Sub_DoG = Diff_of_Gauss(i-iter:i+iter, j-iter:j+iter, k-iter:k+iter, s:s+2);
                        if (Sub_DoG(2,2,2,2) == max(Sub_DoG(:)) || Sub_DoG(2,2,2,2) == min(Sub_DoG(:)))
                            TotalCount = TotalCount + 1;
                            KeyPoint_Candidate(TotalCount, 1:4) = [i j k s]; 
                        end
                     end
                 end
             end
         end
    end
    
    %refine the key points
    for i = (Octave_Count+1):TotalCount
  
        SetFlag = EdgeDetect;
        Iteration = 0;
        NewMaxima = Diff_of_Gauss(KeyPoint_Candidate(i,1), KeyPoint_Candidate(i,2), KeyPoint_Candidate(i,3), KeyPoint_Candidate(i,4)+1);
    
        while (SetFlag == true && Iteration < 3)
        
            temp =  Diff_of_Gauss(KeyPoint_Candidate(i,1)-1:KeyPoint_Candidate(i,1)+1, ... 
                              KeyPoint_Candidate(i,2)-1:KeyPoint_Candidate(i,2)+1, ...
                              KeyPoint_Candidate(i,3)-1:KeyPoint_Candidate(i,3)+1, ...
                              KeyPoint_Candidate(i,4):KeyPoint_Candidate(i,4)+2);
   
            [Adjustment NewMaxima Edge]= TaylorExtremum(temp);
            SetFlag = false;
            Iteration = Iteration + 1;
        
            if abs(Adjustment(1))>.6
                KeyPoint_Candidate(i,1) = KeyPoint_Candidate(i,1) + 1*sign(Adjustment(1));
                SetFlag = true;
            end
            if abs(Adjustment(2))>.6
                KeyPoint_Candidate(i,2) = KeyPoint_Candidate(i,2) + 1*sign(Adjustment(2));
                SetFlag = true;
            end
            if abs(Adjustment(3))>.6
                KeyPoint_Candidate(i,3) = KeyPoint_Candidate(i,3) + 1*sign(Adjustment(3));
                SetFlag = true;
            end
            if abs(Adjustment(4))>.6
                if ((KeyPoint_Candidate(i,4)==1 && sign(Adjustment(4)) == -1) || (KeyPoint_Candidate(i,4)== SamplePerScale && sign(Adjustment(4)) == 1)) 
                    SetFlag = false;
                else
                    KeyPoint_Candidate(i,4) = KeyPoint_Candidate(i,4) + 1*sign(Adjustment(4));
                    SetFlag = true;
                end
            end
        end
        if ((Edge == false) && (abs(NewMaxima) > Ther))
            RobustCount = RobustCount + 1;
            for d=1:4                    %final check for the adjustment to be in the range
                if abs(Adjustment(d))>.5 
                    Adjustment(d)=sign(Adjustment(d))*.5;
                end
            end
      
            KeyPoints(RobustCount,1:4) = [((KeyPoint_Candidate(i,1)-1+Adjustment(2))* 2^(Octave-1+StartingOctave)), ...
                     ((KeyPoint_Candidate(i,2)-1+Adjustment(1))* 2^(Octave-1+StartingOctave)), ...
                     ((KeyPoint_Candidate(i,3)-1+Adjustment(3))* 2^(Octave-1+StartingOctave)), ...
                     (InitialSigma * 2^(Octave-1+StartingOctave+((KeyPoint_Candidate(i,4)+Adjustment(4)-1)/SamplePerScale))) ];
        end
    end
    fprintf('\n%G Candidate founded in Octave: %G', (TotalCount-Octave_Count), Octave-1+StartingOctave);
    fprintf('\n%G KeyPoints Ignored for low contrast in Octave: %G', ((TotalCount-Octave_Count)-(RobustCount-RobustOctaveCount)), Octave-1+StartingOctave);
    fprintf('\n%G KeyPoints founded in Octave: %G\n', ((RobustCount-RobustOctaveCount)), Octave-1+StartingOctave);
    
else
    
    for i=(iter+2):m-(iter+1)
        for j=(iter+2):n-(iter+1)
             for s=1:SamplePerScale
                if abs(Diff_of_Gauss(i,j,s+1)) > (0.8*Ther)
                    Sub_DoG = Diff_of_Gauss(i-iter:i+iter, j-iter:j+iter, s:s+2);
                    if (Sub_DoG(2,2,2) == max(Sub_DoG(:)) || Sub_DoG(2,2,2) == min(Sub_DoG(:)))
                         TotalCount = TotalCount + 1;
                         KeyPoint_Candidate(TotalCount, 1:3) = [i j s]; 
                    end
                end
             end
         end
    end
    
    %refine the key points
    for i = (Octave_Count+1):TotalCount
  
        SetFlag = EdgeDetect;
        Iteration = 0;
        NewMaxima = Diff_of_Gauss(KeyPoint_Candidate(i,1), KeyPoint_Candidate(i,2), KeyPoint_Candidate(i,3)+1);
    
        while (SetFlag == true && Iteration < 3)
        
            temp =  Diff_of_Gauss(KeyPoint_Candidate(i,1)-1:KeyPoint_Candidate(i,1)+1, ... 
                              KeyPoint_Candidate(i,2)-1:KeyPoint_Candidate(i,2)+1, ...
                              KeyPoint_Candidate(i,3):KeyPoint_Candidate(i,3)+2);
   
            [Adjustment NewMaxima Edge]= TaylorExtremum(temp);
            SetFlag = false;
            Iteration = Iteration + 1;
        
            if abs(Adjustment(1))>.6
                KeyPoint_Candidate(i,1) = KeyPoint_Candidate(i,1) + 1*sign(Adjustment(1));
                SetFlag = true;
            end
            if abs(Adjustment(2))>.6
                KeyPoint_Candidate(i,2) = KeyPoint_Candidate(i,2) + 1*sign(Adjustment(2));
                SetFlag = true;
            end
            if abs(Adjustment(3))>.6
                if ((KeyPoint_Candidate(i,3)==1 && sign(Adjustment(3)) == -1) || (KeyPoint_Candidate(i,3)== SamplePerScale && sign(Adjustment(3)) == 1)) 
                    SetFlag = false;
                else
                    KeyPoint_Candidate(i,3) = KeyPoint_Candidate(i,3) + 1*sign(Adjustment(3));
                    SetFlag = true;
                end
            end
        
        
        end
        if ((Edge == false) && (abs(NewMaxima) > Ther))
            RobustCount = RobustCount + 1;
            
            for d=1:3                    %final check for the adjustment to be in the range
                if abs(Adjustment(d))>.5 
                    Adjustment(d)=sign(Adjustment(d))*.5;
                end
            end
                    
            KeyPoints(RobustCount,1:3) = [((KeyPoint_Candidate(i,2)-1+Adjustment(2))* 2^(Octave-1+StartingOctave)), ...
                     ((KeyPoint_Candidate(i,1)-1+Adjustment(1))* 2^(Octave-1+StartingOctave)), ...
                     (InitialSigma * 2^(Octave-1+StartingOctave+((KeyPoint_Candidate(i,3)+Adjustment(3)-1)/SamplePerScale))) ];
        end
    end
    
    fprintf('\n%G Candidate founded in Octave: %G', (TotalCount-Octave_Count), Octave-1+StartingOctave);
    fprintf('\n%G KeyPoints Ignored for low contrast in Octave: %G', ((TotalCount-Octave_Count)-(RobustCount-RobustOctaveCount)), Octave-1+StartingOctave);
    fprintf('\n%G KeyPoints founded in Octave: %G\n', ((RobustCount-RobustOctaveCount)), Octave-1+StartingOctave);
end        
    Octave_Count = TotalCount;
    RobustOctaveCount = RobustCount;
    clear Diff_of_Gauss;
end


fprintf('\n\nFinal valid Keypoints are: %G\n\n', RobustCount);
