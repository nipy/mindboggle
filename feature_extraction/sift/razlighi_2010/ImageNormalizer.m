function Result = ImageNormalizer(Image, NormalizationType, NgL)
 
[m n r] = size(Image);

Fix_Min = double(min(Image(:)));
Fix_Max = double(max(Image(:)));

if (NormalizationType == 1)

    %Min Max Normalizer normalizes the image between its minimum and
    %maximum value based on the number of the gray level (NgL)
    
    Fix_NormFactor = ((NgL-1)/(Fix_Max - Fix_Min));
    Result = uint8(round((double(Image) - Fix_Min) * Fix_NormFactor));
    fprintf('Input Image Normalized from its [Min Max] value to the range of [0 %G] \n', NgL-1);
    fprintf('Original Min Max Value were [%G %G]\n', Fix_Min, Fix_Max);
    
    
elseif (NormalizationType == 2)

    % Maximum Energy Normalizer normalizes the image to the range that it
    % contains 99.9 percent of its energy.
    
    
    Hist = zeros(Fix_Max+1, 1);
    for i=1:m
        for j=1:n
            for k=1:r
                 Hist(Image(i,j,k)+1) =  Hist(Image(i,j,k)+1) + 1;
             
            end
        end
    end

    TotalEnergy = sum(sum(Hist));

    EnergyMin = 0;
    EnergyMax = 0;
    
    for i= 1:1000
    
        EnergyMin = EnergyMin + Hist(Fix_Min+i);
    
        if (EnergyMin/TotalEnergy >=.001)
            New_Min = i-1
            break;
        end
    end

    for i= 1:1000000
    
        EnergyMax = EnergyMax + Hist(Fix_Max+2-i);
    
        if (EnergyMax/TotalEnergy >=.001)
            New_Max = Fix_Max+2-i
            break;
        
        end
        
    end

    %Hist(New_Max-100:New_Max+100)    
    Fix_NormFactor = ((NgL-1)/(New_Max - New_Min));

    Result = uint8(round((double(Image) - New_Min) * Fix_NormFactor));  
    fprintf('Input Image Normalized to the range that contains 99.9 persent of its energy\n');
    fprintf('Original Min Max Value were [%G %G]\n', Fix_Min, Fix_Max);
    fprintf('New Min Max Value is [%G %G] which is again normalized to range [0 %d]\n', New_Min, New_Max, NgL-1);

   
    % Normalized Image to the range of [0 1) return as double based on maximin and minimum image intensity value
      
    
elseif (NormalizationType == 3)
    Fix_NormFactor = ((1)/(Fix_Max - Fix_Min));
    Result = (double(Image) - Fix_Min) * Fix_NormFactor;
    fprintf('Input Image Normalized from its range of [%G %G] to the range of [0 1) \n', Fix_Min, Fix_Max);
 
    
    
end

    