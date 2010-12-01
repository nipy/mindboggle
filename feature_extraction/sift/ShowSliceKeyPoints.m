function [KeyPoints_inSlice KeyPoints_inSlice_1] = ShowSliceKeyPoints(InputImage, KeyPoints1, RotatedPoint_10, Plane, SliceNumber)
Norm_X = InputImage;

Count = 0;

if Plane == 1

    for i =1:length(KeyPoints1)
        if (KeyPoints1(i,1) > (SliceNumber-1.5) && KeyPoints1(i,1) < (SliceNumber+1.5))
            Count = Count +1;
            KeyPoints_inSlice(Count,1) = KeyPoints1(i,2);
            KeyPoints_inSlice(Count,2) = KeyPoints1(i,3);

        end
    end
    Count = 0;
    for i =1:length(RotatedPoint_10)
        if (RotatedPoint_10(i,1) > (SliceNumber-1.5) && RotatedPoint_10(i,1) < (SliceNumber+1.5))
            Count = Count +1;
            KeyPoints_inSlice_1(Count,1) = RotatedPoint_10(i,2);
            KeyPoints_inSlice_1(Count,2) = RotatedPoint_10(i,3);

        end
    end
   imshow((imrotate(fliplr(squeeze(Norm_X(SliceNumber,:,:))),90)),[min(Norm_X(:)) max(Norm_X(:))]);
    hold on
 
elseif Plane == 2
    for i =1:length(KeyPoints1)
        if (KeyPoints1(i,2) > (SliceNumber-1.5) && KeyPoints1(i,2) < (SliceNumber+1.5))
            Count = Count +1;
            KeyPoints_inSlice(Count,1) = KeyPoints1(i,1);
            KeyPoints_inSlice(Count,2) = KeyPoints1(i,3);

        end
    end
    Count = 0;
    for i =1:length(RotatedPoint_10)
        if (RotatedPoint_10(i,2) > (SliceNumber-1.5) && RotatedPoint_10(i,2) < (SliceNumber+1.5))
            Count = Count +1;
            KeyPoints_inSlice_1(Count,1) = RotatedPoint_10(i,1);
            KeyPoints_inSlice_1(Count,2) = RotatedPoint_10(i,3);

        end
    end
   imshow((imrotate(fliplr(squeeze(Norm_X(:,SliceNumber,:))),90)),[min(Norm_X(:)) max(Norm_X(:))]);
    hold on

elseif Plane == 3
    
   for i =1:length(KeyPoints1)
        if (KeyPoints1(i,3) > (SliceNumber-1.5) && KeyPoints1(i,3) < (SliceNumber+1.5))
            Count = Count +1;
            KeyPoints_inSlice(Count,1) = KeyPoints1(i,1);
            KeyPoints_inSlice(Count,2) = KeyPoints1(i,2);

        end
    end
    Count = 0;
    for i =1:length(RotatedPoint_10)
        if (RotatedPoint_10(i,3) > (SliceNumber-1.5) && RotatedPoint_10(i,3) < (SliceNumber+1.5))
            Count = Count +1;
            KeyPoints_inSlice_1(Count,1) = RotatedPoint_10(i,1);
            KeyPoints_inSlice_1(Count,2) = RotatedPoint_10(i,2);

        end
    end
   imshow((imrotate(fliplr(squeeze(Norm_X(:,:,SliceNumber))),90)),[min(Norm_X(:)) max(Norm_X(:))]);
    hold on
    
end

    scatter(KeyPoints_inSlice(:,1), KeyPoints_inSlice(:,2),'o','y')
    scatter(KeyPoints_inSlice_1(:,1), KeyPoints_inSlice_1(:,2),'+','r')


