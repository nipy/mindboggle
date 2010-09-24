function Volume = VolumeShow(InputVolume)


[M N R] = size(InputVolume);

Min = min(InputVolume(:));
Max = max(InputVolume(:));


figure(round(rand*100));
subplot (3, 3, 1); imshow(imrotate(squeeze(InputVolume(round(M/3), :, :)),90),[Min Max]);
title('Sagittal View')
xlabel('Y')
ylabel('Z')

subplot (3, 3, 2); imshow(imrotate(squeeze(InputVolume(:, round(N/3), :)),90),[Min Max]);
title('Coronal View')
xlabel('X')
ylabel('Z')

subplot (3, 3, 3); imshow(imrotate(squeeze(InputVolume(:, :, round(R/3))),90),[Min Max]);
title('Transverse View')
xlabel('X')
ylabel('Y')

subplot (3, 3, 4); imshow(imrotate(squeeze(InputVolume(round(M/2), :, :)),90),[Min Max]);
xlabel('Y')
ylabel('Z')
subplot (3, 3, 5); imshow(imrotate(squeeze(InputVolume(:, round(N/2), :)),90),[Min Max]);
xlabel('X')
ylabel('Z')
subplot (3, 3, 6); imshow(imrotate(squeeze(InputVolume(:, :, round(R/2))),90),[Min Max]);
xlabel('X')
ylabel('Y')
subplot (3, 3, 7); imshow(imrotate(squeeze(InputVolume(round(2*M/3), :, :)),90),[Min Max]);
xlabel('Y')
ylabel('Z')
subplot (3, 3, 8); imshow(imrotate(squeeze(InputVolume(:, round(2*N/3), :)),90),[Min Max]);
xlabel('X')
ylabel('Z')
subplot (3, 3, 9); imshow(imrotate(squeeze(InputVolume(:, :, round(2*R/3))),90),[Min Max]);
xlabel('X')
ylabel('Y')

 
 
 % figure(2)
% imshow(squeeze(RotatedVolume(:,:,20)), [min(RotatedVolume(:)) max(RotatedVolume(:))])
% 
% figure(1)
% imshow(squeeze(InputVolume(:,50,:)),[min(inputVolume(:)) max(InputVolume(:))])
% figure(2)
% imshow(squeeze(RotatedVolume(:,50,:)), [min(RotatedVolume(:)) max(RotatedVolume(:))])
% 
% figure(1)
% imshow(squeeze(InputVolume(50,:,:)),[min(inputVolume(:)) max(InputVolume(:))])
% figure(2)
% imshow(squeeze(RotatedVolume(50,:,:)), [min(RotatedVolume(:)) max(RotatedVolume(:))])