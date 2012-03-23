function a = outputPits(lh_pits,subjectId,side)

% Read nifti volume to get its transformation to native space
%subjectId = 1;
subject = ['CUMC12_' num2str(subjectId)]
filename_data = '';
filename_label = '';

%filename_data = ['../../data/CUMC12/Heads/m' num2str(subjectId) '.nii.gz'];
filename_data = ['/Users/yrjo/Public/2011_HBM/CUMC12/Heads/m' num2str(subjectId) '.nii'];
% Get transformation matrix for all subjects (they differ)
command = ['/Applications/freesurfer/bin/mri_info --cras ' filename_data];
[status, result] = system(command);

vec = str2num(result);
crasT = [1 0 0 vec(1); 0 1 0 vec(2); 0 0 1 vec(3); 0 0 0 1]

% Load the image data to get voxel dimensions
struct_nii = load_nifti(filename_data);
data = struct_nii.vol;

%struct_nii = load_nifti(filename_label);
%label = struct_nii.vol;

voxelT = struct_nii.vox2ras;



% Convert from freesurfer tkRas space to native volume space
lh_pitsInVoxCRS = inv(voxelT)*crasT*[lh_pits ones(size(lh_pits,1),1)]';
lh_pitsInVoxCRS = int32(lh_pitsInVoxCRS(1:3,:)');

% Save nifti volume for quality check
mask = zeros(size(data),'uint8');
pSize = 0;
% LH
for ii = 1:size(lh_pitsInVoxCRS,1)
    
    c=lh_pitsInVoxCRS(ii,1);
    r=lh_pitsInVoxCRS(ii,2);
    s=lh_pitsInVoxCRS(ii,3);
    
    mask(c-pSize:c+pSize,r-pSize:r+pSize,s-pSize:s+pSize) = 255;
end
% RH
% for ii = 1:size(rh_pitsInVoxCRS,1)
%     
%     c=rh_pitsInVoxCRS(ii,1);
%     r=rh_pitsInVoxCRS(ii,2);
%     s=rh_pitsInVoxCRS(ii,3);
%     
%     mask(c-pSize:c+pSize,r-pSize:r+pSize,s-pSize:s+pSize) = 255;
% end
struct_nii.vol = mask;
if (side == 1)
    filename_pits_mask = ['./experiments/lh_' subject 'pits_yrjo_hame_mask.nii.gz'];
else
    filename_pits_mask = ['./experiments/rh_' subject 'pits_yrjo_hame_mask.nii.gz'];
end
save_nifti(struct_nii, filename_pits_mask);

a = 0;
