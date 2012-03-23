function a = convertFBdata2Native()

addpath('/Applications/freesurfer/matlab/');
input_file = '/Users/yrjo/yrjo_work/mindboggle2/2011_HBM/data/CUMC12/Heads/m';

for subjectId = 1:12
    
    
    filename_data = [input_file num2str(subjectId) '.nii.gz'];
    
    % Get transformation matrix for all subjects (they differ)
    command = ['/Applications/freesurfer/bin/mri_info --cras ' filename_data];
    [status, result] = system(command);
    vec = str2num(result);
    crasT = [1 0 0 vec(1); 0 1 0 vec(2); 0 0 1 vec(3); 0 0 0 1]
    
    % Load the image data to get voxel dimensions
    struct_nii = load_nifti(filename_data);
    data = struct_nii.vol;
    
    voxelT = struct_nii.vox2ras;
    
    lh_fundiInVoxCRS = inv(voxelT)*crasT*[lh_fundi ones(size(lh_fundi,1),1)]';
    lh_fundiInVoxCRS = int32(lh_fundiInVoxCRS(1:3,:)');
    
    
    % Convert the 3D coordinates of an existing VTK file to native volume space
    filename_fundi_vtk_lh = ['../../experiments/' subject '/features/lh.fundi.gang.li.native.vtk']
    filename_fundi_vtk_rh = ['../../experiments/' subject '/features/rh.fundi.gang.li.native.vtk']
    
    %LH
    myLine = lh_C{5}
    numPoints = str2num(myLine(8:end-6))
    
    fid = fopen(filename_fundi_vtk_lh, 'wt');
    jj = 1
    for ii = 1:length(lh_C)
        if(ii > 5 && ii <= (5+numPoints))
            fprintf(fid,'%f %f %f\n', lh_fundiInVoxCRS(jj,1),lh_fundiInVoxCRS(jj,2),lh_fundiInVoxCRS(jj,3));
            jj = jj+1;
        else
            fprintf(fid, lh_C{ii});
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    %RH
    myLine = rh_C{5}
    numPoints = str2num(myLine(8:end-6))
    
    fid = fopen(filename_fundi_vtk_rh, 'wt');
    jj = 1
    for ii = 1:length(rh_C)
        if(ii > 5 && ii <= (5+numPoints))
            fprintf(fid,'%f %f %f\n', rh_fundiInVoxCRS(jj,1),rh_fundiInVoxCRS(jj,2),rh_fundiInVoxCRS(jj,3));
            jj = jj+1;
        else
            fprintf(fid, rh_C{ii});
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    
