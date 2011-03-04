function []=DTI_bootstrap_analysys(dti_data_file,bvals_file,bvecs_file,output_basename,options_in)
% Function for the bootstrap analysis of the DTI data.
% Compute the FA map and the FA precision using a bootstrap approach.
% INPUT
%  - dti_data_file: string with the name of the niftii file containing the DTI data
%  - bvals_file: string with the name of the b-value text file (FSL format)
%  - bvecs_file: string with the name of the gradient direction text file (FSL format)
%  - output_basename: string with the basename for the output images
%  - options_in: variable with the options. See the INITIALIZING OPTIONS section to see the available options.
%
% OUTPUT
% Everything will be saved both as niftii and matlab files. 
%  - FA map
%  - ADC map
%  - CI map (95% confidence intervals)
%  - CV map (FA CV)
%  - raw segmentation based only on DTI data.


%                                                      INITIALIZING OPTIONS
% -------------------------------------------------------------------------
% DEFAULT OPTIONS
% Bootstrap approach options
options.subsample_elements = 20;
options.number_of_subsamples = 100;

% Segmentation options
options.CSF_ADC_threshold = 0.0018;
options.WM_GM_FA_threshold = 0.2;

% DTI analysis options
options.parametersDTI.BackgroundTreshold=0; % Analyze only voxel with an S0 intensity grater than this threshold.
options.parametersDTI.WhiteMatterExtractionThreshold=0; % Store only results for voxel with an FA greater than this threshold
options.parametersDTI.textdisplay=false; % Display while computing the FA.

% Output options
[file_path,file_name]=fileparts(dti_data_file);
if strcmp(file_path,'')
    file_path=pwd;
end
options.output_path=file_path;
options.output_basename = [file_name '_bootstrap_analysis'];

% File in
options.file_dti_data = dti_data_file;

% USER DEFINED OPTIONS
if nargin>=4
    [file_path,file_name]=fileparts(output_basename);
    if not(strcmp(file_path,''))
        options.output_path=file_path;
    end
    if not(strcmp(file_name,''))
        options.output_basename = file_name;
    end
    if nargin>=5
        if isfield(options_in,'subsample_elements')
            options.subsample_elements=options_in.subsample_elements;
        end
        if isfield(options_in,'number_of_subsamples')
            options.number_of_subsamples=options_in.number_of_subsamples;
        end
        if isfield(options_in,'CSF_ADC_threshold')
            options.CSF_ADC_threshold=options_in.CSF_ADC_threshold;
        end
        if isfield(options_in,'WM_GM_FA_threshold')
            options.WM_GM_FA_threshold=options_in.WM_GM_FA_threshold;
        end
        if isfield(options_in,'parametersDTI.BackgroundTreshold')
            options.parametersDTI.BackgroundTreshold=options_in.parametersDTI.BackgroundTreshold;
        end
        if isfield(options_in,'parametersDTI.WhiteMatterExtractionThreshold')
            options.parametersDTI.WhiteMatterExtractionThreshold=options_in.parametersDTI.WhiteMatterExtractionThreshold;
        end
        if isfield(options_in,'parametersDTI.textdisplay')
            options.parametersDTI.textdisplay=options_in.parametersDTI.textdisplay;
        end
        if isfield(options_in,'MaskFile')
            options.MaskFile=options_in.MaskFile;
        end
    end
end
clear options_in

%                                                            OPTION SUMMARY
% -------------------------------------------------------------------------
display(' BOOTSTRAP DTI DATA ANALYSIS')
display(' ')
display(['DTI data file:   ',dti_data_file])
display(['B value file:    ',bvals_file])
display(['Gradient file:   ',bvecs_file])
if isfield(options,'MaskFile')
    display(['Mask file:   ',options.MaskFile])
else
    display('Mask:   not provided')
end
display(['Output folder:   ',options.output_path])
display(['Output basename: ',options.output_basename])
display(' ')
display(['Element n# in each subsample: ',num2str(options.subsample_elements)])
display(['Total subsamples: ',num2str(options.number_of_subsamples)])
display(' ')
display(['ADC threshold for segmentation: ',num2str(options.CSF_ADC_threshold)])
display(['FA threshold for segmentation:  ',num2str(options.WM_GM_FA_threshold)])
display(' ')

%                                                              DATA LOADING
% -------------------------------------------------------------------------
display('Loading data ...')
[DTI_data_unmasked,result_check]= DTI_data_loading(dti_data_file,bvals_file,bvecs_file,options);
if (~result_check)
    display(' WARNING: error while loading data...')
    exit
end

%                                                              MASKING DATA
% -------------------------------------------------------------------------
display('Masking data ...')
[DTI_mask,DTI_data,result_check] = mask_DTI_data(DTI_data_unmasked,options);
clear DTI_data_unmasked

%                                                 GENERATING THE SUBSAMPLES
% -------------------------------------------------------------------------
display('Generating the subsamples...')
[subsample_list,result_check]= generate_subsamples(DTI_data,options);
if (~result_check)
    display(' WARNING: error while generating the subsamples...')
    exit
end

%                                                        SUBSAMPLE ANALYSIS
% -------------------------------------------------------------------------
display('Subsample DTI analysis...')
%waitbar_obj = waitbar(0,'Subsample DTI analysis...');
for cont = 1:options.number_of_subsamples
    % Analysis of subsample cont-th
    subsample_DTI_data=DTI_data(subsample_list(cont,:));
    [subsample_ADC,subsample_FA,subsample_VectorF,subsample_DifT]=DTI(subsample_DTI_data,options.parametersDTI);
    subsamples_results(1,cont)=struct('ADC',subsample_ADC,'FA',subsample_FA);
    %waitbar(cont/options.number_of_subsamples,waitbar_obj)
end
%close(waitbar_obj)
clear subsample_ADC subsample_FA subsample_VectorF subsample_DifT cont subsample_DTI_data waitbar_obj

%                                                 BOOTSTRAP RESULT ANALYSIS
% -------------------------------------------------------------------------
display('Subsample result analysis...')
[bootstrap_results,result_check] = analyze_subsample_results(subsamples_results,DTI_mask,options);

%                                                           ROI DEFINITIONS
% -------------------------------------------------------------------------
display('Definition of white and grey matter ROIs...')
[ROI_maps,result_check]=find_ROIs(bootstrap_results,DTI_mask,options);

%                                                                MAP SAVING
% -------------------------------------------------------------------------
display('Saving the results...')
check_results = save_niftii_files(bootstrap_results,ROI_maps,options);
check_results = save_matlab_files(bootstrap_results,ROI_maps,options);


%                           ROI ANALYSIS AND RESULT PRESENTATION (optional)
% -------------------------------------------------------------------------
GM_indices = find(ROI_maps.GM);
WM_indices = find(ROI_maps.WM);
GM_CI_vett = bootstrap_results.CI_map(GM_indices);
WM_CI_vett = bootstrap_results.CI_map(WM_indices);
xVett = [0:0.001:1];
GM_yVett = hist(GM_CI_vett,xVett)./length(GM_CI_vett);
WM_yVett = hist(WM_CI_vett,xVett)./length(WM_CI_vett);
[GM_fmax_CI,GM_fmax_pos]=max(GM_yVett);
GM_CI_max=xVett(GM_fmax_pos);
[WM_fmax_CI,WM_fmax_pos]=max(WM_yVett);
WM_CI_max=xVett(WM_fmax_pos);

figure
plot(xVett,WM_yVett,'r-',xVett,GM_yVett,'b-')
title(strrep([options.output_basename,'- CI histogram'],'_',' '))
xlabel('95% CI length (FA)')
ylabel('CI frequency')
legend('WM','GM')
peak_text = {['WM peak: ',num2str(0.001*round(1000*WM_CI_max)),' - ',num2str(0.001*round(1000*WM_fmax_CI))];['GM peak: ',num2str(0.001*round(1000*GM_CI_max)),' - ',num2str(0.001*round(1000*GM_fmax_CI))]};
text(0.6,0.05,peak_text)
saveas(gcf,[options.output_path,filesep,options.output_basename,'_CIplot'],'png')
xlim([0 0.5])
text(0.3,0.05,peak_text)
saveas(gcf,[options.output_path,filesep,options.output_basename,'_CIplot_zoom'],'png')
clear xVett GM_indices GM_yVett GM_CI_vett WM_indices WM_yVett WM_CI_vett








% -------------------------------------------------------------------------
% |                           SAVE_MATLAB_FILES                           |
% -------------------------------------------------------------------------
function check_results = save_matlab_files(bootstrap_results,ROI_maps,options)
save_path = options.output_path;
basename  = options.output_basename;

FA_map  = bootstrap_results.FA_map;
ADC_map = bootstrap_results.ADC_map;
CV_map  = bootstrap_results.CV_map;
CI_map  = bootstrap_results.CI_map;


IDfile = [basename, '_maps.mat'];
save([save_path,filesep,IDfile],'IDfile','FA_map','ADC_map','CI_map','CV_map','ROI_maps')
check_results = true;







% -------------------------------------------------------------------------
% |                           SAVE_NIFTII_FILES                           |
% -------------------------------------------------------------------------
function check_results = save_niftii_files(bootstrap_results,ROI_maps,options)
save_path = options.output_path;
basename  = options.output_basename;

% load untouched Niftii image do get the and the data structure
niftii_file_obj = load_untouch_nii(options.file_dti_data);

% Maps
% FA
file_name = [basename,'_FA'];
[nR,nC,nS] = size(bootstrap_results.FA_map);
fa_niftii_file_obj                    = niftii_file_obj;
fa_niftii_file_obj.img                = single(bootstrap_results.FA_map);
fa_niftii_file_obj.fileprefix         = file_name;
fa_niftii_file_obj.hdr.dime.dim       = [3, nR, nC, nS, 1, 1, 1, 1];
fa_niftii_file_obj.hdr.dime.pixdim(5) = 1;
fa_niftii_file_obj.hdr.dime.datatype  = 16;
fa_niftii_file_obj.hdr.dime.bitpix    = 32;
save_untouch_nii(fa_niftii_file_obj,[save_path,filesep,file_name])
clear fa_niftii_file_obj file_name nR nC nS

% ADC
file_name = [basename,'_ADC'];
[nR,nC,nS] = size(bootstrap_results.ADC_map);
ADC_niftii_file_obj                    = niftii_file_obj;
ADC_niftii_file_obj.img                = single(bootstrap_results.ADC_map);
ADC_niftii_file_obj.fileprefix         = file_name;
ADC_niftii_file_obj.hdr.dime.dim       = [3, nR, nC, nS, 1, 1, 1, 1];
ADC_niftii_file_obj.hdr.dime.pixdim(5) = 1;
ADC_niftii_file_obj.hdr.dime.datatype  = 16;
ADC_niftii_file_obj.hdr.dime.bitpix    = 32;
save_untouch_nii(ADC_niftii_file_obj,[save_path,filesep,file_name])
clear ADC_niftii_file_obj file_name nR nC nS

% CI
file_name = [basename,'_CI'];
[nR,nC,nS] = size(bootstrap_results.CI_map);
CI_niftii_file_obj                    = niftii_file_obj;
CI_niftii_file_obj.img                = single(bootstrap_results.CI_map);
CI_niftii_file_obj.fileprefix         = file_name;
CI_niftii_file_obj.hdr.dime.dim       = [3, nR, nC, nS, 1, 1, 1, 1];
CI_niftii_file_obj.hdr.dime.pixdim(5) = 1;
CI_niftii_file_obj.hdr.dime.datatype  = 16;
CI_niftii_file_obj.hdr.dime.bitpix    = 32;
save_untouch_nii(CI_niftii_file_obj,[save_path,filesep,file_name])
clear CI_niftii_file_obj file_name nR nC nS

% CV
file_name = [basename,'_CV'];
[nR,nC,nS] = size(bootstrap_results.CI_map);
CV_niftii_file_obj                    = niftii_file_obj;
CV_niftii_file_obj.img                = single(bootstrap_results.CV_map);
CV_niftii_file_obj.fileprefix         = file_name;
CV_niftii_file_obj.hdr.dime.dim       = [3, nR, nC, nS, 1, 1, 1, 1];
CV_niftii_file_obj.hdr.dime.pixdim(5) = 1;
CV_niftii_file_obj.hdr.dime.datatype  = 16;
CV_niftii_file_obj.hdr.dime.bitpix    = 32;
save_untouch_nii(CV_niftii_file_obj,[save_path,filesep,file_name])
clear CV_niftii_file_obj file_name nR nC nS

% ROI GM
file_name = [basename,'_GM_roi'];
[nR,nC,nS] = size(bootstrap_results.CI_map);
GM_niftii_file_obj                    = niftii_file_obj;
GM_niftii_file_obj.img                = int16(ROI_maps.GM);
GM_niftii_file_obj.fileprefix         = file_name;
GM_niftii_file_obj.hdr.dime.dim       = [3, nR, nC, nS, 1, 1, 1, 1];
GM_niftii_file_obj.hdr.dime.pixdim(5) = 1;
GM_niftii_file_obj.hdr.dime.datatype  = 4;
GM_niftii_file_obj.hdr.dime.bitpix    = 16;
save_untouch_nii(GM_niftii_file_obj,[save_path,filesep,file_name])
clear GM_niftii_file_obj file_name nR nC nS

% ROI WM
file_name = [basename,'_WM_roi'];
[nR,nC,nS] = size(bootstrap_results.CI_map);
WM_niftii_file_obj                    = niftii_file_obj;
WM_niftii_file_obj.img                = int16(ROI_maps.WM);
WM_niftii_file_obj.fileprefix         = file_name;
WM_niftii_file_obj.hdr.dime.dim       = [3, nR, nC, nS, 1, 1, 1, 1];
WM_niftii_file_obj.hdr.dime.pixdim(5) = 1;
WM_niftii_file_obj.hdr.dime.datatype  = 4;
WM_niftii_file_obj.hdr.dime.bitpix    = 16;
save_untouch_nii(WM_niftii_file_obj,[save_path,filesep,file_name])
clear WM_niftii_file_obj file_name nR nC nS

% ROI CSF
file_name = [basename,'_CSF_roi'];
[nR,nC,nS] = size(bootstrap_results.CI_map);
CSF_niftii_file_obj                    = niftii_file_obj;
CSF_niftii_file_obj.img                = int16(ROI_maps.CSF);
CSF_niftii_file_obj.fileprefix         = file_name;
CSF_niftii_file_obj.hdr.dime.dim       = [3, nR, nC, nS, 1, 1, 1, 1];
CSF_niftii_file_obj.hdr.dime.pixdim(5) = 1;
CSF_niftii_file_obj.hdr.dime.datatype  = 4;
CSF_niftii_file_obj.hdr.dime.bitpix    = 16;
save_untouch_nii(CSF_niftii_file_obj,[save_path,filesep,file_name])
clear CSF_niftii_file_obj file_name nR nC nS

check_results=true;








% -------------------------------------------------------------------------
% |                               FIND_ROIS                               |
% -------------------------------------------------------------------------
function [ROI_maps,result_check]=find_ROIs(bootstrap_results,DTI_mask,options)
ADC_thr = options.CSF_ADC_threshold;
FA_thr = options.WM_GM_FA_threshold;
FA_map = bootstrap_results.FA_map;
ADC_map = bootstrap_results.ADC_map;
% Computing ROIs
% CSF
mask_indices = find(DTI_mask);
ROI = zeros(size(DTI_mask));
CSF_indices = intersect(mask_indices,find(ADC_map>ADC_thr));
ROI(CSF_indices)=1;
ROI_maps.CSF=ROI;
% WM
mask_indices = intersect(find(DTI_mask),find(ROI_maps.CSF==0));
ROI = zeros(size(DTI_mask));
WM_indices = intersect(mask_indices,find(FA_map>FA_thr));
ROI(WM_indices)=1;
ROI_maps.WM=ROI;
% GM
mask_indices = intersect(find(DTI_mask),find(ROI_maps.CSF==0));
ROI = zeros(size(DTI_mask));
GM_indices = intersect(mask_indices,find(FA_map<=FA_thr));
ROI(GM_indices)=1;
ROI_maps.GM=ROI;
result_check = true;








% -------------------------------------------------------------------------
% |                        ANALYZE_SUBSAMPLE_RESULTS                      |
% -------------------------------------------------------------------------
function [bootstrap_results,result_check] = analyze_subsample_results(subsamples_results,DTI_mask,options)
% Analyze the results from the subsample DTI analysis
[nR,nC,nS]=size(subsamples_results(1).FA);
nSubsamples=length(subsamples_results);
% Initializing result variables
mean_FA_map = zeros(nR,nC,nS);
mean_ADC_map = zeros(nR,nC,nS);
CI_map = zeros(nR,nC,nS);
CV_map = zeros(nR,nC,nS);
% Analyzing results
mask_indices = find(DTI_mask);
nVoxel = length(mask_indices);
%waitbar_obj = waitbar(0,'Analyzing subsample results ...');
for contI = 1:nVoxel
    voxelInd = mask_indices(contI);
    % Analysis of the voxel
    % Extracting voxel FA and ADC samples
    FA_vett = zeros(1,nSubsamples);
    ADC_vett = zeros(1,nSubsamples);
    for cont = 1:nSubsamples
        FA_vett(cont)=subsamples_results(cont).FA(voxelInd);
        ADC_vett(cont)=subsamples_results(cont).ADC(voxelInd);
    end
    [muh_FA,sigma_FA,CI]=normfit(FA_vett);
    mean_FA_map(voxelInd)=muh_FA;
    mean_ADC_map(voxelInd)=mean(ADC_vett);
    CI_map(voxelInd)=CI(2)-CI(1);
    CV_map(voxelInd)=sigma_FA/muh_FA;
    %waitbar(contI/nVoxel,waitbar_obj)
end
%close(waitbar_obj)
bootstrap_results.FA_map = mean_FA_map;
bootstrap_results.ADC_map = mean_ADC_map;
bootstrap_results.CI_map = CI_map;
bootstrap_results.CV_map = CV_map;
result_check = true;








% -------------------------------------------------------------------------
% |                             MASK_DTI_DATA                             |
% -------------------------------------------------------------------------
function [mask,DTI_data,result_check] = mask_DTI_data(DTI_data_unmasked,options)
% Compute the mask and mask the DTI data

% Look for mask file in the options
if isfield(options,'MaskFile')
    % There is a mask file provided
    if exist(options.MaskFile,'file')
        display('Mask provided in a file')
        file_obj = load_untouch_nii(options.MaskFile);
        mask = double(file_obj.img);
        already_mask=true;
    else
        display('WARNING: mask file expeted, but file not found')
        display(options.MaskFile)
        already_mask = false;
    end
else
    already_mask= false;
end
% If a file is not provided the mask is computed from the data
if not(already_mask)
    display('Mask computed from data...')
    [nR,nC,nS]=size(DTI_data_unmasked(1,1).VoxelData);
    nVol = length(DTI_data_unmasked);
    % Computing the summed volume
    non_weighted_sum = zeros(nR,nC,nS);
    nNonWeighted=0;
    for cont=1:nVol
        if (DTI_data_unmasked(cont).Bvalue == 0)
            non_weighted_sum = non_weighted_sum + double(DTI_data_unmasked(cont).VoxelData);
            nNonWeighted=nNonWeighted+1;
        end
    end
    non_weighted_sum = non_weighted_sum/nNonWeighted;
    % Computing the mask
    % TO FINISH
    mask_threshold = 600;
    mask = non_weighted_sum>=mask_threshold;
end
% Masking the data
DTI_data=DTI_data_unmasked;
nVol = length(DTI_data_unmasked);
out_mask_voxels = find(mask==0);
for cont=1:nVol
    DTI_data(cont).VoxelData(out_mask_voxels)=0;
end
result_check=true;








% -------------------------------------------------------------------------
% |                          GENERATE_SUBSAMPLES                          |
% -------------------------------------------------------------------------
function [subsample_list,result_check]=generate_subsamples(DTI_data,options)
% Generate the subsamples for the bootstrap analysis
% - all non weighted volumes are preserved
% - no check is performed for the distributions of the bvecs

% Looking for the non weighted volumes
nVol = length(DTI_data);
bvals_vect = zeros(1,nVol);
for cont = 1:nVol
    bvals_vect(cont)=DTI_data(cont).Bvalue;
end
index_weighted = find(bvals_vect);
index_non_weighted = find(bvals_vect==0);
% Checking the data
nWeighted = length(index_weighted);
nNonWeighted = length(index_non_weighted);
nElements = options.subsample_elements;
nSubsamples = options.number_of_subsamples;
if (nElements<=nWeighted) && (nchoosek(nWeighted,nElements)>=nSubsamples)
    % Permutations can be computed
    exitFlag = false;
    permutation_table = zeros(nSubsamples,nElements+nNonWeighted);
    pos = 1;
    while ~exitFlag
        % Generating a candidate
        permutation = randperm(nWeighted);
        candidate = [index_non_weighted, index_weighted(permutation(1:nElements))];
        % Candidate check
        candidate_check = true;
        for cont = 1:pos-1
            if length(unique([permutation_table(cont,:),candidate]))<=nElements
                candidate_check = false;
                breack
            end
        end
        % Adding the candidate if accepted
        if candidate_check
            % Candidate accepted
            permutation_table(pos,:)=candidate;
            pos = pos +1;
            if pos>nSubsamples
                % All required permutations have been computed
                exitFlag=true;
            end
        end
    end
    subsample_list=permutation_table;
    result_check=true;
else
    % There are not enough elements to compute the required subsamples
    display('ERROR: cannot compute the permutations for bootstrap analysis')
    display(['non weighted volumes:  ',num2str(nWeighted)])
    display(['elements in each subsample:  ',num2str(nElements)])
    display(['number of subsamples: ',num2str(nSubsamples)])
    subsample_list=0;
    result_check=false;
end








% -------------------------------------------------------------------------
% |                           DTI_DATA_LOADING                            |
% -------------------------------------------------------------------------
function [DTI_data,resultFlag]= DTI_data_loading(dti_data_file,bvals_file,bvecs_file,options)
% Read the data from the input files and store them as required by the DTI
% quantification tool.

if (exist(dti_data_file,'file') && exist(bvals_file,'file') && exist(bvecs_file,'file'))
    % read DTI Data
    file_obj = load_untouch_nii(dti_data_file);
    volumes_data = file_obj.img;
    [nR,nC,nS,nVol] = size(volumes_data);
    clear file_obj;
    % read Bvals
    bvals_values = textread(bvals_file,'', 'delimiter',' ','emptyvalue', NaN);
    [nR1,nBvals] = size(bvals_values);
    % read Bvecs
    bvecs_values = textread(bvecs_file,'', 'delimiter',' ','emptyvalue', NaN);
    [nR2,nBvecs] = size(bvecs_values);
    if (nR1==1) && (nR2==3) && (nBvecs==nBvals) && (nBvecs==nVol)
        % Data are consistent
        for cont = 1:nVol
            DTI_data(1,cont)=struct('VoxelData',volumes_data(:,:,:,cont),'Gradient',bvecs_values(:,cont)','Bvalue',bvals_values(1,cont));
        end
        resultFlag=true;
    else
        % Data are not consistent
        display('ERROR: data are not consistent')
        display(['DTI data:  ',num2str(nR),'x',num2str(nC),'x',num2str(nS),'x',num2str(nVol),'x'])
        display(['b values:  ',num2str(nR1),'x',num2str(nBvals)])
        display(['b vectors: ',num2str(nR2),'x',num2str(nBvecs)])
        resultFlag = false;
        DTI_data = 0;
    end
else
    % Files not found
    display('ERROR: file not found')
    display(dti_data_file)
    display(bvals_file)
    display(bvecs_file)
    resultFlag = false;
    DTI_data = 0;
end