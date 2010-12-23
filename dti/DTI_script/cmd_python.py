# Common commands in the script

# How to import subject data from min_burst, mri_bank, etc...
python 
/mind_cart/DTI/DTI_script/import_data.py # Script
/mind_cart/DTI/50548/ 
/data/BI/mindburst/data_pieces/mri_atropos_crop/2/71/40/brain_n3.nii.gz # Cropped Brain
/data/BI/mindburst/data_pieces/mri_mask/2/71/42/mask.nii.gz # mask
/data/BI/human/MRIBANK/50548/raw/s006a001_TENSOR_ASSET_25_Dir.nii.gz #DTI data

# e.g.
python /mind_cart/DTI/DTI_script/import_data.py /mind_cart/DTI/50548/ /data/BI/mindburst/data_pieces/mri_atropos_crop/2/71/40/brain_n3.nii.gz /data/BI/mindburst/data_pieces/mri_mask/2/71/42/mask.nii.gz /data/BI/human/MRIBANK/50548/raw/s006a001_TENSOR_ASSET_25_Dir.nii.gz




# -------------------------------------------------------------------------------------------------------------
# Pre elaboration Steps
# To perform the whole pre-elaboration
python /mind_cart/DTI/DTI_script/full_preparation_FSL.py /mind_cart/DTI/50192

# Single steps:
# Coregistration and mask step
python /mind_cart/DTI/DTI_script/prepare_data.py /mind_cart/DTI/50192

# Formatting and Eddy current correction step
python /mind_cart/DTI/DTI_script/format_data_FSL.py /mind_cart/DTI/50192

# FA map step
python /mind_cart/DTI/DTI_script/fa_map_FSL.py /mind_cart/DTI/50192

# BedpostX step
python /mind_cart/DTI/DTI_script/preparation4tract_FSL.py /mind_cart/DTI/50192




# -------------------------------------------------------------------------------------------------------------
# Tractography Step
python 
/mind_cart/DTI/DTI_script/fsl_tractography.py # Script
/mind_cart/DTI/50440/ # Sub_dir
/mind_cart/DTI/50440/SulcalPits/pits_dilate3.nii.gz # Seeds mask
/mind_cart/DTI/50440/grant_fsl_data.sul_pits_tracts/ # Output dir

#e.g.
python /mind_cart/DTI/DTI_script/fsl_tractography.py /mind_cart/DTI/50440/ /mind_cart/DTI/50440/SulcalPits/pits_dilate3.nii.gz /mind_cart/DTI/50440/grant_fsl_data.sul_pits_tracts/




# -------------------------------------------------------------------------------------------------------------
# Graph step
python 
/mind_cart/DTI/DTI_script/make_graph.py # Script
/mind_cart/DTI/50246/fsl_data.sul_pits_tracts/ # Tract location
/mind_cart/DTI/50246/fsl_data.sul_pits_graph/ # Output folder
50246 # Sub ID

# e.g.
python /mind_cart/DTI/DTI_script/make_graph.py /mind_cart/DTI/50246/fsl_data.sul_pits_tracts/ /mind_cart/DTI/50246/fsl_data.sul_pits_graph/ 50246