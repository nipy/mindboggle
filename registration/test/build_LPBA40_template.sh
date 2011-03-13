# Run inside the LPBA40 brain-extracted folder:
# -c 0 (non-parallel, 1 for sge)

bash buildtemplateparallel.sh -d 3 -m 30x100x10 -t GR -s CC -c 0 -j 4 -o LPBA40 -z InitialTemplate.nii.gz  *.nii.gz