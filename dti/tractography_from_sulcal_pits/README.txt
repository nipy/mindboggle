TRACTOGRAPHY FROM SULCAL PITS

The method (tract_sulcal_pits.py) perform the tractography using the sulcal pits as seed and target.
It moves the sulcal pits from T1 to DTI space and dilate them to make a seed.
Then it perform tractography using fsl.