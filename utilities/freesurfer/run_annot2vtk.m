
subjdir = '/Users/arno/Documents/Projects/mindboggle/data/ManualSurfandVolLabels/subjects/';

subjects = {'KKI2009-16';
'KKI2009-33';
'OAS1_0101_MR1';
'OAS1_0395_MR1';
'KKI2009-01';
'KKI2009-17';
'KKI2009-34';
'OAS1_0111_MR1';
'KKI2009-01_fs5.0';
'KKI2009-18KKI2009-35';
'OAS1_0117_MR1';
'plos_CJ_700_3_1';
'KKI2009-02';
'KKI2009-19';
'KKI2009-36';
'OAS1_0145_MR1';
'plos_Ding_1009_3_1';
'KKI2009-03';
'KKI2009-20';
'KKI2009-37';
'OAS1_0150_MR1';
'plos_Ding_1035_3_1';
'KKI2009-04';
'KKI2009-21';
'KKI2009-38';
'OAS1_0156_MR1';
'plos_Ding_1042_3_1';
'KKI2009-05';
'KKI2009-22';
'KKI2009-39';
'OAS1_0191_MR1';
'plos_Ding_1050_3_1';
'KKI2009-06';
'KKI2009-23';
'KKI2009-40';
'OAS1_0202_MR1';
'plos_Ding_1052_3_1';
'KKI2009-07';
'KKI2009-24';
'KKI2009-41';
'OAS1_0230_MR1';
'plos_Ding_756_3_1';
'KKI2009-08';
'KKI2009-25';
'KKI2009-42';
'OAS1_0236_MR1';
'plos_Language_587_3_1';
'KKI2009-09';
'KKI2009-26';
'MMR2011A_3T';
'OAS1_0239_MR1';
'plos_Language_795_3_1';
'KKI2009-10';
'KKI2009-27';
'MMR2011A_7T';
'OAS1_0249_MR1';
'plos_Language_820_3_1';
'KKI2009-11';
'KKI2009-28';
'MMR2011B_3T';
'OAS1_0285_MR1';
'plos_Language_886_3_1';
'KKI2009-12';
'KKI2009-29';
'MMR2011B_7T';
'OAS1_0353_MR1';
'plos_Language_888_3_1';
'KKI2009-13';
'KKI2009-30';
'OAS1_0061_MR1';
'OAS1_0368_MR1';
'KKI2009-14';
'KKI2009-31';
'OAS1_0080_MR1';
'KKI2009-15';
'KKI2009-32';
'OAS1_0092_MR1';
'OAS1_0379_MR1'};

hemis = {'lh', 'rh'};

for isubjects = 1:length(subjects)
    subject = subjects{isubjects}; 
    for ihemis = 1:length(hemis)
      hemi = hemis{ihemis};
      annot = [subjdir,subject,'/label/',hemi,'.aparcNMMjt.annot'];
      if(exist(annot)>0)
        annot
        
        pial = [subjdir,subject,'/surf/',hemi,'.pial'];
        white = [subjdir,subject,'/surf/',hemi,'.white'];
        inflated = [subjdir,subject,'/surf/',hemi,'.inflated'];
        sphere = [subjdir,subject,'/surf/',hemi,'.sphere'];

        if(exist(pial)>0)
          annot2vtk(annot, pial);
        end
        if(exist(white)>0)
          annot2vtk(annot, white);
        end
        if(exist(inflated)>0)
          annot2vtk(annot, inflated);
        end
        if(exist(sphere)>0)
          annot2vtk(annot, sphere);
        end
      end
    end
end
