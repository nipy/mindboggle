README for the Mindboggle repository

In the directory tree below, ? stands for multiple entries,
such as ``?h'' for left and right hemispheres: [lh, rh]

mindboggle/
    README.txt   # this file
    pipeline.py  # main program (nipype pipeline)
    atlases.py   # functions for multi-atlas surface labeling
    features.py  # functions for feature-based surface labeling
    label_volume.py  # functions for label filling a volume
    data/
        templates/
            templates_freesurfer/  # FreeSurfer-ready templates
        atlases/
            labels.txt  # list of atlas labels
            info/  # atlas scanning and participant information
            subjects/
                <subject>/
                    originalT1.nii.gz
                    ?h.pial.[depth, curv…].vtk
                    ?.inflated.vtk
                    ?h.sphere_to_?_template.reg
                    aparcNMMjt+aseg.vtk
                    aparcNMMjt+aseg.nii.gz
            copy_contents_to_freesurfer_subjects/ # (see below)
                <subject>/
                    label/
                        ?h.aparcNMMjt.annot  # manual label files
                    surf/
                        ?h.[pial, inflated]  # surfaces
                        ?h.sphere_to_?_template.reg  # template transforms
