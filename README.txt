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
                    labels.manual.vtk
                    labels.manual.nii.gz
                    labels.ants.subcortex.nii.gz
                    labels.freesurfer.subcortex.nii.gz
            subjects_freesurfer/ # (see below)
                <subject>/
                    label/
                        ?h.labels.manual.annot  # manual label files
                    surf/
                        ?h.[pial, inflated]  # surfaces
                        ?h.sphere_to_?_template.reg  # template transforms
