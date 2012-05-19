Node: Convert_surfaces (utility)
================================

 Hierarchy : pipeline.Convert_surfaces
 Exec ID : Convert_surfaces.a0

Original Inputs
---------------

* fs_surface_files : ['/projects/mindboggle/data/KKI2009-11/surf/lh.pial', '/projects/mindboggle/data/KKI2009-11/surf/rh.pial']
* function_str : S'def convert_to_vtk(fs_surface_files):\n    """Measure\n\n    measure_()\n    """\n    #import subprocess as sp\n    import os\n    if type(fs_surface_files) is str:\n        dummy = []\n        dummy.append(fs_surface_files)\n        fs_surface_files = dummy\n    surface_files = []\n    for fs_surface_file in fs_surface_files:\n        surface_file = fs_surface_file + \'.vtk\'\n        surface_files.append(surface_file)\n        cmd = [\'mris_convert\', fs_surface_file, surface_file]\n        os.system(\' \'.join(cmd))\n    #proc = sp.Popen(\' \'.join(cmd))\n    #o, e = proc.communicate()\n    #if proc.returncode > 0 :\n    #    raise Exception(\'\\n\'.join([cmd + \' failed\', o, e]))\n    return surface_files\n'
.
* ignore_exception : False

Execution Inputs
----------------

* fs_surface_files : ['/projects/mindboggle/data/KKI2009-11/surf/lh.pial', '/projects/mindboggle/data/KKI2009-11/surf/rh.pial']
* function_str : S'def convert_to_vtk(fs_surface_files):\n    """Measure\n\n    measure_()\n    """\n    #import subprocess as sp\n    import os\n    if type(fs_surface_files) is str:\n        dummy = []\n        dummy.append(fs_surface_files)\n        fs_surface_files = dummy\n    surface_files = []\n    for fs_surface_file in fs_surface_files:\n        surface_file = fs_surface_file + \'.vtk\'\n        surface_files.append(surface_file)\n        cmd = [\'mris_convert\', fs_surface_file, surface_file]\n        os.system(\' \'.join(cmd))\n    #proc = sp.Popen(\' \'.join(cmd))\n    #o, e = proc.communicate()\n    #if proc.returncode > 0 :\n    #    raise Exception(\'\\n\'.join([cmd + \' failed\', o, e]))\n    return surface_files\n'
.
* ignore_exception : False

Execution Outputs
-----------------

* surface_files : ['/projects/mindboggle/data/KKI2009-11/surf/lh.pial.vtk', '/projects/mindboggle/data/KKI2009-11/surf/rh.pial.vtk']

Runtime info
------------

* duration : 5.27620220184
* hostname : bijli.local

Environment
~~~~~~~~~~~

* AFNIHOME : /Users/arno/Software/AFNI
* ANTSPATH : /Users/arno/Software/ANTS_1.9/bin/
* Apple_PubSub_Socket_Render : /tmp/launch-0p6yZ1/Render
* BOOST_ROOT : /Users/arno/Software/boost_1_43_0
* C3DDIR : /Users/arno/Software/c3D
* CAMINO : /Users/arno/Software/camino
* CAMINO_HEAP_SIZE : 3000
* COMMAND_MODE : unix2003
* DISPLAY : /tmp/launch-jzcUz6/org.x:0
* EDITOR : /usr/bin/emacs
* FIX_VERTEX_AREA : 
* FMRI_ANALYSIS_DIR : /Applications/freesurfer/fsfast
* FREESURFER_HOME : /Applications/freesurfer
* FSFAST_HOME : /Applications/freesurfer/fsfast
* FSF_OUTPUT_FORMAT : nii.gz
* FSLCONFDIR : /usr/local/fsl/config
* FSLDIR : /usr/local/fsl
* FSLLOCKDIR : 
* FSLMACHINELIST : 
* FSLMACHTYPE : apple-darwin10-gcc4.2
* FSLMULTIFILEQUIT : TRUE
* FSLOUTPUTTYPE : NIFTI_GZ
* FSLREMOTECALL : 
* FSLTCLSH : /usr/local/fsl/bin/fsltclsh
* FSLWISH : /usr/local/fsl/bin/fslwish
* FSL_BIN : /usr/local/fsl/bin
* FS_FREESURFERENV_NO_OUTPUT : 1
* FS_OVERRIDE : 0
* FUNCTIONALS_DIR : /Applications/freesurfer/sessions
* GITDIR : /usr/local/git
* GRIN_ARGS : -C 2 --no-skip-dirs
* HOME : /Users/arno
* IMAGEDIR : /Software/ImageMagick-6.6.9/bin
* LANG : en_US.UTF-8
* LOCAL_DIR : /Applications/freesurfer/local
* LOGNAME : arno
* LSCOLORS : ExgxfxfxCxDxDxCxCxaCaC
* MACOSX_DEPLOYMENT_TARGET : 10.5
* MANPATH : /Users/arno/Software/camino/man:/Users/arno/Software/camino/man:
* MINC_BIN_DIR : /Applications/freesurfer/mni/bin
* MINC_LIB_DIR : /Applications/freesurfer/mni/lib
* MKL_NUM_THREADS : 1
* MNI_DATAPATH : /Applications/freesurfer/mni/data
* MNI_DIR : /Applications/freesurfer/mni
* MNI_PERL5LIB : /Applications/freesurfer/mni/lib/../System/Library/Perl/5.8.6
* OLDPWD : /projects/mindboggle/mindboggle/measure
* OS : Darwin
* PATH : /Library/Frameworks/Python.framework/Versions/Current/bin:/Library/Frameworks/EPD64.framework/Versions/Current/bin:/Library/Frameworks/EPD64.framework/Versions/Current/bin:/Users/arno/Software/boost_1_43_0:/Users/arno/Software/VTK-build/bin:/Users/arno/Software/VTK-build:/Applications/freesurfer/bin:/Applications/freesurfer/fsfast/bin:/Applications/freesurfer/tktools:/usr/local/fsl/bin:/Applications/freesurfer/bin/freeview.app/Contents/MacOS/:/Applications/freesurfer/mni/bin:/Applications/freesurfer:/usr/local/fsl/bin:/Software/ImageMagick-6.6.9/bin:/Users/arno/Software/camino:/Users/arno/Software/camino/bin:/Users/arno/Software/c3D/bin:/Users/arno/Software/AFNI:/usr/local/git/bin:/Library/Frameworks/Python.framework/Versions/Current/bin:/Library/Frameworks/EPD64.framework/Versions/Current/bin:/Library/Frameworks/EPD64.framework/Versions/Current/bin:/Users/arno/Software/boost_1_43_0:/Users/arno/Software/VTK-build/bin:/Users/arno/Software/VTK-build:/Applications/freesurfer/bin:/Applications/freesurfer/fsfast/bin:/Applications/freesurfer/tktools:/usr/local/fsl/bin:/Applications/freesurfer/bin/freeview.app/Contents/MacOS/:/Applications/freesurfer/mni/bin:/Applications/freesurfer:/usr/local/fsl/bin:/Software/ImageMagick-6.6.9/bin:/Users/arno/Software/camino:/Users/arno/Software/camino/bin:/Users/arno/Software/c3D/bin:/Users/arno/Software/AFNI:/usr/local/git/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/local/git/bin:/usr/texbin:/usr/X11/bin
* PERL5LIB : /Applications/freesurfer/mni/lib/../System/Library/Perl/5.8.6
* PWD : /projects/mindboggle/mindboggle
* PYTHONPATH : :/Library/Frameworks/Python.framework/Versions/7.1/bin:/Library/Frameworks/Python.framework/Versions/7.1/bin
* SECURITYSESSIONID : 31f171
* SHELL : /bin/bash
* SHLVL : 1
* SSH_AUTH_SOCK : /tmp/launch-tXUrWc/Listeners
* SUBJECTS_DIR : /Applications/freesurfer/subjects
* TERM : xterm-color
* TERM_PROGRAM : Apple_Terminal
* TERM_PROGRAM_VERSION : 273.1
* TMPDIR : /var/folders/mR/mR3bjLb-HhCu36omEdpWSU+++TQ/-Tmp-/
* USER : arno
* VTK : /Users/arno/Software/VTK-build
* _ : /Library/Frameworks/Python.framework/Versions/Current/bin/python
* __CF_USER_TEXT_ENCODING : 0x1F7:0:0

