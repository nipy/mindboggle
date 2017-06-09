.. _singularity:

------------------------------------------------------------------------------
 How do I create and use a Singularity image rather than Docker?
------------------------------------------------------------------------------

To convert the Mindboggle Docker image to a Singularity image, run::

  singularity pull mindboggle.img docker://nipy/mindboggle

After the conversion, you'll need 3 mounts (data, working directory, and
home directory) for it to work. Here is an example of the full command
together with some of the optional arguments (thanks, Satra!), using the same
environment variables as in the `README <http://mindboggle.info/software.html>`_::

    singularity run \
      -B $HOST:$DOCK:ro \
      -B $PWD:$DOCK \
      -B $PWD/jovyan:/home/jovyan \ 
      -e mindboggle.img \
      $DOCK/example_mri_data/T1.nii.gz \
      --id arno \
      --fs_T2image $DOCK/example_mri_data/T2.nii.gz \
      --plugin MultiProc --plugin_args "dict(n_procs=2)" \
      --fs_openmp 5 --ants_num_threads 5 --mb_num_threads 10


