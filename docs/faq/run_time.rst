.. _run_time:

------------------------------------------------------------------------------
 How long does it take to run Mindboggle?
------------------------------------------------------------------------------

We have tested the software most extensively with Python 3.5.1 on Ubuntu Linux 14.04.
Running Mindboggle on the Docker installation on a Macbook Pro (2.6 GHz Intel Core i7
with 16 GB memory; macOS 10.12) took about 100 minutes, of which 20 minutes were spent
optionally computing Laplace-Beltrami spectra and Zernike moments.

When only the surface shapes of gyrus labels were computed, 
without Laplace-Beltrami spectra or Zernike moments, 
Mindboggle took less than 7 minutes using the following command::

    mindboggle $FREESURFER_SUBJECT --out $MINDBOGGLED \
     --no_volumes --no_sulci --no_moments --no_spectra
