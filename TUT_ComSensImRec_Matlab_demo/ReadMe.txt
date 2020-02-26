------------------------------------------------------------------------------------------

    Demo software for Compressed Sensing Image Reconstruction
              Public release ver. 1.1 (23 July 2015)

------------------------------------------------------------------------------------------

The software reproduces the experiments published in the paper

 K. Egiazarian, A. Foi, and V. Katkovnik,
 "Compressed Sensing Image Reconstruction via Recursive Spatially Adaptive Filtering,"
 Proc. IEEE ICIP 2007, San Antonio (TX), USA, pp. 549-552, Sept. 2007.
 DOI http://dx.doi.org/10.1109/ICIP.2007.4379013

------------------------------------------------------------------------------------------

authors:               Karen Egiazarian
                       Alessandro Foi

web page:              http://www.cs.tut.fi/~comsens

contact:               firstname.lastname@tut.fi

------------------------------------------------------------------------------------------
Copyright (c) 2006-2011 Tampere University of Technology.
All rights reserved.
This work should be used for nonprofit purposes only.
------------------------------------------------------------------------------------------

Contents
--------

demo_TUT_ComSensImRec_single.m  main script using single-precision mex-file
demo_TUT_ComSensImRec_double.m  main script using double-precision mex-file (Win32 only)
subscript_figure_iter.m         subscript for updating figures
subscript_text_iter.m           subscript for text printout
LineMaskLimitedAngle.m          function to generate radial masks
function_BM3D_HT_double.m       function wrapper for BM3D (double precision)
function_BM3D_HT_single.m       function wrapper for BM3D (single precision)
bm3d_thr_double.dll             BM3D hard-thresholding with double-precision floats
                                (compiled for Windows 32-bit)
bm3d_thr.*                      BM3D hard-thresholding with single-precision floats
                                (compiled for most platforms)


Requirements
------------
This demo is designed for Matlab (ver. 7.4 and above)


Description
-----------

This demo software reproduces the following four experiments from the paper

 K. Egiazarian, A. Foi, and V. Katkovnik,
 "Compressed Sensing Image Reconstruction via Recursive Spatially Adaptive Filtering,"
 Proc. IEEE ICIP 2007, San Antonio (TX), USA, pp. 549-552, Sept. 2007.
 DOI http://dx.doi.org/10.1109/ICIP.2007.4379013


Experiments 1-3
Tomographic reconstruction from sparse projections
(approximating Radon projections as radial lines in FFT domain)
 We consider three illustrative inverse problems of compressed sensing for computerized
 tomography. In particular, we show reconstruction examples of the Shepp-Logan phantom
 from sparse Radon projections, with 22 and 11 radial lines in FFT-domain
 (i.e., available Radon projections), and reconstruction from limited-angle projections,
 with a reduced subset of 61 projections within a 90 degrees aperture. The available
 portions of the spectrum and the initial back-projection estimates are shown below.

  experiment_case=1  - Sparse projections: 22 radial lines
  experiment_case=2  - Sparse projections: 11 radial lines
  experiment_case=3  - Limited-angle (90 degrees aperture, 61 radial lines)


Experiment 4
Reconstruction from low-frequency portion of Fourier spectrum
 Reconstruction of the Cameraman image (256x256 pixels) from the low-frequency portion
 of its Fourier spectrum (a 128×128 square centered at the DC).

  experiment_case=4



For all the above experiments, the block-matching and 3D filtering algorithm (BM3D)
serves as the spatially adaptive filter used in the recursive stochastic approximation
algorithm. The separable 3D Haar wavelet decomposition is adopted as the transform
utilized internally by the BM3D algorithm.

The general version of the BM3D algorithm is available at
http://www.cs.tut.fi/~foi/GCF-BM3D/

 K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian,
 "Image denoising by sparse 3D transform-domain collaborative filtering,"
 IEEE Trans. Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
 DOI http://dx.doi.org/10.1109/TIP.2007.901238


------------------------------------------------------------------------------------------

Disclaimer
----------

Any unauthorized use of these routines for industrial or profit-oriented activities is
expressively prohibited. By downloading and/or using any of these files, you implicitly
agree to all the terms of the TUT limited license (included in the file Legal_Notice.txt).
------------------------------------------------------------------------------------------
