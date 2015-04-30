# dlt
MATLAB code to compute the Discrete Legendre transform (DLT)

This repo contains MATLAB code related to the paper 
  N. Hale and A. Townsend, "A fast FFT-based discrete Legendre transform" (submitted to IMAJNA - preprint [here](http://dip.sun.ac.za/~nhale/publications/)).
In particular, it enables computation of  the discrete Legendre transform (DLT) in O( N log(N)^2 / loglog(N) ). The repo also contains all necessary codes to reproduce the results contained within the paper.

The dlt, idlt, and ndct routines are also contained in Chebfun (www.chebfun.org & https://github.com/chebfun/chebfun/) where they will be more carefully maintained.
