####################
Multi-precision math
####################

The nature of the calculations involved in HJCFIT especially in the root
finding means that the code is susceptible to issues related to the limited
precision of the floating point numbers use to represent values within the 
code. To work around these issues the code is able to use GMP and MPFR to
perform specific calculations at higher precision. As multi precision math is 
not implemented in hardware it comes with very significant run time overhead, 
typically orders of magnitude slower we have added support for multi-precision 
as a fall back mechanism. Currently there is only support for fallback in 
``find_eigs_bound`` which is known to be problematic. The pattern that we 
follow is to do the calculations with regular precision floating point. If that
fails we rerun  the calculation with multi-precision. Currently only at 128 
bits but that could easily be extended via an iterative process to arbitrary
precision.


The code can be compilled with and without this feature. Currently CMake is
setup to automatically download and build ``GMP`` and ``MPFR`` on ``OSX`` and
linux if not found. However, this feature is not enabled on Windows by default.
``GMP`` is not supported on Windows but it should be possible to build and use
``MPIR`` as a drop-in replacement for ``GMP`` on windows which ``MPFR`` can be
linked against. However, there is currently no support for building this
automatically on windows. To control if this option should be enabled the flag
``DCPROGS_USE_MPFR`` can be set to ``on`` or ``off``. 
