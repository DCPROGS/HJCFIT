Python API Reference
====================

.. toctree::
   :maxdepth: 2

   python/qmatrix.rst
   python/log10likelihood.rst
   python/missed_eventsG.rst
   python/idealG.rst
   python/exact_survivor.rst
   python/approx_survivor.rst
   python/determinanteq.rst
   python/roots.rst


Extras
------

Internal Numpy Type
+++++++++++++++++++

Many of the python bindings accept both scalars and values that can be converted to numpy arrays.
The latter are first converted to a specific numpy type, depending on how this package was compiled.

.. py:data:: dcprogs.internal_dtype

   Numpy type used internally

If the package is configured and compiled with ``-DCPROGS_LONG_DOUBLE`` (see wiki), then long
doubles are used through out. The exact size of long doubles is platform dependent
and is typically 80-bit on Linux and OSX and 64-bit on Windows. See `Wikipedia
<https://en.wikipedia.org/wiki/Long_double>`__ for more information.
Otherwise, the package defaults to garden variety 64-bit doubles.


Linear Algebra
++++++++++++++

Numpy only provides a subset of its utilities for arrays consisting of long doubles. As a
result, this package exposes some of Eigen_'s capabilities, as needed. Their interface is generally
reminiscent of the numpy utility they mirror.


.. currentmodule:: dcprogs.likelihood
.. autofunction:: eig
.. autofunction:: inv
.. autofunction:: svd
.. autofunction:: det
.. autofunction:: expm
