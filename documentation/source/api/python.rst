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


Internal Type
-------------

Many of the python bindings accept both scalars and values that can be converted to numpy arrays.
The latter are first converted to a specific numpy type, depending on how this package was compiled.

.. py:data:: dcprogs.internal_dtype

   Numpy type used internally

If the package is configured and compiled with ``-DCPROGS_LONG_DOUBLE`` (see wiki), then long
doubles of 128bits are used through out. Otherwise,  the package defaults to garden variety 64bit
doubles.
