Ideal Likelihood :math:`\mathcal{G}(t)`
=======================================

A wrapper around :py:class:`~dcprogs.likelihood.QMatrix` is provided which allows the calculation of
the ideal likelihood:

.. math::

  \mathcal{G}(t) = \left(
  \begin{eqnarray}
     e^{-t\mathcal{Q}_{af}} & 0 \\
     0& e^{-t\mathcal{Q}_{fa}} \\
  \end{eqnarray}
  \right)


This object can be initialized directly from a :py:class:`QMatrix`.

:python:

  .. literalinclude:: ../../code/idealG.py
     :language: python
     :lines: 2-11


:c++11:

  .. literalinclude:: ../../code/idealG.cc
     :language: c++
     :lines: 1-20, 30-

It provides the ideal likelihood as a function of time, as well as the laplace transforms:

:python:

  .. literalinclude:: ../../code/idealG.py
     :language: python
     :lines: 1, 12-

  .. note::

     This package can be compiled to use real number of 128bits. However, despite providing arrays
     for reals of that size, numpy does not include the corresponding linear algebra functions. A
     few useful functions, such as ``expm`` in this example, are provided to remediate to this
     situation.

:c++11:

  .. literalinclude:: ../../code/idealG.cc
     :language: c++
     :lines: 22-28
