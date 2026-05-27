.. _manual_vectors:

Equilibrium vectors
===========

The start and end equilibrium vectors can be computed in two ways depending on the burst/cluster nature. They can be the equilibrium
vectors determined from the equation: 

.. math::

  \phi_A = \phi_A {}^e\mathcal{G}_{AF} {}^e\mathcal{G}_{FA},\\
  \phi_F = {}^e\mathcal{G}_{FA} {}^e\mathcal{G}_{AF} \phi_F,\\

subject to the constraints,

.. math::

 \sum_i [\phi_A]_i = 1,\\
 \sum_i [\phi_F]_i = 1.

Where :math:`{}^e\mathcal{G}_{AF}` and :math:`{}^e\mathcal{G}_{FA}` are the laplacians of the
missed-events transition densities :math:`{}^e\mathcal{G}(t)` (or equivalently, the ideal transition densities) for :math:`s=0`, and
:math:`[a]_i` indicates the :math:`i^{th}` component of vector :math:`a`.

On the other hand equilibrium vectors can be computed as CHS vectors, e.g. equation 5.11 and 5.8 from :cite:`colquhoun:1996`. 

The vectors are accessed differently in :ref:`C++ <cpp_vectors_api>` and in Python.

:Python:

  The equilibrium vectors are accessed as properties of :py:class:`~HJCFIT.likelihood.IdealG`
  and :py:class:`~HJCFIT.likelihood.MissedEventsG` instances. The CHS vectors are functions of
  these same classes that take as arguments the critical time.

   
  .. literalinclude:: ../../code/vectors.py
     :language: python


:C++11:

  All vectors are accessed via function calls acting on the
  :cpp:class:`IdealG` and :cpp:class:`MissedEventsG`.

  .. literalinclude:: ../../code/vectors.cc
     :language: c++

In C++, the vectors are kept outside of the classes because computing these values is outside
the pure remit of the classes (which is to compute the likelihood). However, in Python, practicality
beats purity, and it makes practical sense to keep likelihood and equilibrium vectors together.

.. note::

   :py:func:`~HJCFIT.likelihood.Log10Likelihood` uses equilibrium vectors depending on the
   value of its attribute :py:attr:`~HJCFIT.likelihood.Log10Likelihood.tcritical`:

   - if it is ``None``, ``numpy.NaN``, or negative, then the equilibrium vectors are used;
   - if it is a strictly positive real number, then the CHS vectors are computed.
