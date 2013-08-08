.. _manual_occupancies:

Occupancies
===========

The start and end occupancies can be computed one of two way. They can be the equilibrium
occupancies determined from the equation: 

.. math::

  \phi_A = \phi_A {}^e\mathcal{G}_{AF} {}^e\mathcal{G}_{FA},\\
  \phi_F = {}^e\mathcal{G}_{FA} {}^e\mathcal{G}_{AF} \phi_F,\\

subject to the constraints,

.. math::

 \sum_i [\phi_A]_i = 1,\\
 \sum_i [\phi_F]_i = 1.

Where :math:`{}^e\mathcal{G}_{AF}` and :math:`{}^e\mathcal{G}_{FA}` are the laplacians of the
missed-events likelihoods (or equivilantly, the ideal likelihoods) for :math:`s=0`, and
:math:`[a]_i` indicates the :math:`i^{th}` component of vector :math:`a`.

Or they can be computed as CHS vectors, e.g. equation 5.11 and 5.8 from :cite:`colquhoun:1996`. 

The occupancies are accessed differently in :ref:`c++ <cpp_occupancies_api>` and in python.

:python:

  The equilibrium occupancies are accessed as properties of :py:class:`~dcprogs.likelihood.IdealG`
  and :py:class:`~dcprogs.likelihood.MissedEventsG` instances. The CHS vectors are functions of
  these same classes that take as arguments the critical time.

   
  .. literalinclude:: ../../code/occupancies.py
     :language: python


:c++11:

  Both equilibrium and CHS occupancies are accessed via function calls acting on the
  :cpp:class:`IdealG` and :cpp:class:`MissedEventsG`.

  .. literalinclude:: ../../code/occupancies.cc
     :language: c++

In c++, the occupancies are kept outside of the classes because computing these values is outside
the pure remit of the classes (which is to compute the likelihood). However, in python, practicality
beats purity, and it makes practical sense to keep likelihood and occupancies together.

.. note::

   :py:func:`~dcprogs.likelihood.Log10Likelihood` uses equilibrium occupancies depending on the
   value of its attribute :py:attr:`~dcprogs.likelihood.Log10Likelihood.tcritical`:

   - if it is ``None``, ``numpy.NaN``, or negative, then the equilibrium occupancies are used
   - if it a strictly positive real number, then the CHS vectors are computed
