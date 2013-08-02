.. _cpp_log10_section:

Log10Likelihood
---------------

The purpose of this class is to provide an interface for makinimizing the likelihood. It computes,
for a given set of observed open and shut intervals, the likelihood :math:`\frac{\ln(L(Q))}{ln 10}`,
where :math:`L(Q)` is declared in :ref:`the likelihood equation <log10likelihood_equation>`. 

.. doxygenclass:: DCProgs::Log10Likelihood
  :members:
