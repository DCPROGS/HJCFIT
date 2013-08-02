Likelihood :math:`L(Q)`: 
========================

The callable objects :math:`L_{\log10}(Q)` are an interface designed for obtaining the
log10-likelihood as a function of the rates, from the knowledge of the observed open and shut time
intervals and the (type of) initial and final occupancies. This object is especially suited for
maximizing the likelihood of a given mechanism. It is defined as 

.. math:: 

  L_{\log10}(Q) = \frac{1}{\mathrm{ln} 10}\sum_b \mathrm{ln} L_b(Q)

with :math:`L_b(Q)` defined :ref:`here <log10likelihood_equation>`, the index :math:`b` refers
to the individual bursts, and :math:`\mathrm{ln}` is the natural logarithm.

The purpose of this class is to provide an interface for maximizing the likelihood. It computes,
for a given set of observed open and shut intervals, the likelihood :math:`\frac{\ln(L(Q))}{ln 10}`,
where :math:`L(Q)` is declared in :ref:`the likelihood equation <log10likelihood_equation>`. 

A callable object :math:`L(Q)` exists in both :ref:`c++ <cpp_log10_section>` and :ref:`python
<python_log10_section>`. It can be initialized as follows

:python: 

  .. literalinclude:: code/init_log10.py
     :language: python
     :lines: 2-15


:c++11:

  .. literalinclude:: code/init_log10.cc
     :language: c++
     :lines: 1-14, 27-

  The initialisation of `bursts` above is done in using two newer c++11 coding techniques: 
  `initializer lists <initializerlist_>`_, and `uniform initialization <uniforminit_>`_.
  It may not be available from all compilers just yet...

Once the objects are initialized, the input attributes can be accessed (and modified) through the
'.' operator: `likelihood.nopen = 2`. 

It is required that the bursts have been pre-filtered so that there are no intervals smaller than
the resolution :math:`\tau`. This can be done using :c:func:`time_filter` (:py:func:`time_filter` or
:cpp:func:`interval_filter` (:cpp:func:`interval_filter`).


The likelihood for any Q-matrix can then be computed by calling the `likelihood` objects as though
they were function. The following snippets are inserted at the tail end of the previous code.

:python:

  .. literalinclude:: code/init_log10.py
     :language: python
     :lines: 15, 15-25 
  
  The function can take any sort square matrix, whether using standard python lists or a numpy
  array. It can only take one matrix at a time. 

:c++11:

  .. literalinclude:: code/init_log10.cc
     :language: c++
     :lines: 14-25

  
The return is the log-likelihood associated with the bursts and the input Q-matrix. In both python
and c++, the functions accepts either a matrix or an actual :cpp:class:`QMatrix` (:py:class:`QMatrix`)
object. In the former case, the number of open states is set to `nopen`.

Finally, it should be noted that the python the bursts are accessed in python directly from the
likelihood using normal sequence operations. Only a small subset of sequence operations where
implemented.

:python:

  .. literalinclude:: code/init_log10.py
     :language: python
     :lines: 1, 26-

.. _initializerlist: https://en.wikipedia.org/wiki/C++11#Initializer_lists
.. _uniforminit: https://en.wikipedia.org/wiki/C++11#Uniform_initialization
