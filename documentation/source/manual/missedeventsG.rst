.. _manual_eG:

Missed-Events Likelihood :math:`{}^eG(t)`
=========================================

The callable object :cpp:class:`DCProgs::MissedEventsG` provides an interface to compute the
likelihood :math:`{}^eG(t)` of open and shut events as a function of their lengths, for a fixed
:math:`Q`-matrix. It has the ability to compute both exact and approximate missed-events likelihood,
returning one or the other depending on a given time cutoff.

The asymptotic expression of the likelihood can be computed from the knowledge of the roots of a
specific equations. On the one hand, root-finding can be a fairly difficult numerical operation. On
the other, it would be more convenient if we can initialize :cpp:class:`DCProgs::MissedEventsG`
directly from a :math:`Q`-matrix object. As such, there are several means to initialize the functor:

- from the knowledge of the roots and the determinant equations
- directly from a :math:`Q`-matrix, using the default root-finding mechanism
- from a :math:`Q`-matrix, using a custom root-finding mechanism (c++ only)


Initialization from a :math:`Q`-matrix
""""""""""""""""""""""""""""""""""""""

:python: 

  .. literalinclude:: ../../code/missedeventsG.py
     :language: python
     :lines: 1-11, 19


:c++11:

  .. literalinclude:: ../../code/missedeventsG.cc
     :language: c++
     :lines: 1-18, 41, 55-

A fair amount of work goes on behind the scene. First reasonable upper and lower bounds for the
roots obtained (:cpp:func:`DCProgs::find_lower_bound_for_roots`, and
:cpp:func:`DCProgs::find_upper_bound_for_roots`). Then intervals for each roots are computed
(:cpp:func:`DCProgs::find_root_intervals`). And finally, the roots themselves are obtained
(:cpp:func:`DCProgs::brentq`). All this work is done automatically in the case of this particular
instantiation. A few extra parameters to control the root-finding process can be passed to the c++
and python constructors.


Initialization from the roots and determinant equations
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

:python: 

  .. literalinclude:: ../../code/missedeventsG.py
     :language: python
     :lines: 13-16


:c++11:

  .. literalinclude:: ../../code/missedeventsG.cc
     :language: c++
     :lines: 20-31

In this case, it is expected the roots are known, somehow, as well as their multiplicity. This
format allows external root-finding methods to be interfaced with the packages.


Initialization from the :math:`Q`-matrix and a root finding function
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Given a root-finding function, it is possible to instantiate  :math:`{}^eG`. The root finding
function should take a determinant equation as input, and return a vector of
:cpp:class:`DCProgs::Root` as output. In the code below, we show how the prior initialization could
be recreated.

.. literalinclude:: ../../code/missedeventsG.cc                
   :language: c++                                              
   :lines: 35-38

This is mostly a convenience function, to make it slightly easier to interface with other
root-finding methods in c++. This interface is not explicitly available in python, although it can
be created with ease.



Computing :math:`{}^eG_{AF}(t)` and :math:`{}^eG_{FA}(t)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The computation of the likelihood matrices is illustrated by comparing the three initialization
methods (two in python). It is clear that all three yield the same function.

:python: 

  .. literalinclude:: ../../code/missedeventsG.py
     :language: python
     :lines: 22-


:c++11:

  .. literalinclude:: ../../code/missedeventsG.cc
     :language: c++
     :lines: 45-54


The python bindings accept any input that can be transformed to a numpy array of reals. If the input
is a scalar, then the AF and FA blocks are returned. If the input is an array, then an array of
similar shape is returned, where each component is a matrix. 

The :cpp:class:`DCProgs::MissedEventsG` provides further functionality. For instance, the cutoff
point between exact and asymptotic calculations can be set explicitly (it defaults to :math:`t <
3\tau`). And the likelihood can be computed in Laplace space (see
:cpp:member:`DCProgs::MissedEventsG::laplace_af` and
:cpp:member:`DCProgs::MissedEventsG::laplace_fa`). We invite users to turn to the :ref:`python
<python_eG_api>` and the :ref:`c++ <cpp_eG_api>` API for more details.
