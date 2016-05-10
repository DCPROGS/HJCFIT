Searching for Roots
===================

The default procedure for finding roots is a three step process:

1. Searching for sensible upper and lower bounds bracketing all roots
2. Bisecting the bracket above to obtain bracket with a single root each
3. Using a standard root-finding method to search for the root within that bracket

The first step is carried out by iteratively computing the eigenvalues of :math:`H(s)` and setting
the candidate lower(upper) boundary below(above) the smallest(largest) eigenvalue. Additionally, the
upper boundary is set to a value such that :math:`\mathop{det}H(s)` is strictly positive. There are
two convenience functions to encapsulate this functionality, :py:func:`find_upper_bound_for_roots`
and :py:func:`find_lower_bound_for_roots`. 

:python: 

  .. literalinclude:: ../../code/roots.py
     :language: python
     :lines: 1-19

  .. note:: 

     This package can be compiled to use 128bit reals. However, numpy does not provide all of its
     linear algebra utilities for such type. As consequence, this package exposes some of the
     functionality that it needs for its reals.

:c++11:

  .. literalinclude:: ../../code/roots.cc
     :language: c++
     :lines: 1-33, 59-


The snippets above check that upper and lower bounds are indeed upper and lower bound, as
advertised. It is possible that overflow errors make it difficult to find an upper or lower bound.
It is possible in most function and classes to pass actual values for the upper and lower bounds
(rather than the default, ``None``).

The second step of the process is to bisect the input bracket until intervals are found which
contain a single root (e.g. a single eigenvalue of H(s)).

:python: 

  .. literalinclude:: ../../code/roots.py
     :language: python
     :lines: 31

:c++11:

  .. literalinclude:: ../../code/roots.cc
     :language: c++
     :lines: 44, 45

The third step is performed by calling (by default) the :py:func:`brentq` subroutine. This method is
copied straight from Scipy, with some modifications to allow it to cope with 128bit reals, if need
be.

:python: 

  .. literalinclude:: ../../code/roots.py
     :language: python
     :lines: 32-35

:c++11:

  .. literalinclude:: ../../code/roots.cc
     :language: c++
     :lines: 48-52

The whole procedure is encapsulated within the function :py:func:`find_roots`. On top of the
parameters in the snippets below, :py:func:`find_roots` can take variety of parameters to control
the root-finding procedure. Most notably, it accepts ``lower_bound`` and ``upper_bound`` keywords,
allowing users to by-pass the first step if need be.

:python: 

  .. literalinclude:: ../../code/roots.py
     :language: python
     :lines: 38-39

:c++11:

  .. literalinclude:: ../../code/roots.cc
     :language: c++
     :lines: 55-58
