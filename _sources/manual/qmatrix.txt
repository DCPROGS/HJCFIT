.. _manual_qmatrix:

The :math:`Q`-matrix
====================

The :math:`Q`-matrix represents the most basic data structure of this library. It holds the rates
between different states in the mechanism, for open and shut states. In practice, it can be
represented by a matrix of rates and the number of open states. As such, the initialization takes
both of these input:

:python: 

  The input array must be square. It can be any type that is convertible to a numpy array of reals.
  
  .. literalinclude:: ../../code/qmatrix.py
     :language: python
     :lines: 1-12
  
:c++11:

  .. literalinclude:: ../../code/qmatrix.cc
     :language: c++
     :lines: 1-19, 40-


In both cases, the open-open transitions *must* be located in the top-left corner of the input
matrix. 

:cpp:class:`QMatrix` is a structure instance, without private members. Hence both the matrix and
number of open states can be accessed directly. In the following, we show access to `nopen` and to
the rates. Among other things, the code below verifies that the rates of the CH82 model do sum to
zero over rows (this condition is not enforced by the :cpp:class:`QMatrix` object). 

:python: 

  The matrix is a numpy array. As such, it can be used within any of the efficient numerical
  algorithms provided. There are some restrictions if this package is compiled with 128 bit reals.

  .. literalinclude:: ../../code/qmatrix.py
     :language: python
     :lines: 14-26

:c++11:

  The matrix is an eigen_ array. Many linear algebra operations can be performed very efficiently
  using Eigen_'s interface.

  .. literalinclude:: ../../code/qmatrix.cc
     :language: c++
     :lines: 21-33


Finally, the blocks of the  :math:`Q` can be accessed efficiently via functions (or properties in
python). These are wrappers around the rate matrix. These wrappers are numerically efficient. They
can and *should* be used within numerical algorithms.

:python: 

  .. literalinclude:: ../../code/qmatrix.py
     :language: python
     :lines: 28-

:c++11:

  .. literalinclude:: ../../code/qmatrix.cc
     :language: c++
     :lines: 35-39

