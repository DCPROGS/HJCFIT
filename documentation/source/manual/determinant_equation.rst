Determinant Equation :math:`W(s) = sI - H(s)`
=============================================


The function :math:`H(s)` is an integral defined as:

.. math::

   H(s) = sI - \mathcal{Q}_{AA} - \mathcal{Q}_{AF}\
     \int_0^\tau e^{-st}e^{\mathcal{Q}_{FF}t}\partial\,t\ \mathcal{Q}_{FA}

It is possible the function :math:`H` as well as its determinant using the
:py:class:`~dcprogs.likelihood.DeterminantEq` objects. This is the object used when solving for the
approximate missed-events likelihood. The determinant equation is initialized in one of two ways,
either from a matrix or :py:class:`~dcprogs.likelihood.QMatrix`.

:python: 

  .. literalinclude:: ../../code/determinanteq.py
     :language: python
     :lines: 3-17


:c++11:

  .. literalinclude:: ../../code/determinanteq.cc
     :language: c++
     :lines: 1-21, 42-


With an object in hand, it is possible to compute :math:`\mathop{det}W(s)` for any :math:`s`. In the
following we demonstrate that the two initialization methods are equivalent and that the determinant
is zero at the roots of :math:`W(s)`, per definition. 

:python: 

  The python bindings accept both scalars and arrays as input. 

  .. literalinclude:: ../../code/determinanteq.py
     :language: python
     :lines: 19-23


:c++11:

  .. literalinclude:: ../../code/determinanteq.cc
     :language: c++
     :lines: 23-31


There exists a convenience function to transform a determinant equation into its "transpose", e.g.
one where A states become F states and F states become A states:

:python: 

  .. literalinclude:: ../../code/determinanteq.py
     :language: python
     :lines: 25-28

  .. note::
     
     Here we choose to create an input which has same internal type as the dcprogs package. This may
     result in faster code since no conversion are required.

:c++11:

  .. literalinclude:: ../../code/determinanteq.cc
     :language: c++
     :lines: 32-36


Finally, it is possible to compute :math:`H(s)` directly, as well as :math:`\frac{\partial
W(s)}{\partial s}`, as demonstrated below.

:python: 

  .. literalinclude:: ../../code/determinanteq.py
     :language: python
     :lines: 30-

:c++11:

  .. literalinclude:: ../../code/determinanteq.cc
     :language: c++
     :lines: 38-42


