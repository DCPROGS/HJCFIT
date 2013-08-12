Exact Survivor Function :math:`R_A(t)`
======================================


The exact survivor function is defined as:

.. math::
  \mathcal{R}_A(t) = \sum_{m=0}^\infty (-1)^m M_m(t-m\tau)

with,

.. math::
  M_m(t)=\left\{\begin{split}
    0&\quad\text{if}\quad t < 0\\
    \sum_{i=0}^k B_{im}(t) e^{-\lambda_i t}&\quad\text{if}\quad t \geq 0
    &
  \end{split}\right.,\\
  B_{im}(t) = \sum_{r=0}^mC_{imr}t^r,

where the matrices :math:`C_{imr}` are defined as a recursion:

.. math::
  \left\{\begin{eqnarray}
    C_{imm} &=& \frac{1}{m}D_iC_{i(m-1)(m-1)}\\
    C_{im(l\leq m-1)} &=& \frac{1}{l}D_iC_{i(m-1)(l-1)}
        -\sum_{j\neq i}^k
         \sum_{r=l}^{m-1}\frac{r!}{l!(\lambda_i-\lambda_j)^{r-l+1}}D_jC_{i(m-1)r}\\
    C_{im0} &=& \sum_{j\neq i}^k\sum_{r=0}^{m-1}
       \frac{r!}{(\lambda_i-\lambda_j)^{r+1}}
       \left[D_iC_{j(m-1)r}-D_jC_{i(m-1)r}\right]
  \end{eqnarray}\right.

The initial values are :math:`C_{i00} = A_{iAA}`, with :math:`A_i` the
spectral decomposition of the :math:`\mathcal{Q}`-matrix, and :math:`\lambda_i` are its eigenvalues.
Finally, the matrices :math:`D_i` are defined as:

.. math::
  D_i = A_{iAF}e^{\mathcal{Q}_{FF}\tau}\mathcal{Q}_{FA}


.. note::
  This recursion is implemented in the header ``likelihood/recursion.h`` in such a way that it can act
  upon a variety of objects. This makes testing it somewhat easier, since we can defined the
  :math:`D_i`, for instance, as scalars rather than matrices. 


The survivor function can be initialized from a :math:`\mathcal{Q}`-matrix and the resolution
:math:`\tau`:

:python: 

  .. literalinclude:: ../../code/exact_survivor.py
     :language: python
     :lines: 1-13


:c++11:

  .. literalinclude:: ../../code/exact_survivor.cc
     :language: c++
     :lines: 1-19, 39-

The open and shut time survivor likelihood can be computed using a single call:

:python: 

  The python bindings accept both scalars and array inputs.

  .. literalinclude:: ../../code/exact_survivor.py
     :language: python
     :lines: 15-19


:c++11:

  .. literalinclude:: ../../code/exact_survivor.cc
     :language: c++
     :lines: 21-30


The details of the recursions, i.e. the :math:`C_{iml}` matrices, can be accessed directly as shown
below.

:python:

  .. literalinclude:: ../../code/exact_survivor.py
     :language: python
     :lines: 23-

:c++11:

  .. literalinclude:: ../../code/exact_survivor.cc
     :language: c++
     :lines: 32-38
