Approximate Survivor Function :math:`R_A(t)`
============================================


The approximate survivor function is defined, and used, for :math:`t > 2\tau`:

.. math::
  \mathcal{R}_A(u) \backsim
    \sum_{i=1}^{k_A}\frac{\vec{c}_i\vec{r}_i}{\vec{r}_iW'(s_i)\vec{c}_i}  e^{u s_i}

where :math:`c_i` and :math:`r_i` are the column and row eigenvectors of a
:math:`H(s_i)` (defined below), :math:`s_i` are the roots of :math:`\mathop{det}W(s)`.

.. math::
   W(s) = sI - H(s)\\

.. math::
   H(s) = \mathcal{Q}_{AA} + \mathcal{Q}_{AF}\left[sI-\mathcal{Q}_{FF}\right]^{-1}
          \left(I-e^{-(sI-\mathcal{Q}_{FF})\tau}\right)\mathcal{Q}_{FA}

.. math::

   W'(s) = I+\mathcal{Q}_{AF}\left[\bar{S}_{FF}(s)\left(sI-\mathcal{Q}_{FF}\right)^{-1}
      -\tau \left(I -\bar{\mathbf{S}}_{FF}(s)\right)\right]\bar{\mathcal{G}}_{FA}(s)

.. math::

   \bar{\mathbf{S}}_{FF}(s) = \int_0^\tau\!\mathrm{d}\,t\ e^{-st}e^{\mathcal{Q}_{FF}t}

.. math::
   \bar{\mathcal{G}}_{FA}(s) = \left[sI-\mathcal{Q}_{FF}\right]^{-1}\mathcal{Q}_{FA}




The approximate survivor function can be intialized from a :math:`\mathcal{Q}`-matrix and the resolution
:math:`\tau`:

:python: 

  .. literalinclude:: ../../code/approx_survivor.py
     :language: python
     :lines: 1-13

:c++11:

  .. literalinclude:: ../../code/approx_survivor.cc
     :language: c++
     :lines: 1-19, 36-

.. note::
   
   The initialization shown above hides many details of the implementation, most notably how the
   roots are obtained. However, this object can also be initialed in exactly the same way as
   :ref:`MissedEventsG <manual_eG>`, with the same sets parameters (except for ``nmax``).

The open and shut time survivor likelihood can be computed using a single call:

:python: 

  The python bindings accept both scalars and array inputs.

  .. literalinclude:: ../../code/approx_survivor.py
     :language: python
     :lines: 15-19


:c++11:

  .. literalinclude:: ../../code/approx_survivor.cc
     :language: c++
     :lines: 21-30


The coefficient and the exponents of the exponentials that make up the asymptotic expression are
exposed as shown below. 

:python:

  .. literalinclude:: ../../code/approx_survivor.py
     :language: python
     :lines: 23-

:c++11:

  .. literalinclude:: ../../code/approx_survivor.cc
     :language: c++
     :lines: 32-35
