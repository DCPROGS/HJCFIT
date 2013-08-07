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

Where :math:`{}^e\mathcal{G}_{AF}` and :math:`{}^e\mathcal{G}_{FA}` are the laplacian of the
missed-events likelihood (or equivilantly, the ideal likelihood) for :math:`s=0`, and :math:`[a]_i`
indicates the :math:`i^{th}` component of vector :math:`a`.
