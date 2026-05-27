**********
User Guide
**********

The likelihood :math:`L` of a sequence of observed open and shut events :math:`\{t_{oi}, t_{si}\}`
can be computed as a product of open-shut,
:math:`{}^e\mathcal{G}_{AF}(t_{oi})` and shut-open transition densities,
:math:`{}^e\mathcal{G}_{FA}(t_{oi})`, both accounting for missed events :cite:`colquhoun:1996`: 

.. _log10likelihood_equation:
.. math::

  L(Q, \{t_{oi}, t_{si}\})
    = \phi_A 
      \left[\prod_i^{n-1}{}^e\mathcal{G}_{AF}(t_{oi}){}^e\mathcal{G}_{AF}(t_{si})\right]
      {}^e\mathcal{G}_{AF}(t_{on}) \phi_e

where :math:`Q` is the transition rate matrix, and :math:`\phi_A` and :math:`\phi_e` are the initial
and final equilibrium vectors. All these objects -- as well as their components -- can be accessed both from
Python and from C++. In the following, we try and show how they can be created and manipulated from
either language.


.. toctree::
   :maxdepth: 2


   manual/qmatrix.rst
   manual/likelihood_of_Q.rst
   manual/missedeventsG.rst
   manual/idealG.rst
   manual/vectors.rst
   manual/exact_survivor.rst
   manual/approx_survivor.rst
   manual/determinant_equation.rst
   manual/roots.rst

How to read this manual
-----------------------

Each topic below is illustrated by an example in C++, and another in Python. These examples can be
found in the source code of the package within the directory ``documentation/code``. In any case,
one can copy and paste.

The user guide should read using the left eye only, and keeping the right eye on the api, or
vice-versa. Most of the C++ functions and Python bindings offer more functionality or parameters
than mentioned in the manual. 
