Type Hierarchy
--------------

It is often convenient to express all types in a software from a few starting points. This is called
a type hierarchy. The advantage is that it is then fairly easy to swich, say, from using `double` to
`long double`. Indeed, doing so takes all of one line in this package, as can be seen in the file
:file:`DCProgsConfig.h.in` at the root of the package. We describe below the type hierarchy of the
package. All the types are prefixe with ``t_`` to indicate that it is a type. Where instructive, the
names of more complex types will start with ``i``, ``r``, ``c``, ``b`` for integer, real, complex,
and bool.


Basic types
+++++++++++

.. doxygentypedef:: DCProgs::t_real
.. doxygentypedef:: DCProgs::t_int
.. doxygentypedef:: DCProgs::t_uint

Eigen/Math types
++++++++++++++++

.. doxygentypedef:: DCProgs::t_complex

.. doxygentypedef:: DCProgs::t_rvector
.. doxygentypedef:: DCProgs::t_initvec
.. doxygentypedef:: DCProgs::t_bmatrix

.. doxygentypedef:: DCProgs::t_rmatrix
.. doxygentypedef:: DCProgs::t_cmatrix

.. doxygentypedef:: DCProgs::t_Bursts
.. doxygentypedef:: DCProgs::t_Burst


Global Data
+++++++++++

.. cpp:var: DCProgs::quiet_nan

   Holds an alias to `NaN` representation in :cpp:type:`DCProgs::t_real`. The code will fail to
   compile if `NaN` cannot be represented in this type.
