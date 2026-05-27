Type Hierarchy
--------------

It is often convenient to express all types in a software from a few starting points. This is called
a type hierarchy. The advantage is that it is then fairly easy to switch, say, from using `double` to
`long double`. Indeed, doing so takes all of one line in this package, as can be seen in the file
:file:`HJCFITConfig.h.in` at the root of the package. We describe below the type hierarchy of the
package. All the types are prefixed with ``t_`` to indicate that it is a type. Where instructive, the
names of more complex types will start with ``i``, ``r``, ``c``, ``b`` for integer, real, complex,
and bool.


Basic types
+++++++++++

.. doxygentypedef:: HJCFIT::t_real
.. doxygentypedef:: HJCFIT::t_int
.. doxygentypedef:: HJCFIT::t_uint

Eigen/Math types
++++++++++++++++

.. doxygentypedef:: HJCFIT::t_complex

.. doxygentypedef:: HJCFIT::t_rvector
.. doxygentypedef:: HJCFIT::t_initvec
.. doxygentypedef:: HJCFIT::t_bmatrix

.. doxygentypedef:: HJCFIT::t_rmatrix
.. doxygentypedef:: HJCFIT::t_stack_rmatrix
.. doxygentypedef:: HJCFIT::t_cmatrix

.. doxygentypedef:: HJCFIT::t_Bursts
.. doxygentypedef:: HJCFIT::t_Burst

Multi-precision types
+++++++++++++++++++++
.. unfortunately this will trigger a warning with missing symbols if
.. HJCFIT_USE_MPFR == False since the non included section is parsed anyway
.. see https://github.com/sphinx-doc/sphinx/issues/1635
.. ifconfig:: HJCFIT_USE_MPFR==True

    .. doxygentypedef:: HJCFIT::t_mpfr_real
    .. doxygentypedef:: HJCFIT::t_mpfr_complex
    .. doxygentypedef:: HJCFIT::t_mpfr_cvector
    .. doxygentypedef:: HJCFIT::t_mpfr_rmatrix

.. ifconfig:: HJCFIT_USE_MPFR==False

    HJCFIT is build without Multi-precision support

Global Data
+++++++++++

.. c:var:: HJCFIT::quiet_nan

    Holds an alias to `NaN` representation in :cpp:type:`HJCFIT::t_real`. The code will fail to
    compile if `NaN` cannot be represented in this type.
