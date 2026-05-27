.. _cpp_exceptions_api:

Exceptions
----------

Exceptions are located in the file ``likelihood/errors.h``. 

.. doxygenclass:: HJCFIT::errors::Root


General
+++++++

.. doxygenclass:: HJCFIT::errors::Index
.. doxygenclass:: HJCFIT::errors::Runtime
.. doxygenclass:: HJCFIT::errors::NotImplemented


Math
++++

.. doxygenclass:: HJCFIT::errors::Math
.. doxygenclass:: HJCFIT::errors::Mass
.. doxygenclass:: HJCFIT::errors::ComplexEigenvalues
.. doxygenclass:: HJCFIT::errors::NaN
.. doxygenclass:: HJCFIT::errors::Domain
.. doxygenclass:: HJCFIT::errors::MaxIterations
.. doxygenclass:: HJCFIT::errors::NotInvertible


Python 
++++++

.. doxygenclass:: HJCFIT::errors::Python
.. doxygenclass:: HJCFIT::errors::PythonErrorAlreadyThrown
.. doxygenclass:: HJCFIT::errors::PythonTypeError
.. doxygenclass:: HJCFIT::errors::PythonValueError
