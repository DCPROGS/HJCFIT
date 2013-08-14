.. _cpp_exceptions_api:

Exceptions
----------

Exceptions are located in the file ``likelihood/errors.h``. 

.. doxygenclass:: DCProgs::errors::Root


General
+++++++

.. doxygenclass:: DCProgs::errors::Index
.. doxygenclass:: DCProgs::errors::Runtime
.. doxygenclass:: DCProgs::errors::NotImplemented


Math
++++

.. doxygenclass:: DCProgs::errors::Math
.. doxygenclass:: DCProgs::errors::Mass
.. doxygenclass:: DCProgs::errors::ComplexEigenvalues
.. doxygenclass:: DCProgs::errors::NaN
.. doxygenclass:: DCProgs::errors::Domain
.. doxygenclass:: DCProgs::errors::MaxIterations
.. doxygenclass:: DCProgs::errors::NotInvertible


Python 
++++++

.. doxygenclass:: DCProgs::errors::Python
.. doxygenclass:: DCProgs::errors::PythonErrorAlreadyThrown
.. doxygenclass:: DCProgs::errors::PythonTypeError
.. doxygenclass:: DCProgs::errors::PythonValueError
