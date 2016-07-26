******************************
Customising build and install:
******************************

Customising Installation location
=================================

CMake does provide a gui and a ncurse interface. They are not covered here,
but highly recommended nevertheless.


Installing to a particular root directory from the command-line with CMake is
fairly easy:

.. code-block:: bash

  cd /path/to/dcprogs_source/build
  cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/to
  make
  make install

The above will put executable in ``/path/to/install/to/bin``, headers in
``/path/to/install/to/include``, and libraries in ```/path/to/install/to/lib``.

Specific Eigen Installation
===========================

Similarly, compiling with a particular installation of
`Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ in mind can be
done with:

.. code-block:: bash

  cd /path/to/dcprogs_source/build
  cmake .. -DEIGEN3_INCLUDE_DIR=/path/to/include/eigen3


On windows, the path is likely something like `\path\to\eigen\installation\include\eigen3`.


Customising Compiler
====================


It is also possible to set the compiler explicitly. However, this *must* be done at the very
start. If CMake was run in the build directory, then everything in that directory should be
deleted[^delete the build, not the source directory!] before attempting to set the compiler.

.. code-block:: bash

  cd /path/to/dcprogs_source
  mkdir build && build
  export CC=/path/to/ccompiler
  export CXX=/path/to/cppcompiler
  cmake .. -DCMAKE_CXX_COMPILER=/path/to/compiler
  make

Enable Long Double numbers
==========================

Finally, it is possible to compile the code to run with long double numbers
(long doubles are 80bit on most platforms appart from MS Visual Studio where
they are identical to doubles i.e. 64bit). This could alleviate some of the
overflow/underflow errors at the cost of performance. It is not a solution but
a step in the right direction.

.. code-block:: bash

  cd /path/to/dcprogs/build
  cmake .. -DDCPROGS_LONG_DOUBLE=TRUE

At this juncture, functions that return python scalars are still returning real
numbers of 64bit. Functions that return numpy arrays have the correct size, however.


Disable OpenMP
===============

OpenMP support is automatically enabled provided that your compiler supports it.
You can explicitly disable it by doing:

.. code-block:: bash

  cd /path/to/dcprogs/build
  cmake .. -Dopenmp=off


Python Bindings
===============

The Python bindings are automatically enabled but can be disabled by doing:

.. code-block:: bash

  cd /path/to/dcprogs/build
  cmake .. -DpythonBindings=off

Enabling fallback to Multi precision arithmetic
===============================================

Todo.
