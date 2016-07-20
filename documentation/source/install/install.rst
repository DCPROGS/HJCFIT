*******************************
Building and installing HJCFIT:
*******************************


Compiling DCProgs
=================

A couple of desing decisions affect the compilation of DCProgs.

  * `c++11 <http://en.wikipedia.org/wiki/C%2B%2B11>`_ is the new standard for 
    the C++ programming languages. It is almost fully implemented by modern 
    (2013) compilers. However, access to c++11 is now always default, and not 
    always straight-forward. However, c++11 introduces a number of features that 
    simplifies programming (e.g. `move semantics <http://www.cprogramming.com/c++11/rvalue-references-and-move-semantics-in-c++11.html>`_)
    greatly. This is a forward looking solution implying some temporary hassle.
  * [GTest](https://code.google.com/p/googletest/) is the c++ unit-test 
    framework from google. It is required when running DCProgs' unit tests only.
    However, `GTest <https://code.google.com/p/googletest/>`_ must be compiled 
    by the code it is testing. This means it should be shipped with DCProgs, 
    or it should be downloaded automatically by the compilation tools. This is
    the option we have chosen. When compiling tests,
    `CMake <http://www.cmake.org/>`_ will automatically download and compile
    `GTest`_
  * The math is done using `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_,
    an efficient and widely used C++ numerical library. 

Dependencies
------------

#. Modern C++ compiler: 
    * g++-4.6 or greater on any platform,
    * clang with g++-4.2 on Mountain Lion
    * intel 13 with gnu standard libraries from g++ 4.6 (intel provides older 
      gnu standard libraries).
#. `CMake`_
   It is available on any platform, either directly through its website or via 
   a packaging system (yum, apt-get for Linux, `Homebrew <http://brew.sh/>`_ 
   for Mac, `chocolatey <http://chocolatey.org/>`_ for windows)
#. `Eigen`_: Users can either install eigen, or install ``wget``,  ``curl``
   or `mercurial <http://mercurial.selenic.com/>`_ and let the build process 
   download eigen.
   
To compile the python bindings for HJCFIT a few additional dependencies are
needed.

#. A working Python installation. 

Multiple different ways of installing python exist. In general we recommend 
`Anaconda <https://www.continuum.io/downloads>`_ but alternatives should work
as well.  In any case Python along with ``numpy`` and ``scipy`` should be 
installed. HJCFIT supports both Python 2.7 and Python 3

#. `SWIG <http://www.swig.org/>`_ used to generate the wrappings beween C++ and
   Python.
#. Behave. A behavior driven developement framework for Python.

Compilation on Mac:
===================

Make sure that you have the needed dependencies installed. We recomend using
`Homebrew`_ to manage depedencies on OSX. Follow the Homebrew documentation to
install it. The instructions below assumes that you are running a working 
python installation with ``Numpy`` and ```Scipy``. Behave is not currently in
Anaconda and should be installed using ``pip`` in any case.

Installing depdendencies using Homebrew:

.. code-block:: bash

  brew install swig cmake
  pip install behave

Then configure and build the code:

.. code-block:: bash
  
  cd /path/to/DCProgs
  mkdir build && cd build
  cmake ..
  make
  make test

Assuming everything works as expected we can now install HJCFIT

  .. code-block:: bash

  make install

Compilation on Legion:
======================

For any compiler, do:

.. code-block:: bash

   cd /path/to/DCProgs
   mkdir build && cd build


Then, if the preferred compiler is g++, the module `compilers/gnu/4.6.3` 
should first be loaded. A working python version should be loaded. Currently
the Legion python modules does not have behave installed so it's recommended
to install it using ``pip --user``. We also need to module load ``SWIG``.

.. code-block:: bash

  pip install --user behave
  module load python3/recommended
  module load swig/3.0.7/gnu-4.9.2
  module swap compilers compilers/gnu/4.9.2

Then:

.. code-block:: bash

  cmake ..
  make
  make test

Assuming everything works as expected we can now install HJCFIT

  .. code-block:: bash

  make install

Compilation on Windows:
=======================


Several different ways of building and installing on Windows exist. It should
be possible to build the code with both MS Visual Studio and MinGW. Currently 
we recommend building using MS Visual Studio 2015. The free `Community edition 
of Visual Studio <https://www.visualstudio.com/en-us/products/visual-studio-community-vs.aspx>`_
is sufficient to build HJCFIT. Note that older versions of Visual Studio did not
ship 64 bit compilers in the free version. This is no longer an issue with the 
2015 version. Python 3.5 in normally build with Visual Studio 2015 where as 
older versions are build with older versions of Visual Studio so to reduce
any issues it is recommended to use Visual Studio 2015 and Python 3.5.

Visual Studio 2015:
-------------------

First ensure Visual Studio is installed. Make sure to select the C++ components
during the installation.

You then need to install the dependencies Swig, CMake. You can install curl from
Anaconda to enable automatic download of Eigen. It's recommended to install 
CMake and Swig from their respective homepages. Make sure that you select add
to path when installing CMake. Following this open a command prompt.

To put the relevant Microsoft compilers on Path you should run the relevant
bat sript. On most systems it should be something like:

.. code-block:: bash

  "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" amd64

You can verify that Visual Studio is correctly loaded by execution ``cl`` as


To install Eigen we need wget, curl or mercurial. Curl can be installed directly
from conda. To run the Python tests we need to install behave.
  
.. code-block:: bash

  conda install curl
  pip install behave


We can now build the code and run the tests.

.. code-block:: bash

  cd /path/to/HJCFIT
  mkdir build && cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release -G "NMake Makefiles" 
  nmake
  nmake test

Assuming everything works as expected we can now install HJCFIT

.. code-block:: bash

  nmake install
  
