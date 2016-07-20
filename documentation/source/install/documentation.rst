***********************
Building documentation:
***********************

The documentation is written using `doxygen <http://www.doxygen.org>`_ for c++, 
`sphinx <http://sphinx-doc.org/>`_ for python, and 
`breathe <https://pypi.python.org/pypi/breathe>`_ to tie them together. 
It yields `this <http://dcprogs.github.io/HJCFIT/>`_ an automatic build of the 
documentation from the development branch is 
`here <http://dcprogs.github.io/HJCFITdevdocs/>`__.

Pre-requisites
==============

Installing doxygen
------------------

* Linux: apt-get, yum, or other package manager
* Mac: ``brew install doxygen``
* Windows: binaries can be found 
  `here <http://www.stack.nl/~dimitri/doxygen/download.html>`__

Installing Sphinx
-----------------

Sphinx is a regular python package that can be installed as any other python
package.

* Linux: 
    apt-get, yum, pip or conda depending on your setup.
* Mac: ``pip install sphinx`` or ``conda install sphinx``
* Windows: Depending on your setup
    - ``conda.bat install sphinx``
    - or, ``pip install sphinx``

.. warning::
    If using a virtual environment, it is recommended to run 
    ``pip install --upgrade sphinx`` (within the virtual environment) 
    so that sphinx is set up with the right path. Then one should make sure 
    that ``SPHINX_EXECUTABLE`` in ```build/CMakeCache.txt`` points to the 
    correct sphinx executable.

Installing Breathe
------------------

In all cases, do ``pip install Breathe``. If you are using conda it might 
first be necessary to install the pip package via conda.

Installing Sphinx plugin for citation
-------------------------------------

``sphinxcontrib.bibtext`` makes it possible to use a bibtex file, 
``documentation/source/bibliography.bib``, to hold reference to papers. 
It can be installed with ``pip install sphinxcontrib-bibtex``

.. note:: 
  On Mac, with python2.7 installed from brew, I've had to add an empty 
  ``__init__.py`` file in ``site-packages/sphinxcontrib`` for its ``bibtex`` 
  subpackage to be found. This may have been fixed in a later release?

Compiling the Documentation for Python
======================================

The code must be installed and accessible before trying to compile the 
documentation, at least when the python bindings are compiled.
This means:

1. The library is in the ``PATH`` (windows), ``DYLD_LIBRARY_PATH`` (Mac), 
   or the ``LD_LIBRARY_PATH`` (Linux)
1. The python bindings are in the ``sys.path`` 
   (e.g. ``python -c "import dcprogs.likelihood"`` does not fail)

The reason for this is that python documentation software will interrogate 
the package to find out what it contains. Hence the package needs to be 
available and loadable. Please refer to the documentation 
for instructions on installing HJCFIT.

Once HJCFIT is installed and available, do:

.. code-block:: bash
  
  cd /path/to/buid/
  make documentation


The documentation should be available at 
``path/to/build/documentation/sphinx/index.html``.


Compiling the documentation without Python bindings
===================================================

A fairly bare documentation of the c++ api is available. 
It can be obtained by running

.. code-block:: bash

  cd /path/to/build/
  make doxydocs


The documents are then available at ``build/documentation/html/index.html``.
