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

#. The library is in the ``PATH`` (windows), ``DYLD_LIBRARY_PATH`` (Mac), 
   or the ``LD_LIBRARY_PATH`` (Linux)
#. The python bindings are in the ``sys.path`` 
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

Including Jupyter Notebooks in documentation
============================================

``make documentation`` will attempt to convert and include all Jupyter 
notebooks in ``exploration`` to be included in the documentation. This requires
``pandoc`` and ``nbconvert`` to be installed. If these dependencies are not
found the build will skip this step. ```nbconvert`` can be installed from your
regular python package manager (conda or pip). ``pandoc`` can be installed with
os packages managers. On OSX:

.. code-block:: bash
  
  brew install pandoc

on Linux:

.. code-block:: bash

  apt-get install pandoc
  
as well as from conda.

``make documentation`` will also attempt to execute the notebooks before 
including them. This requires ``ipykernel`` to be installed which can be
installed in the same way as ``nbconvert`` above. In addition it requires 
all dependencies of the notebooks to be installed. These are currently
``matplotlib`` and ``networkx`` 

You can also control the Jupyter kernel used to execute the notebooks by
setting the CMake variable ``JUPYTER_KERNEL`` This defaults to ``python2`` or
``python3`` depending on the python version used but can be useful to override
if you need to use a custom environment.

The execution of notebooks can also explicitly be disabled by setting the 
variable ``executenotebooks`` to off.

Updating the web-page
=====================

The data for the web page resides on the same git repository that the code does 
in a special branch called ``gh-pages``. And conversely, github knows to 
render `here <http://dcprogs.github.io/HJCFITdevdocs/>`__. anything that is in 
the branch ``gh-pages``. 

It is possible to update the data and the web-page with the following commands:

#. Commit any changes to the code that should be kept safe.
#. Go to the build directory
#. Update the docs

.. code-block:: bash

  make documentation 


#. Checkout the gh_pages using one the two lines below:

.. code-block:: bash

    git checkout -t origin/gh-pages # If first time, if the branch does not exist 
    git checkout gh-pages 


At this point, the source directory does not contain code anymore. It contains data for the documentation webpage.

1. Copy the new documentation from the build directory to the source directory:

.. code-block:: bash

  rsync -r documentation/sphinx/* ..

1. Commit the changes to the documentation. If nothing happens, 
   there were likely no changes:

.. code-block:: bash

  git commit -a

At this juncture, the data has been updated on the local computer. All that 
needs to be done is to push it to github, so github knows to render it.

1. Push the changes back to github so the web-site can be updated:

.. code-block:: bash
  
  git push
      
1. Checkout the master branch again

.. code-block:: bash
  
  git checkout master

Compiling the documentation without Python bindings
===================================================

A fairly bare documentation of the c++ api is available. 
It can be obtained by running

.. code-block:: bash

  cd /path/to/build/
  make doxydocs


The documents are then available at ``build/documentation/html/index.html``.
