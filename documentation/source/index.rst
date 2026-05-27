.. HJCFIT documentation master file, created by
   sphinx-quickstart on Wed Jul 31 17:46:57 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###################################
Welcome to HJCFIT's documentation!
###################################

HJCFIT provides full maximum likelihood fitting of a kinetic mechanism directly to the entire sequence of open and shut times, with exact missed events correction.
The package is derived from the DCPROGS_ suite and consists of a C++ implementation of the Likelihood
calculations along with Python wrappers. 

The name of the program is an acronym for Hawkes, Jalali & Colquhoun, whose papers in 1990 and 1992 (HJC92) described the exact solution of the missed event problem, which is the basis of the program. The HJCFIT method was first described by Colquhoun, Hawkes & Srodzinski in 1996 (CHS96).

For a description  of the methods involved, see :cite:`colquhoun:1982`, :cite:`hawkes:1992`,
:cite:`colquhoun:1995a`, :cite:`colquhoun:1995b`, :cite:`colquhoun:1996`.

Contents:

.. toctree::
   :maxdepth: 2

   manual.rst

   api.rst

   examples.rst

   parallel.rst

   mpfr.rst

   install.rst

Indices and tables
==================


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. bibliography:: bibliography.bib

.. _DCPROGS: http://www.ucl.ac.uk/Pharmacology/dcpr95.html
