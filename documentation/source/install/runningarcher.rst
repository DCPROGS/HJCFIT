.. _runningonarcher:

*************************
Running HJCFIT on Archer:
*************************

There's good documentation about ARCHER on `their website
<http://www.archer.ac.uk>`__, but here's an extract of what is useful for
running HJCFIT on Archer.

Filesystems
===========

Once logged in ARCHER, there are two different filesystems we need to be aware of:

* ``/home``: This is the place were we land when logging in ARCHER. It is **only 
  visible from the login nodes and the MOM nodes** (the ones that can run ``aprun``
  and send tasks to the compute nodes). You should use this space for storing job 
  files (PBS files), outputs and the like.
* ``/work``: This is **visible from everywhere**: login, MOM and compute nodes. 
  You should use this for storing executable files that'll be run by the compute 
  nodes and other files that the compute nodes need to have access to when running 
  those executables.

Note we have a shared folder for RSDG and DCProgs team under
``/work/ecse0506/ecse0506/shared`` where we can place test data sets, etc.


.bashrc
=======


To make your life easier, it is a good idea to add in your `.bashrc` two lines 
to create variables to move around the filesystem:

.. code-block:: bash

   export WORK=/work/ecse0506/ecse0506/$USER
   export SHARED=/work/ecse0506/ecse0506/shared

You can also configure aliases, like:
 
.. code-block:: bash

   alias qu="qstat -u $USER"
 
And make the system put you on the $WORK filesystem from start, load a virtual 
environment or source a `script <https://github.com/DCPROGS/HJCFIT/blob/develop/utils/archer/loadArcherModules.sh>`__
to load all necessary modules to work on HJCFIT, so you don't have to do it everytime
you log into ARCHER. See a sample `.bashrc` `here <https://github.com/DCPROGS/HJCFIT/blob/develop/utils/archer/sample_bashrc>`__.

Python Virtual Environment
==========================

To work with Python on ARCHER, we are using a virtual environment, which is the strategy recommended by ARCHER. 

To create it, you can run `this script <https://github.com/DCPROGS/HJCFIT/blob/develop/utils/archer/createVirtualEnv.sh>`__
that will install all the necessary packages to run HJCFIT in a virtual environment called ``dcprogs``. 

You will also need to install any extra packages or projects you need, for example
to work with DCPYPS, you'll need to clone it and then install it:

.. code-block:: bash

   cd $WORK
   pip install git+https://github.com/DCPROGS/DCPYPS.git

Once the virtual environment is ready, you can activate or deactivate it with:

.. code-block:: bash

   source activate dcprogs
   source deactivate dcprogs


Loging in to ARCHER and getting HJCFIT
======================================


Since we all have SAFE accounts configured for the project, we just need to do:

.. code-block:: bash

   ssh $USER@login.archer.ac.uk

To get the code, build it and test it in the login node, do as usual for Unix systems:

.. code-block:: bash

   git clone -b develop https://github.com/DCPROGS/HJCFIT.git
   cd HJCFIT
   mkdir build
   cd build
   cmake ..
   make install
   make test

Note that for this work, your ``.bashrc`` should have loaded all the necessary modules. 
See the ``Environment`` section.

Job files
=========


Job scripts are written on a PBS file and follow a specific structure. Here's a sample job script that runs a hello world bash script :

.. code-block:: bash
   
   #!/bin/bash --login

   #PBS -N hello_archer
   #PBS -l select=1
   #PBS -l walltime=0:0:30
   #PBS -A $BUDGET
   
   # This shifts to the $WORK directory
   cd $WORK
   
   aprun -n 24 ./scripts/hello_archer.sh

Job files need a few parameters to be set in the header of the PBS script:

* ``-N <string>```: Specifies job name
* ``-l select=<number>```: Number of nodes needed
* ``-l walltime=<hours>:<minutes>:<seconds>``: time requested for the job
* ``-A <project_code>``: budget code from where the used time will be subtracted 

The body of the job script ``cd`` to the ``$WORK`` folder where we have our 
executables and then uses ``aprun`` to execute the script in parallel using 24 nodes.


Submitting a job
================

To submit a job in the queue, you can do this:

.. code-block:: bash

   qsub myjobfile.pbs

This will submit it to the general queue, and you can check its status with 

.. code-block:: bash

   qstat -u $USER

Or 

.. code-block:: bash

   checkQueue


To delete a job you have submitted:

.. code-block:: bash

   qdel <job_ID_seen_in_queue>

Use ``man qsub``, etc., for more info.

Queues
======

The **standard queue** takes sometimes too long for jobs to be run.

From 9am to 5pm, Monday to Friday, there is a **short queue** available to run
interactive jobs. You will land on a MOM node once you've launched the job, and
will be able to run ``aprun`` that'll trigger tasks in the compute nodes. This is
very handy for short tests for example when testing project configuration. Jobs
in this queue are restricted to 20 minutes walltime and a maximum of 8 nodes.
This is the command you need to run:

.. code-block:: bash

   qsub -q short -IVl select=1,walltime=0:5:0 -A $BUDGET

* ``-q short``: Indicates we don't want to use the standard queue, but the short one. 
* ``-I`` indicates the job is interactive.
* ``-V`` exports the user's environment (I think it runs ~/.bashrc)
* ``-l`` followed by resource list: 
  * ``select=1`` indicates one node will be used
  * ``walltime:0:10:0`` indicates 10 minutes of time available for our job
  * ``-A $BUDGET`` followed by project code indicates the budget the time/resources allocates should go to.

In a short time, you'll be on a ``MOM`` node and sent to your ``$HOME`` folder. 
Remember to cd to ``$WORK`` again, otherwise you can't run ``aprun``.

Once there, you can do things like running a likelihood test in 4 nodes:

.. code-block:: bash

   aprun -n 4 ./HJCFIT/build/likelihood/test_likelihood


More on ``aprun`` `here <http://www.archer.ac.uk/documentation/user-guide/batch.php#sec-5.4.2>`__.

See other kinds of ARCHER queues `here <http://www.archer.ac.uk/documentation/user-guide/batch.php#sec-5.8>`__.


Copying files to and from ARCHER
================================

You might need to copy files to/from ARCHER. This can be done via ``scp``, for example. 
See `ARCHER's documentation <http://www.archer.ac.uk/documentation/transfer/#ARCHER_scp>`__.

Note that if you are on a **Windows** machine and want to run ``scp`` from the command line, 
you can use `chocolatey <https://chocolatey.org>`__ and install it like this:

.. code-block:: bash

   choco install win32-openssh


Virtualenvs on archer.
======================

When running Anaconda python in a virtual env you may see something like.

.. code-block:: none

   python: error while loading shared libraries: libpython3.5m.so.1.0: cannot open shared object file: No such file or directory

This happens because aprun copies python to the compute node. It breaks because
the r path to ``libpython3.5m.so.1.0`` is is defined as ``$ORIGIN/../lib/`` and
the linker resolves ``$ORIGIN`` to the directory where the executable is
installed. You can prevent aprun from coping the executable by passing -b to it
(see the aprun man page) Alternatively you can set LD_LIBRARY_PATH to help
python find the library. 
