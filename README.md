HJCFIT
======

Full maximum likelihood fitting of a mechanism directly to the entire sequence of open and shut times, with exact missed events correction.

The name of the program is an acronym for Hawkes, Jalali & Colquhoun, whose papers in 1990 and 1992 described the exact solution of the missed event problem, which is the basis of the program.  The HJCFIT method was first described by Colquhoun, Hawkes & Srodzinski in 1996. The properties of the estimates of rate constants obtained by this method have now been evaluated (Colquhoun, Hatton & Hawkes, 2003).

The input for HJCFIT is a list of idealised open and shut time intervals.  A kinetic mechanism is specified with some initial guesses for the rate constants. Fitting is done using the Simplex algorithm to maximise the likelihood.


The documentation for this package can be found [here](http://ucl.github.io/dcprogs/index.html).

Explanations on compiling and installing the package can be found in the
[wiki](https://github.com/UCL/dcprogs/wiki).
