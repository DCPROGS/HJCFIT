HJCFIT
======

Full maximum likelihood fitting of a mechanism directly to the entire sequence of open and shut times, with exact missed events correction.

The name of the program is an acronym for Hawkes, Jalali & Colquhoun, whose papers in 1990 and 1992 (HJC92) described the exact solution of the missed event problem, which is the basis of the program.  The HJCFIT method was first described by Colquhoun, Hawkes & Srodzinski in 1996 (CHS96). The properties of the estimates of rate constants obtained by this method have now been evaluated (Colquhoun, Hatton & Hawkes, 2003).

The input for HJCFIT is a list of idealised open and shut time intervals.  A kinetic mechanism is specified with some initial guesses for the rate constants. Fitting is done using the Simplex algorithm to maximise the likelihood.


The documentation for this package can be found [here](http://ucl.github.io/dcprogs/index.html).

Explanations on compiling and installing the package can be found in the
[wiki](https://github.com/UCL/dcprogs/wiki).

Relevant references

CH82: Colquhoun D, Hawkes AG (1982) On the stochastic properties of bursts of single ion channel openings and of clusters of bursts. Phil Trans R Soc Lond B 300, 1-59.

HJC92: Hawkes AG, Jalali A, Colquhoun D (1992) Asymptotic distributions of apparent open times and shut times in a single channel record allowing for the omission of brief events. Phil Trans R Soc Lond B 337, 383-404.

CH95a: Colquhoun D, Hawkes AG (1995a) The principles of the stochastic interpretation of ion channel mechanisms. In: Single-channel recording. 2nd ed. (Eds: Sakmann B, Neher E) Plenum Press, New York, pp. 397-482.

CH95b: Colquhoun D, Hawkes AG (1995b) A Q-Matrix Cookbook. In: Single-channel recording. 2nd ed. (Eds: Sakmann B, Neher E) Plenum Press, New York, pp. 589-633.

CHS96: Colquhoun D, Hawkes AG, Srodzinski K (1996) Joint distributions of apparent open and shut times of single-ion channels and maximum likelihood fitting of mechanisms. Phil Trans R Soc Lond A 354, 2555-2590.
