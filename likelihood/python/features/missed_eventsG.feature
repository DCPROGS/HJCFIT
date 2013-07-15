Feature: Check Missed Events G functionality

  The missed events G objects should contain most of everything we need to compute likelihoods over
  time series.

  Scenario Outline: Computation of af, fa, compared to exact and approx survivors
    Given a list of 20 random q-matrices
      And a list of 10 random times between 1e-4 and 3e-4
      And a list of 10 random times between 0e-4 and 1e-4
      And a list of 10 random times between 3e-4 and 4e-4
      And a parameter tolerance=1e-8
     When MissedEventsG objects are instantiated with the q-matrices and tau=1e-4 and nmax=2
      And ExactSurvirvor objects are instantiated with the q-matrices and tau=1e-4
      And ApproxSurvivor objects are instantiated with the q-matrices and tau=1e-4
     Then <name> is zero if t is between 0 and 1e-4
      And <name> can be found from ExactSurvivor if t is between 1e-4 and 3e-4
      And <name> can be found from ApproxSurvivor if t is larger than 3e-4

    Examples:

      | name |
      | af   |
      | fa   |




 
  Scenario Outline: Computation of equilibrium occupancy
    Given a list of 100 random missed events with tau=1e-4 and nmax=2
      And a parameter tolerance=1e-8
     When the <name> equilibrium occupancies are computed 
     Then the <name> equilibrium occupancies are the only solution to the equilibrium equations
 
    Examples:
 
      | name    |
      | initial |
      | final   |
