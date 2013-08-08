Feature: Check Missed Events G functionality

  The missed events G objects should contain most of everything we need to compute likelihoods over
  time series.

  Scenario Outline: Computation of af, fa, compared to exact and approx survivors
    Given a list of 20 random q-matrices
      And a list of 10 random times between 1e-4 and 3e-4
      And a list of 10 random times between 0e-4 and 1e-4
      And a list of 10 random times between 3e-4 and 4e-4
      And a parameter tolerance=1e-8
     When MissedEventsG objects are instantiated with the q-matrices and tau=1e-4 and nmax=3
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
    Given a list of 100 random missed-events likelihoods with tau=1e-4 and nmax=3
      And a parameter tolerance=1e-8
     When the <name> equilibrium occupancies are computed 
     Then the <name> equilibrium occupancies are the only solution to the equilibrium equations
      And the components of the <name> equilibrium occupancies sum to one
 
    Examples:
 
      | name    |
      | initial |
      | final   |




  Scenario Outline: Computation of CHS vectors
    Given a list of 10 random missed-events likelihoods with tau=1e-4 and nmax=3
      And a list of 10 random times between 1e-3 and 3e-3
      And a parameter tolerance=1e-6
     When the <name> CHS occupancies are computed 
     Then the <name> CHS occupancies are the solutions to the CHS equations
 
    Examples:
 
      | name    |
      | initial |
      | final   |




  Scenario Outline: Test CHS vector for different models against prior calculation.
    Given the <model> missed-events likelihood
      And a parameter tolerance=1e-10
     When we compute the <which> CHS occupancies with tcrit=<tcrit>
     Then the <which> CHS occupancies compare to <prior> 

    Examples: 

      |  model   |  which  | tcrit | prior                                                     |
      | classic  | initial | 4e-3  | [0.220418113246138, 0.779581886753862]                    |
      | classic  | initial | 1e-3  | [0.220339411420015, 0.779660588579985]                    |
      | classic  | initial | 5e-4  | [0.220101271499783, 0.779898728500217]                    |
      | classic  | initial | 5e-5  | [0.1122948887589, 0.8877051112411]                        |
      | classic  |  final  | 4e-3  | [0.974851711507551, 0.213460487245399, 0.999179485194199] |
      | classic  |  final  | 5e-5  | [0.99960557806125, 1.96870728834576, 0.999999979357141]   |
      |  CH82    | initial |   4   | [0.17394315362718,  0.82605684637282]                     |
      |  CH82    |  final  |   4   | [0.976491211386193, 0.222305380522348, 0.999257244552633] |
      |   CKS    | initial |   4   | [1]                                                       |
      |   CKS    |  final  |   4   | [0.36908082444647, 0.942440306684466]                     |
      |   CB     | initial |   4   | [1]                                                       |
      |   CB     |  final  |   4   | [0.846530054887709, 0.168045183806247, 0.852959014045751] |
