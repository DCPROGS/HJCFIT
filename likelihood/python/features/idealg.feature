Feature: Check IdealG bindings

  IdealG is a simple wrapper around a transition matrix. It gives the ability to compute the ideal
  likelihood. 

  Scenario: Computation of af, fa
    Given a list of 10 random q-matrices
      And a list of 5 random times between 0 and 8*1e-4
      And a parameter tolerance=1e-7
     When IdealG objects are instantiated with the q-matrices
     Then computing af for each time yields exp(t Q_AA) Q_AF
      And computing fa for each time yields exp(t Q_FF) Q_FA

  Scenario: Computation of laplace_af, laplace_fa
    Given a list of 10 random q-matrices
      And a list of 5 random scales between 0 and 1
      And a parameter tolerance=1e-7
     When IdealG objects are instantiated with the matrices and nopens
     Then computing laplace_af for each scale yields (sI - Q_AA)^-1 Q_AF
      And computing laplace_fa for each scale yields (sI - Q_FF)^-1 Q_FA
      
  Scenario Outline: Computation of equilibrium vectors
    Given a list of 100 random ideal likelihoods
      And a parameter tolerance=1e-7
     When the <name> equilibrium occupancies are computed 
     Then the <name> equilibrium occupancies are the only solution to the equilibrium equations
      And the components of the <name> equilibrium occupancies sum to one

    Examples:
 
      | name    |
      | initial |
      | final   |
