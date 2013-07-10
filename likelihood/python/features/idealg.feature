Feature: Check IdealG bindings

  IdealG is a simple wrapper around a transition matrix. It gives the ability to compute the ideal
  likelihood. 

  Scenario: Computation of af, fa
    Given a list of 10 random state matrices
      And a list of random of 5 times between 0 and 8*1e-4
      And a parameter tolerance=1e-8
     When IdealG objects are instantiated with the state matrices
     Then computing af for each time yields exp(t Q_AA) Q_AF
      And computing fa for each time yields exp(t Q_FF) Q_FA

  Scenario: Computation of laplace_af, laplace_fa
    Given a list of 10 random state matrices
      And a list of random of 5 scales between 0 and 1
      And a parameter tolerance=1e-8
     When IdealG objects are instantiated with the matrices and nopens
     Then computing laplace_af for each scale yields (sI - Q_AA)^-1 Q_AF
      And computing laplace_fa for each scale yields (sI - Q_FF)^-1 Q_FA
      
  Scenario: Computation of equilibrium vectors
    Given a list of 500 random state matrices
      And a parameter tolerance=1e-6
     When IdealG objects are instantiated with the state matrices
     Then the initial occupancies exists and is the kernel of I - laplace_af * laplace_fa
      And the final occupancies exists and is the kernel of I - laplace_fa * laplace_af
