Feature: Check that asymptotic expression makes sense.

  The following tests the expression for the asymptotic approximation to the missed event R.
  In order to make things easier and testable, the asymptotic approximations are instantiated with a
  single root at a time. 

  Scenario: Checks that asymptotic expression is made of sum of projection matrix
    Given a list of 500 random determinant equations
     When the roots are computed for each
     Then for each root, applying H - root * I to asymptote at t=0    yields zero (<1e-6)
      And for each root, applying H - root * I to asymptote at t=1e-4 yields zero (<1e-6)
      And for each root, applying H - root * I to asymptote at t=2e-4 yields zero (<1e-6)
      And for each root, applying H - root * I to asymptote at t=3e-4 yields zero (<1e-6)
