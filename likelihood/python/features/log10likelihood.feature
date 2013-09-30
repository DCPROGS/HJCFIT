Feature: Check log10likelihood

  The log10-likelihood object returns the likelihood of a list of bursts as a
  function of an input matrix. It should be able to handle numbers with
  small/large exponents without producing NaN or infinities (or not too often...).

  Scenario Outline: Compute very small likelihoods
    Given a list of bursts <bursts>
      And a dictionary of parameters
          """ 
          {'itermax': 100, 'lower_bound': -1000000.0,
           'nmax': 2, 'nopen': 3, 'rtol': 1e-12, 'tau': 4e-05, 
           'tcritical': 0.01, 'upper_bound': 0, 'xtol': 1e-12}
          """
     When The log10-likelihood of the classic matrix is computed
     Then the result is defined and finite

    Examples:

      |   bursts                            |
      | [[2e-2]*501]                        |
      | [[2e-4]*3]                          |
      | [[1e-2]*501]                        |
      | [[1e-2]*4 + [2e-4]*501, [1e-2]*501] |
