Feature: Determinantal Equation bindings

  The following tests the c++ to python bindings of DeterminantEq. 

  Scenario Outline: Initialise the determinantal equation with basic objects
    Given a <matrix>, <nopen>, <tau>, and <doopen>
    When  a determinantal equation is instantiated
    Then  instantiation did not throw
    And   equation's tau is <tau>

    Examples:

      | matrix   | nopen | tau   | doopen  |
      | classic  | 1     |  0.4  |  True   |
      | classic  | 1     |  0.1  |  False  |
      | classic  | 4     |  0.3  |  True   |
      | classic  | 4     |  0.1  |  False  |
      | classic  | 2     |  0.2  |  True   |
      | classic  | 2     |  0.1  |  False  |

  Scenario Outline: Initialise the determinantal equation with a state matrix
    Given a StateMatrix instantiated with <matrix> and <nopen>
    When  a determinantal equation is instantiated from a state matrix, <tau>, and <doopen>
    Then  instantiation did not throw
    And   equation's tau is <tau>

    Examples:

      | matrix               | nopen |  tau  | doopen  |
      | classic              |   1   |  0.5  |  True   |
      | classic              |   1   |  0.1  |  False  |
      | classic              |   4   |  0.1  |  True   |
      | classic              |   4   |  0.6  |  False  |
      | classic              |   2   |  0.1  |  True   |
      | classic              |   2   |  0.7  |  False  |
      | complex eigenvalues  |   2   | 1e-4  |  True   | 

  Scenario Outline: Expected initialization failure with basic objects
    Given a <matrix>, <nopen>, <tau>, and <doopen>
    When  a determinantal equation is instantiated
    Then  instantiation threw <Exception>

    Examples:

      | matrix               | nopen | tau   | doopen | Exception       |
      | classic              |   0   |  0.4  |  True  | ValueError      |
      | classic              |   5   |  0.1  |  False | ValueError      |
      | empty                |   0   |  0.3  |  True  | ValueError      |
      | spam                 |   1   |  0.1  |  False | TypeError       |
      | rectangle            |   1   |  0.1  |  True  | ValueError      |


  Scenario Outline: Computes the determinant for open states of classic case
    Given the transition matrix below with 2 open states
          """
              -3050,        50,  3000,      0,    0
              2./3., -1502./3.,     0,    500,    0
                 15,         0, -2065,     50, 2000
                  0,     15000,  4000, -19000,    0
                  0,         0,    10,      0,  -10
          """
          
    When  the open-state determinant is computed for s=(<s>) and tau=1e-4
    Then  the result is close to zero (1e-2)

    Examples:

      | s                                                | 
      | -3045.285776037674                               |
      |  -162.92946543451328                             |
      | -3045.285776037674, -162.92946543451328          |
      | [-3045.285776037674, -162.92946543451328]        |
      | array([-3045.285776037674, -162.92946543451328]) |



  Scenario Outline: Computes the determinant for closed states of classic case
    Given the transition matrix below with 2 open states
          """
              -3050,        50,  3000,      0,    0
              2./3., -1502./3.,     0,    500,    0
                 15,         0, -2065,     50, 2000
                  0,     15000,  4000, -19000,    0
                  0,         0,    10,      0,  -10
          """
          
    When  the close-state determinant is computed for s=(<s>) and tau=1e-4
    Then  the result is close to zero (1e-2)

    Examples:

      | s                                                                      | 
      | -17090.192769236815                                                    |
      | -2058.0812921673496                                                    |
      | -0.24356535498785126                                                   |
      | ( -17090.1927692368, -2058.08129216735, -0.24356535498785126 )         |
      | [ -17090.1927692368, -2058.08129216735, -0.24356535498785126 ]         |
      | array([ -17090.1927692368, -2058.08129216735, -0.24356535498785126 ] ) |

