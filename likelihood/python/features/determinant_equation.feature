Feature: Determinantal Equation bindings

  The following tests the c++ to python bindings of DeterminantEq. 

  Scenario Outline: Initialise the determinant equation with basic objects
    Given a <matrix>, <nopen>, and <tau>
    When  a determinant equation is instantiated
    Then  instantiation did not throw
    And   equation's tau is <tau>

    Examples:

      | matrix   | nopen | tau   |
      | classic  | 1     |  0.4  |
      | classic  | 1     |  0.1  |
      | classic  | 4     |  0.3  |
      | classic  | 4     |  0.1  |
      | classic  | 2     |  0.2  |
      | classic  | 2     |  0.1  |

  Scenario Outline: Initialise the determinant equation with a state matrix
    Given a QMatrix instantiated with <matrix> and <nopen>
    When  a determinant equation is instantiated from a state matrix and <tau>
    Then  instantiation did not throw
    And   equation's tau is <tau>

    Examples:

      | matrix               | nopen |  tau  |
      | classic              |   1   |  0.5  |
      | classic              |   1   |  0.1  |
      | classic              |   4   |  0.1  |
      | classic              |   4   |  0.6  |
      | classic              |   2   |  0.1  |
      | classic              |   2   |  0.7  |
      | complex eigenvalues  |   2   | 1e-4  |

  Scenario Outline: Expected initialization failure with basic objects
    Given a <matrix>, <nopen>, and <tau>
    When  a determinant equation is instantiated
    Then  instantiation threw <Exception>

    Examples:

      | matrix               | nopen | tau   | Exception               |
      | classic              |   0   |  0.4  | ValueError              |
      | classic              |   5   |  0.1  | ValueError              |
      | empty                |   0   |  0.3  | ValueError              |
      | spam                 |   1   |  0.1  | (ValueError, TypeError) |
      | rectangle            |   1   |  0.1  | ValueError              |


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
          
    When  the shut-state determinant is computed for s=(<s>) and tau=1e-4
    Then  the result is close to zero (1e-2)

    Examples:

      | s                                                                      | 
      | -17090.192769236815                                                    |
      | -2058.0812921673496                                                    |
      | -0.24356535498785126                                                   |
      | ( -17090.1927692368, -2058.08129216735, -0.24356535498785126 )         |
      | [ -17090.1927692368, -2058.08129216735, -0.24356535498785126 ]         |
      | array([ -17090.1927692368, -2058.08129216735, -0.24356535498785126 ] ) |


  Scenario Outline: Computes the H matrix using different input formats
    Given a determinant equation instantiated with the <model> Q matrix
      And a parameter tolerance=1e-8
      And the output matrices below
      """
      q0 = [[ -2.06112149e+03,   5.00016694e+01,   2.00000000e+03],
            [  4.00013355e+03,  -1.82684655e+04,   0.00000000e+00],
            [  1.00000000e+01,   0.00000000e+00,  -1.00000000e+01]]

      q1 = [[ -2.06093128e+03,   5.00017832e+01,   2.00000000e+03],
            [  4.00014266e+03,  -1.82309594e+04,   0.00000000e+00],
            [  1.00000000e+01,   0.00000000e+00,  -1.00000000e+01]]

      q2 = [[  8.17100439e+07,   7.32058151e+04,   2.00000000e+03],
            [  5.85646521e+06,   1.73486370e+10,   0.00000000e+00],
            [  1.00000000e+01,   0.00000000e+00,  -1.00000000e+01]]
       
      q3 = [[  5.47068987e+03,   5.63388151e+01,   2.00000000e+03],
            [  4.50710521e+03,   1.56014336e+06,   0.00000000e+00],
            [  1.00000000e+01,   0.00000000e+00,  -1.00000000e+01]]
      """
     When H is computed from the scalar or sequence <input>
     Then H returns the matrix or sequence of matrices <output>
   
    Examples:

      |      model        | input                                | output           |
      | transpose classic | 0                                    | q0               |    
      | transpose classic | 0, -1000, -2e5                       | [q0, q1, q2]     |
      | transpose classic | (0, -1000, -2e5)                     | [q0, q1, q2]     |
      | transpose classic | [0, -2e5, -1000, -1e5]               | [q0, q2, q1, q3] |
      | transpose classic | array([0, -2e5, -1000, -1e5])        | [q0, q2, q1, q3] |
      | transpose classic | array([0, -1000, -2e5, -1e5])[::-2]  | [q3, q1]         |
