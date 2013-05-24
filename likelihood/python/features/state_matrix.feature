Feature: StateMatrix bindings

  The following tests the c++ to python bindings of StateMatrix. Most of the tests are about making
  sure that c++ objects are translated correctly, or signal somehow that translation could not
  happen. Sometimes the bindings are shallow wrappers around c++ objects. Othertimes they may copies
  of c++ objects. In any case, values modified in python should be reflected in C and vice-versa.

  Scenario: Initialize empty matrix
    Given StateMatrix is accessible
    When  we instantiate StateMatrix without arguments
    Then  nopen is 0
    And   matrix is a numpy array
    And   matrix is empty



  Scenario: Initialize empty matrix with arguments
    Given StateMatrix is accessible
    When  we instantiate StateMatrix with empty and 0
    Then  nopen is 0
    And   matrix is a numpy array
    And   matrix is empty




  Scenario Outline: Initialize with non-empty arguments
    Given StateMatrix is accessible
    When  we instantiate StateMatrix with <matrix> and <nopen>
    Then  instantiation did not throw
    And   nopen is <nopen>
    And   matrix is a numpy array
    And   matrix is <matrix>
    And   aa is given by rows=":nopen" and cols=":nopen" of <matrix>, and nopen=<nopen>
    And   af is given by rows=":nopen" and cols="nopen:" of <matrix>, and nopen=<nopen>
    And   fa is given by rows="nopen:" and cols=":nopen" of <matrix>, and nopen=<nopen>
    And   ff is given by rows="nopen:" and cols="nopen:" of <matrix>, and nopen=<nopen>
 
    Examples: 
 
      | matrix   | nopen |
      | Qmatrix  | 0     |
      | Qmatrix  | 1     |
      | Qmatrix  | 2     |
      | Qmatrix  | 3     |
      | Qmatrix  | 4     |
      | Qmatrix  | 5     |



  Scenario Outline: Initialization throws
    Given StateMatrix is accessible
    When  we instantiate StateMatrix with <matrix> and <nopen>
    Then  instantiation threw <exception> 
 
    Examples:
 
      | matrix     | nopen | exception  |
      | Qmatrix    | -1    | ValueError |
      | Qmatrix    | 6     | ValueError |
      | spam       | 0     | TypeError  |
      | numpy_spam | 0     | TypeError  |


  Scenario Outline: Modify matrix in-place

    Given a StateMatrix instantiated with Qmatrix and 0
    When  item <item> is set to <value>
    Then  item <item> is <value>

    Examples:
      | item   | value |
      | (0, 0) |     0 |
      | (0, 0) |     1 |
      | (1, 3) |  3.14 |

  Scenario: Modify matrix in-place

    Given a StateMatrix instantiated with Qmatrix and 0
    When  nopen is set to 2
    Then  nopen is 2
