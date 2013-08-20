Feature: QMatrix bindings

  The following tests the c++ to python bindings of QMatrix. Most of the tests are about making
  sure that c++ objects are translated correctly, or signal somehow that translation could not
  happen. Sometimes the bindings are shallow wrappers around c++ objects. Othertimes they may copies
  of c++ objects. In any case, values modified in python should be reflected in C and vice-versa.
  Below, a "numpy array" is a python object which represents a matrix. 

  Scenario: Initialize empty matrix
    Given QMatrix is accessible
    When  we instantiate QMatrix without arguments
    Then  nopen is 0
    And   matrix is a numpy array
    And   matrix is empty



  Scenario: Initialize empty matrix with arguments
    Given QMatrix is accessible
    When  we instantiate QMatrix with empty and 0
    Then  nopen is 0
    And   matrix is a numpy array
    And   matrix is empty




  Scenario Outline: Initialize with non-empty arguments
    Given QMatrix is accessible
    When  we instantiate QMatrix with <matrix> and <nopen>
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
      | classic  | 0     |
      | classic  | 1     |
      | classic  | 2     |
      | classic  | 3     |
      | classic  | 4     |
      | classic  | 5     |



  Scenario Outline: Initialization throws
    Given QMatrix is accessible
    When  we instantiate QMatrix with <matrix> and <nopen>
    Then  instantiation threw <exception> 
 
    Examples:
 
      | matrix     | nopen | exception                |
      | classic    | -1    | ValueError               |
      | classic    | 6     | ValueError               |
      | spam       | 0     | (ValueError, TypeError)  |
      | numpy_spam | 0     | TypeError                |


  Scenario Outline: Modify matrix in-place

    Given a QMatrix instantiated with classic and 0
    When  item <item> is set to <value>
    Then  item <item> is <value>

    Examples:
      | item   | value |
      | (0, 0) |     0 |
      | (0, 0) |     1 |
      | (1, 3) |  3.14 |

  Scenario: Modify matrix in-place

    Given a QMatrix instantiated with classic and 0
    When  nopen is set to 2
    Then  nopen is 2
