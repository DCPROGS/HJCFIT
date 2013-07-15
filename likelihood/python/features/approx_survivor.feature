Feature: Check approximate survivor implementation.

  The following tests some results for the approximate survivor function.

  Scenario Outline: At t=0, the survivor function should be one
    Given a list of 20 random approximate survivor functions with tau=1e-4
      And a parameter tolerance=1e-8
     When the <name> block of the survivor function is computed for t=0
     Then the value is the identity

    Examples:

      | name |
      | af   |
      | fa   |

