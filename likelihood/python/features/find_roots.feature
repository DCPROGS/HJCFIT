Feature: Figure out roots for a state matrix

  Given a transition matrix and the number of open states, we want to figure out its roots.

  Scenario: Find intervals of classic matrix
    Given the open-states determinantal equation "classic" with tau=1e-4
    When the root intervals are computed 
    Then no exception was thrown
     And there are 2 intervals 
     And interval 0 contains -3045.285776
     And interval 1 contains -162.929465
