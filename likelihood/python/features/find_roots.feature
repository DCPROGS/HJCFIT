Feature: Figure out roots for a state matrix

  Given a transition matrix and the number of open states, we want to figure out its roots.
  There are several aspectes to this:
  
    - figuring out the bounds where all roots resides
    - figuring out the interval where each root resides
    - dealing with the possibility of complex eigenvalues in the transition matrix.


  Scenario: Find intervals of classic matrix
    Given the open-states determinantal equation "classic" with tau=1e-4
    When the root intervals are computed 
    Then no exception was thrown
     And there are 2 intervals 
     And interval 0 contains -3045.285776
     And interval 1 contains -162.929465
