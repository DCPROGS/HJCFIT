Feature: Figure out roots for a state matrix

  Given a transition matrix and the number of open states, we want to figure out its roots.
  There are several aspectes to this:
  
    - figuring out the bounds where all roots resides
    - figuring out the interval where each root resides

  Scenario: Find intervals of classic open-state matrix
    Given the open-states determinantal equation "classic" with tau=1e-4
     When the root intervals are computed 
     Then no exception was thrown
      And there are 2 intervals 
      And interval 0 contains -3045.285776
      And interval 1 contains -162.929465
    
  Scenario: Find intervals of classic closed-state matrix
    Given the closed-states determinantal equation "classic" with tau=1e-4
     When the root intervals are computed 
     Then no exception was thrown
      And there are 3 intervals 
      And interval 0 contains -17090.1927692368
      And interval 1 contains -2058.08129216735
      And interval 2 contains -0.243565355

  Scenario Outline: Compare brute force and sensible search for intervals
    Given the <type>-states determinantal equation "<equation>" with tau=1e-4
     When the root intervals are computed
      And a brute force search for roots is perfomed with resolution=<resolution>
     Then the intervals larger than <resolution> do overlap
      And roots are found within each overlap for which sign changes

     Examples: 

      | type   | equation        | resolution | 
      | open   | classic         |  10        |
      | closed | classic         |  10        |
      | open   | singular matrix |  10        |
      | closed | singular matrix |  10        |
      | open   | random          |  10        |
      | open   | random          |  10        |
      | open   | random          |  10        |
      | open   | random          |  10        |


  Scenario Outline: Check roots in prospective intervals
    Given the <type>-states determinantal equation "<equation>" with tau=1e-4
     When the root intervals are computed
     Then roots are found within each interval for which sign changes

     Examples: 

      | type   | equation        | resolution | 
      | open   | classic         |  10        |
      | closed | classic         |  10        |
      | open   | singular matrix |  10        |
      | closed | singular matrix |  10        |
      | open   | random          |  10        |
      | open   | random          |  10        |
      | open   | random          |  10        |
      | open   | random          |  10        |
