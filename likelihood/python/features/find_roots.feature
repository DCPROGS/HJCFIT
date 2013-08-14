Feature: Figure out roots for a state matrix

  Given a transition matrix and the number of open states, we want to figure out its roots.
  There are several aspectes to this:
  
    - figuring out the bounds where all roots resides
    - figuring out the interval where each root resides
    - An irreducible(?) reversible Markove process has n _real_ roots, where n is the number of open
      states.

  At present, I'm not sure how to create random matrices that strictly respect point 3 above. Hence
  the failure allowance in some of the scenario below.

  Most of the tests below are designed to ignore failures do to issues such as convergence problems
  when computing eigenvalues. This is to make it easier to detect failure due to (not) finding
  roots.

  Scenario: Find intervals of classic open-state matrix
    Given the open-states determinant equation "classic" with tau=1e-4
     When the root intervals are computed 
     Then no exception was thrown
      And there are 2 intervals 
      And interval 0 contains -3045.285776037674
      And interval 1 contains -162.92946543451328
    
  Scenario: Find intervals of classic closed-state matrix
    Given the closed-states determinant equation "classic" with tau=1e-4
     When the root intervals are computed 
     Then no exception was thrown
      And there are 3 intervals 
      And interval 0 contains -17090.192769236815
      And interval 1 contains -2058.0812921673496
      And interval 2 contains -0.24356535498785126

  Scenario Outline: Compare brute force and sensible search for intervals
    Given the <type>-states determinant equation "<equation>" with tau=1e-4
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
    Given the <type>-states determinant equation "<equation>" with tau=1e-4
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


  Scenario: Check roots add up for many random cases
    Given a list of 500 random determinant equations
      And allowing for 2% failure in tests below
     When the roots are computed for each
     Then roots are roots indeed, to within tolerance=1e-5 * variation of det W
      And the multiplicity of the real roots add up to the number of open states
