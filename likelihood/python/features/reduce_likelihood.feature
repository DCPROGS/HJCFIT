Feature: Reduce number of components to optimize in likelihood

  The likelihood functions provided by dcprogs take the full Q-matrix. However, a fair number of the
  components should be fixed, or depend on others. In other words, we need optimize only on a subset
  of components. The function `dcprogs.likelihood.optimization.reduce_likelihood` creates  a
  function that take the minimum number of components only.


  Scenario: Check that fixed components are fixed and variable components variable
    Given the graph matrix below
          """
          [["V", "V", "V",   0,  0.1],
           ["V", "V",   0, "V",    0],
           ["V",   0, "V", "V",  "V"],
           [  0, "V", "V", "V",    0],
           [  1.1,   0, "V",   0,  "V"]] 
          """
     When the reduced likelihood is created 
      And the reduced likelihood called with a random vector 
     Then the sum of rows of the resulting matrix is zero
      And the fixed components of the resulting matrix are those above
      And the variable components of the resulting matrix are from the random matrix

  Scenario: Check that expression components work
    Given the graph matrix below
          """
          [[  "V", "V",                       "V",           0,  0.1],
           [  "V", "V",                         0, "2*q[0, 2]",    0],
           [  "V",   0,                       "V",         "V",  "V"],
           [    0, "V",                       "V",         "V",    0],
           [  1.1,   0, "q[0, 2] + 3*cos(q[2,4])",           0,  "V"]] 
          """
     When the reduced likelihood is created 
      And the reduced likelihood called with a random vector 
     Then the sum of rows of the resulting matrix is zero
      And the fixed components of the resulting matrix are those above
      And the variable components of the resulting matrix are from the random matrix
      And the two expression components come out as expected
