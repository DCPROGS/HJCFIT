#include <iostream>

#include <likelihood/determinant_equation.h>
#include <likelihood/root_finder.h>
#include <likelihood/brentq.h>
 
int main() {

  // Define parameters.
  DCProgs::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;
  DCProgs::QMatrix qmatrix(matrix, /*nopen=*/2);
  DCProgs::DeterminantEq det(qmatrix, 1e-4);

  // Find upper and lower bound
  DCProgs::t_real upper_bound = DCProgs::find_upper_bound_for_roots(det);
  DCProgs::t_real lower_bound = DCProgs::find_lower_bound_for_roots(det);


  // computes eigenvalues of H(s) for given s
  auto get_eigenvalues = [&det](DCProgs::t_real _s) -> DCProgs::t_rvector {
    return Eigen::EigenSolver<DCProgs::t_rmatrix>(det.H(_s)).eigenvalues().real();
  };

  // Checks bounds are correct.
  if((get_eigenvalues(lower_bound).array() < lower_bound).any()) 
    throw DCProgs::errors::Runtime("Incorrect lower bound.");
  if((get_eigenvalues(upper_bound).array() > upper_bound).any()) 
    throw DCProgs::errors::Runtime("Incorrect upper bound.");

  std::cout << "Root Determination\n"
               "==================\n\n" 
            << "  * Interval containing roots: " << upper_bound << ", " << lower_bound << "\n"
            << "  * Eigenvalues of H at lower bound: " 
            << get_eigenvalues(lower_bound).transpose() << "\n"
            << "  * Eigenvalues of H at upper bound: " 
            << get_eigenvalues(upper_bound).transpose() << "\n\n";

  // Figure out bracket for each root.
  std::vector<DCProgs::RootInterval> intervals
    = DCProgs::find_root_intervals(det, lower_bound, upper_bound);

  // Find root for each interval
  for(DCProgs::RootInterval const& interval: intervals) {
    auto brentq_result = DCProgs::brentq(det, interval.start, interval.end);
    std::cout << "  * Root interval: [" << interval.start << ", " << interval.end << "]\n" 
              << "    Corresponding root: " << std::get<0>(brentq_result) << "\n\n";
  }

  // Look for roots in one go.
  std::vector<DCProgs::Root> roots = DCProgs::find_roots(det);
  std::cout <<  "  * All roots: ";
  for(DCProgs::Root const &root: roots) std::cout << root.root << " ";
  std::cout << "\n";

  return 0;
}
