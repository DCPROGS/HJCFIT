#include <iostream>

#include <likelihood/determinant_equation.h>
#include <likelihood/root_finder.h>
#include <likelihood/brentq.h>
 
int main() {

  // Define parameters.
  HJCFIT::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;
  HJCFIT::QMatrix qmatrix(matrix, /*nopen=*/2);
  HJCFIT::DeterminantEq det(qmatrix, 1e-4);

  // Find upper and lower bound
  HJCFIT::t_real upper_bound = HJCFIT::find_upper_bound_for_roots(det);
  HJCFIT::t_real lower_bound = HJCFIT::find_lower_bound_for_roots(det);


  // computes eigenvalues of H(s) for given s
  auto get_eigenvalues = [&det](HJCFIT::t_real _s) -> HJCFIT::t_rvector {
    return Eigen::EigenSolver<HJCFIT::t_rmatrix>(det.H(_s)).eigenvalues().real();
  };

  // Checks bounds are correct.
  if((get_eigenvalues(lower_bound).array() < lower_bound).any()) 
    throw HJCFIT::errors::Runtime("Incorrect lower bound.");
  if((get_eigenvalues(upper_bound).array() > upper_bound).any()) 
    throw HJCFIT::errors::Runtime("Incorrect upper bound.");

  std::cout << "Root Determination\n"
               "==================\n\n" 
            << "  * Interval containing roots: " << upper_bound << ", " << lower_bound << "\n"
            << "  * Eigenvalues of H at lower bound: " 
            << get_eigenvalues(lower_bound).transpose() << "\n"
            << "  * Eigenvalues of H at upper bound: " 
            << get_eigenvalues(upper_bound).transpose() << "\n\n";

  // Figure out bracket for each root.
  std::vector<HJCFIT::RootInterval> intervals
    = HJCFIT::find_root_intervals(det, lower_bound, upper_bound);

  // Find root for each interval
  for(HJCFIT::RootInterval const& interval: intervals) {
    auto brentq_result = HJCFIT::brentq(det, interval.start, interval.end);
    std::cout << "  * Root interval: [" << interval.start << ", " << interval.end << "]\n" 
              << "    Corresponding root: " << std::get<0>(brentq_result) << "\n\n";
  }

  // Look for roots in one go.
  std::vector<HJCFIT::Root> roots = HJCFIT::find_roots(det);
  std::cout <<  "  * All roots: ";
  for(HJCFIT::Root const &root: roots) std::cout << root.root << " ";
  std::cout << "\n";

  return 0;
}
