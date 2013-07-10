#include <DCProgsConfig.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "errors.h"
#include "approx_survivor.h"


namespace DCProgs {

  ApproxSurvivor :: ApproxSurvivor(DeterminantEq const &_af, std::vector<Root> const &_roots_af, 
                                   DeterminantEq const &_fa, std::vector<Root> const &_roots_fa ) {

    asymptotes_af_.reset(new Asymptotes(_af, _roots_af));
    if(not asymptotes_af_.get()) throw errors::Runtime("Could not initialize unique_ptr");
    asymptotes_fa_.reset(new Asymptotes(_fa, _roots_fa));
    if(not asymptotes_fa_.get()) throw errors::Runtime("Could not initialize unique_ptr");
  }
 
  // Function to create approximate missed event survivor function.
  ApproxSurvivor::ApproxSurvivor(QMatrix const &_qmatrix, t_real _tau, t_RootFinder const &_findroots) {
    // First creates determinantal equations.
    DeterminantEq determinant_af(_qmatrix, _tau, true);
    DeterminantEq determinant_fa(_qmatrix, _tau, false);
    // Then finds roots
    std::vector<Root> roots_af = _findroots(determinant_af);
    std::vector<Root> roots_fa = _findroots(determinant_fa);
    // Then creates Asymptotes object
    asymptotes_af_.reset(new Asymptotes(determinant_af, roots_af));
    if(not asymptotes_af_.get()) throw errors::Runtime("Could not initialize unique_ptr");
    asymptotes_fa_.reset(new Asymptotes(determinant_fa, roots_fa));
    if(not asymptotes_fa_.get()) throw errors::Runtime("Could not initialize unique_ptr");
  }
}
