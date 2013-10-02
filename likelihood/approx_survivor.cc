/***********************
    DCProgs computes missed-events likelihood as described in
    Hawkes, Jalali and Colquhoun (1990, 1992)

    Copyright (C) 2013  University College London

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
************************/

#include <DCProgsConfig.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "errors.h"
#include "approx_survivor.h"
#include "root_finder.h"


namespace DCProgs {

  ApproxSurvivor :: ApproxSurvivor(DeterminantEq const &_af, std::vector<Root> const &_roots_af, 
                                   DeterminantEq const &_fa, std::vector<Root> const &_roots_fa ) {

    asymptotes_af_.reset(new Asymptotes(_af, _roots_af));
    if(not asymptotes_af_.get()) throw errors::Runtime("Could not initialize unique_ptr");
    asymptotes_fa_.reset(new Asymptotes(_fa, _roots_fa));
    if(not asymptotes_fa_.get()) throw errors::Runtime("Could not initialize unique_ptr");
  }
 
# ifdef DCPROGS_MACRO
#   error DCPROGS_MACRO already defined
# endif
  // Macro is used to fake constructor delegation...
# define DCPROGS_MACRO(FINDROOTS)                                                          \
    /* First creates determinant equations. */                                             \
    DeterminantEq determinant_af(_qmatrix, _tau);                                          \
    DeterminantEq determinant_fa(determinant_af.transpose());                              \
    /* Then finds roots */                                                                 \
    std::vector<Root> roots_af = FINDROOTS(determinant_af);                                \
    std::vector<Root> roots_fa = FINDROOTS(determinant_fa);                                \
    /* Then creates Asymptotes object */                                                   \
    asymptotes_af_.reset(new Asymptotes(determinant_af, roots_af));                        \
    if(not asymptotes_af_.get()) throw errors::Runtime("Could not initialize unique_ptr"); \
    asymptotes_fa_.reset(new Asymptotes(determinant_fa, roots_fa));                        \
    if(not asymptotes_fa_.get()) throw errors::Runtime("Could not initialize unique_ptr"); 

  // Function to create approximate missed event survivor function.
  ApproxSurvivor::ApproxSurvivor(QMatrix const &_qmatrix, t_real _tau, t_RootFinder const &_findroots) {
    DCPROGS_MACRO(_findroots);
  }


  ApproxSurvivor::ApproxSurvivor( QMatrix const &_qmatrix, t_real _tau,
                                  t_real _xtol, t_real _rtol, t_uint _itermax,
                                  t_real _lowerbound, t_real _upperbound ) 
# ifdef HAS_CXX11_CONSTRUCTOR_DELEGATE
    : ApproxSurvivor( _qmatrix, _tau,
                      [_xtol, _rtol, _itermax, _lowerbound, _upperbound](DeterminantEq const &_c) {
                        return find_roots(_c, _xtol, _rtol, _itermax, _lowerbound, _upperbound); 
                      }) {}
# else
  {
    auto findroots = [_xtol, _rtol, _itermax, _lowerbound, _upperbound](DeterminantEq const &_c) {
      return find_roots(_c, _xtol, _rtol, _itermax, _lowerbound, _upperbound);  
    };
    DCPROGS_MACRO(findroots);
  }
# endif
  
  //! Dumps Survivor function equation to stream
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, ApproxSurvivor const &_self) {

    t_rvector af_roots(_self.nb_af_components()), fa_roots(_self.nb_fa_components());
    for(t_rvector::Index i(0); i < af_roots.size(); ++i) 
      af_roots(i) = std::get<1>(_self.get_af_components(i));
    for(t_rvector::Index i(0); i < fa_roots.size(); ++i) 
      fa_roots(i) = std::get<1>(_self.get_fa_components(i));

    return _stream << "Approximate Survivor function:\n"
                   << "==============================\n\n" 
                   << "  * af roots: " << af_roots.transpose() << "\n"
                   << "  * fa roots: " << fa_roots.transpose() << "\n";
  }
}
