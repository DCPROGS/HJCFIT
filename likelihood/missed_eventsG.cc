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

#include <iostream>
#include "missed_eventsG.h"
#include "root_finder.h"

namespace DCProgs {

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, MissedEventsG const & _self) {
    
    return _stream << "Missed-Event Likelihood:\n"
                   << "=======================\n\n" 
                   << "  * Transition Rate matrix:\n" << numpy_io(_self.get_qmatrix().matrix) << "\n"
                   << "  * Number of 'A' states: " << _self.get_qmatrix().nopen << "\n"
                   << "  * Resolution time tau: " << _self.get_tau() << "\n"
                   << "  * Exact events computed for: t < " << _self.get_nmax() << " tau\n";
  }
  
  // CHS matrices \f$H_{FA}\f$
  t_rmatrix CHS_matrix_Hfa(MissedEventsG const &_g, t_real _tcrit) {
    return partial_CHS_matrix(*_g.asymptotes_fa_, _g.get_tau(), _tcrit) * _g.get_fa_factor();
  }
  // CHS matrices \f$H_{FA}\f$
  t_rmatrix CHS_matrix_Haf(MissedEventsG const &_g, t_real _tcrit) {
    return partial_CHS_matrix(*_g.asymptotes_af_, _g.get_tau(), _tcrit) * _g.get_af_factor();
  }

  MissedEventsG create_missed_eventsG( QMatrix const &_matrix, t_real _tau, t_int _nmax,
                                       t_real _xtol, t_real _rtol, t_int _itermax) {

    DeterminantEq const determinant_af(_matrix, _tau);
    DeterminantEq const determinant_fa = determinant_af.transpose(); 
    std::vector<Root> const roots_af = find_roots(determinant_af, _xtol, _rtol, _itermax);
    std::vector<Root> const roots_fa = find_roots(determinant_fa, _xtol, _rtol, _itermax);
    return MissedEventsG(determinant_af, roots_af, determinant_fa, roots_fa);
  }
}
