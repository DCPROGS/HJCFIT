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
                   << "========================\n\n" 
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

  MissedEventsG::MissedEventsG( QMatrix const &_qmatrix, t_real _tau,
                                t_uint _nmax, t_real _xtol, t_real _rtol, t_uint _itermax,
                                t_real _lowerbound, t_real _upperbound )
# ifdef HAS_CXX11_CONSTRUCTOR_DELEGATE
    : MissedEventsG( _qmatrix, _tau,
                     [_xtol, _rtol, _itermax, _lowerbound, _upperbound](DeterminantEq const &_c) {
                        return find_roots(_c, _xtol, _rtol, _itermax, _lowerbound, _upperbound); 
                     }, _nmax ) {}
# else
    : ExactSurvivor(_qmatrix, _tau),
      ApproxSurvivor(_qmatrix, _tau,
          [_xtol, _rtol, _itermax, _lowerbound, _upperbound](DeterminantEq const &_c) {
            return find_roots(_c, _xtol, _rtol, _itermax, _lowerbound, _upperbound); 
          }),
      laplace_a_(new LaplaceSurvivor(_qmatrix)),
      laplace_f_(new LaplaceSurvivor(_qmatrix.transpose())),
      nmax_(_nmax), tmax_(_tau*t_real(_nmax-1)),
      af_factor_(_qmatrix.af() * (_tau * _qmatrix.ff()).exp()),
      fa_factor_(_qmatrix.fa() * (_tau * _qmatrix.aa()).exp()) {}
# endif
}
