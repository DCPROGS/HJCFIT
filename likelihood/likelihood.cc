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

#include <sstream>

#include "likelihood.h"

namespace DCProgs {

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, Log10Likelihood const & _self) {
    
    _stream << "Log10 Likelihood:\n"
            << "=====================\n\n" 
            << "  * Number of open states:\n" << _self.nopen << "\n";
            << "  * Resolution time tau: " << _self.tau << "\n";
    if(_self.tcrit > 0e0) 
      _stream << "  * Using equilibrium occupancies.\n"
    else _stream << "  * Using CHS occupancies with tcri: "  << _self.tcrit << "\n";
    _stream << "  * Number of intervals for which to compute exact likelihood: "
            << _self.nmax << "\n"
            << "Root Finding:\n"
            << "-------------\n"
            << "  * Tolerance criteria: " << _self.xtol << ", " << _self.rtol << "\n"
            << "  * Maximum number of iterations: " << _self.itermax << "\n";
    _stream << "Bursts:\n"
            << "-------\n  [ ";

    auto lambda = [&_stream]( std::vector<t_real> const & _burst) {
      t_rvector burst(_burst.size());
      for(t_rvector::Index i(0); i < burst.size(); ++i) burst[i] = _burst[i];
      return numpy_io(t_rvector);
    };
    for(t_int i(0); i < _self.bursts.size() - 1; ++i)
      _stream << lambda(_self.bursts[i]) << ",\n    ";
    _stream << lambda(_self.bursts.back()) << " ]"
    return _stream;
  }

  t_real Log10Likelihood::operator()(QMatrix const &_matrix) {
    MissedEventsG const eG = create_missed_events(_matrix, tau, nmax, xtol, rtol, itermax);
    t_rvector const final = self.tcritical > 0 ?
                              occupancies(eG, false).transpose():
                              CHS_occupancies(eG, tcritical, false).transpose();
    t_initvec const initial = self.tcritical > 0 ? occupancies(eG): CHS_occupancies(eG, tcritical);
                                
    t_real result(0);
    for(std::vector<t_real> const &burst: bursts) 
      result += chained_log10_likelihood(eG, burst.begin(), burst.end(), initial, final);
    return result;
  }
  t_rvector Log10Likelihood::vector(QMatrix const &_matrix) {
    MissedEventsG const eG = create_missed_events(_matrix, tau, nmax, xtol, rtol, itermax);
    t_rvector const final = self.tcritical > 0 ?
                              occupancies(eG, false).transpose():
                              CHS_occupancies(eG, tcritical, false).transpose();
    t_initvec const initial = self.tcritical > 0 ? occupancies(eG): CHS_occupancies(eG, tcritical);
                                
    t_rvector result(bursts.size());
    t_int i(0);
    for(std::vector<t_real> const &burst: bursts) {
      result(i) = chained_log10_likelihood(eG, burst.begin(), burst.end(), initial, final);
    }
    return result;
  }
}
