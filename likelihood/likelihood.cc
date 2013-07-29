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
#include "occupancies.h"
#include "missed_eventsG.h"

namespace DCProgs {

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, t_Bursts const & _self) {
    _stream << "Bursts:\n"
            << "-------\n"
            << "  [ ";
    for(t_int i(0); i < _self.size() - 1; ++i) _stream << _self[i] << ",\n    ";
    _stream << _self.back() << " ]\n";
    return _stream;
  }

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, t_Burst const & _self) {
    Eigen::Map<const t_initvec> vector(&(_self[0]), _self.size());
    return _stream << vector.format(Eigen::IOFormat(Eigen::FullPrecision, 0, ",", ",\n",
                                                    "", "", "[", "]" )); 
  }

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, Log10Likelihood const & _self) {
    
    _stream << "Log10 Likelihood:\n"
            << "=================n\n" 
            << "  * Number of open states: " << _self.nopen << "\n"
            << "  * Resolution time tau: " << _self.tau << "\n";
    if(_self.tcritical <= 0e0) _stream << "  * Using equilibrium occupancies.\n";
    else _stream << "  * Using CHS occupancies with tcrit: "  << _self.tcritical << "\n";
    _stream << "  * Exact events computed for: t < "
            << _self.nmax << " tau\n\n"
            << _self.bursts
            << "\nRoot Finding:\n"
            << "-------------\n"
            << "  * Tolerance criteria: " << _self.xtol << ", " << _self.rtol << "\n"
            << "  * Maximum number of iterations: " << _self.itermax << "\n";
    return _stream;
  }

  t_real Log10Likelihood::operator()(QMatrix const &_matrix) const {
    MissedEventsG const eG = create_missed_eventsG(_matrix, tau, nmax, xtol, rtol, itermax);
    t_rvector const final = tcritical > 0 ?
                              CHS_occupancies(eG, tcritical, false).transpose():
                              occupancies(eG, false).transpose();
    t_initvec const initial = tcritical > 0 ? CHS_occupancies(eG, tcritical): occupancies(eG);
                                
    t_real result(0);
    for(t_Burst const &burst: bursts) 
      result += chained_log10_likelihood(eG, burst.begin(), burst.end(), initial, final);
    return result;
  }
  t_rvector Log10Likelihood::vector(QMatrix const &_matrix) const {
    MissedEventsG const eG = create_missed_eventsG(_matrix, tau, nmax, xtol, rtol, itermax);
    t_rvector const final = tcritical > 0 ?
                              CHS_occupancies(eG, tcritical, false).transpose():
                              occupancies(eG, false).transpose();
    t_initvec const initial = tcritical > 0 ? CHS_occupancies(eG, tcritical): occupancies(eG);
                                
    t_rvector result(bursts.size());
    t_int i(0);
    for(t_Burst const &burst: bursts) {
      result(i++) = chained_log10_likelihood(eG, burst.begin(), burst.end(), initial, final);
    }
    return result;
  }
}
