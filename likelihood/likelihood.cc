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
    for(t_Bursts::size_type i(0); i + 1 < _self.size(); ++i) _stream << _self[i] << ",\n    ";
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
            << "=================\n\n" 
            << "  * Number of open states: " << _self.nopen << "\n"
            << "  * Resolution time tau: " << _self.tau << "\n";
    if(DCPROGS_ISNAN(_self.tcritical)  or _self.tcritical <= 0e0)
         _stream << "  * Using equilibrium occupancies.\n";
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
    verify_qmatrix(_matrix);
    MissedEventsG const eG = MissedEventsG( _matrix, tau, nmax, xtol, rtol, itermax,
                                            lower_bound, upper_bound );

    bool const eq_vector = DCPROGS_ISNAN(tcritical) or tcritical <= 0;

    t_rvector final;

    if(eq_vector)
        final = t_rvector::Ones(_matrix.nshut(),1);
    else
        final = CHS_occupancies(eG, tcritical, false).transpose();

    t_initvec const initial = eq_vector ? occupancies(eG): CHS_occupancies(eG, tcritical);
                                
    t_real result(0);
    for(t_Burst const &burst: bursts) 
      result += chained_log10_likelihood(eG, burst, initial, final, omp_num_threads);
    return result;
  }
  t_rvector Log10Likelihood::vector(QMatrix const &_matrix) const {
    verify_qmatrix(_matrix);
    MissedEventsG const eG = MissedEventsG( _matrix, tau, nmax, xtol, rtol, itermax,
                                            lower_bound, upper_bound );
    bool const eq_vector = DCPROGS_ISNAN(tcritical) or tcritical <= 0;

    t_rvector final;

    if(eq_vector)
        final = t_rmatrix::Ones(_matrix.nshut(),1);
    else
        final = CHS_occupancies(eG, tcritical, false).transpose();

    t_initvec const initial = eq_vector ? occupancies(eG): CHS_occupancies(eG, tcritical);
                                
    t_rvector result(bursts.size());
    t_int i(0);
    for(t_Burst const &burst: bursts) {
      result(i++) = chained_log10_likelihood(eG, burst, initial, final, omp_num_threads);
    }
    return result;
  }
}
