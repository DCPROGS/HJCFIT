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

#ifndef DCPROGS_EQUILIBRIUM_STATES_H
#define DCPROGS_EQUILIBRIUM_STATES_H

#include <DCProgsConfig.h>

namespace DCProgs {

  class IdealG;
  class MissedEventsG;

  //! \brief Solves the occupancy equation for initial and final states
  //! \details The equilibrium equation is \f\[\phi = \phi M\f\], \f\[\sum_i \phi_i = 1\f\], where
  //!          \f$M$\f is for initial states \f$\mathcal{G}_{AF}(s=0) \mathcal{G}_{FA}(s=0)\f$. The
  //!          problem is solved using Eigen's linear least-square utility, adding an extra row to
  //!          the matrix to impose the second condition.
  t_initvec MSWINDOBE occupancies(IdealG const &, bool _initial = true);

  //! \brief Solves the occupancy equation for initial and final states
  //! \details The equilibrium equation is \f\[\phi = \phi M\f\], \f\[\sum_i \phi_i = 1\f\], where
  //!          \f$M$\f is for initial states \f${}^e\mathcal{G}_{AF}(s=0) {}^e\mathcal{G}_{FA}(s=0)\f$. The
  //!          problem is solved using Eigen's linear least-square utility, adding an extra row to
  //!          the matrix to impose the second condition.
  t_initvec MSWINDOBE occupancies(MissedEventsG const &, bool _initial = true);
  
  //! \brief Solves the CHS occupancy equation for initial and final states
  //! \details For `_initial=true` this is Eq 5.11 of \cite colquhoun:1996, whereas for
  //! `_initial=false`, it is Eq 5.8.
  //! \note The initial and final vectors are not quite as satifyingly symmetric as for other
  //! occupancies. However, to simplify the API, this function still deals with both cases.
  t_initvec MSWINDOBE CHS_occupancies(MissedEventsG const &, t_real, bool _initial = true);
}
#endif 
