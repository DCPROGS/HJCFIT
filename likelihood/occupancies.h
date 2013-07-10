#ifndef DCPROGS_EQUILIBRIUM_STATES_H
#define DCPROGS_EQUILIBRIUM_STATES_H

#include <DCProgsConfig.h>

namespace DCProgs {

  class IdealG;

  //! \brief Solves the occupancy equation for initial and final states
  //! \details The equilibrium equation is \f\[\phi = \phi M\f\], \f\[\sum_i \phi_i = 1\f\], where
  //!          \f$M$\f is for initial states \f$\mathcal{G}_{AF}(s=0) \mathcal{G}_{FA}(s=0)\f$. The
  //!          problem is solved using Eigen's linear least-square utility, adding an extra row to
  //!          the matrix to impose the second condition.
  t_initvec MSWINDOBE occupancies(IdealG const &, bool _initial = true);
}
#endif 
