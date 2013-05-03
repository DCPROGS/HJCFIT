#ifndef DCPROGS_LIKELIHOOD_IDEALG
#define DCPROGS_LIKELIHOOD_IDEALG

#include <DCProgsConfig.h>

#include <type_traits>

#include "state_matrix.h"
#include "errors.h"

//! General namespace for all things DCProgs.
namespace DCProgs {

  //! \brief Ideal transition matrix of open and shut intervals
  //! \details Given a transition matrix $Q$ it is possible to figure out the evolution of any given
  //! system. 
  class IdealG : protected StateMatrix {

    public:
      //! Constructor
      IdealG() : StateMatrix() {}
      //! Destructor 
      virtual ~IdealG() {}; 
  
  
      //! \brief Sets Q matrix and the number of open states.
      //! \details Enforces \f[Q{ii} = -\sum_{j\neqi} Q{ij}]\f.
      //!          It is expected that open states are in the top-most corner.
      void set(t_rmatrix const &Q, t_int const &_nopen);
      //! Gets Q matrix. 
      t_rmatrix const & get_Q() const { return this->matrix; }
      //! Gets the number of open states
      t_int const & get_nopen() const { return this->nopen; }

      //! Open to open transitions.
      t_rmatrix af(t_real _t0) const -> decltype( matrix.topLeftCorner(nopen, nopen) ) 
//       { return matrix.topLeftCorner(nopen, nopen); }
//     //! Open to shut transitions.
//     t_rmatrix af(t_real) const;
//     //! Shut to open transitions.
//     t_rmatrix fa(t_real) const;
//     //! Shut to shut transitions.
//     t_rmatrix ff(t_real) const;
  };
}

#endif
