#ifndef DCPROGS_LIKELIHOOD_IDEALG
#define DCPROGS_LIKELIHOOD_IDEALG

#include <DCProgsConfig.h>
#include "errors.h"

//! General namespace for all things DCProgs.
namespace DCProgs {

  //! \brief Ideal transition matrix of open and shut intervals
  //! \details Given a transition matrix $Q$ it is possible to figure out the evolution of any given
  //! system. 
  class IdealG {

    public:
      //! Constructor
      IdealG() : Q_(0, 0), nopen_(0) {}
      //! Destructor 
      virtual ~IdealG() {}; 
  
  
      //! \brief Sets Q matrix and the number of open states.
      //! \details Enforces \f[Q_{ii} = -\sum_{j\neqi} Q_{ij}]\f.
      //!          It is expected that open states are in the top-most corner.
      void set(t_rmatrix const &Q, t_int const &_nopen);
      //! Gets Q matrix. 
      t_rmatrix const & get_Q() const { return Q_; }
      //! Gets the number of open states
      t_int const & get_nopen() const { return nopen_; }

      //! Open to open transitions.
      t_rmatrix aa() const { return Q_.topLeftCorner(nopen_, nopen_); }
      //! Open to shut transitions.
      t_rmatrix af() const { return Q_.topRightCorner(nopen_, nopen_); }
      //! Shut to open transitions.
      t_rmatrix fa() const { return Q_.bottomLeftCorner(nopen_, nopen_); }
      //! Shut to shut transitions.
      t_rmatrix ff() const { return Q_.bottomRightCorner(nopen_, nopen_); }
  
    protected:
      //! Transition matrix Q
      t_rmatrix Q_;
      //! Number of open states
      t_int nopen_;
  };
}

#endif
