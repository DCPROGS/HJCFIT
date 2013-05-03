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
      IdealG() : Q_(0, 0)  { open_.resize(0);} 
      //! Destructor 
      virtual ~IdealG() {}; 
  
  
      //! \brief Sets Q matrix and its partition
      //! \details Enforces \f[Q_{ii} = -\sum_{j\neqi} Q_{ij}]\f.
      void set(t_rmatrix const &Q, t_partition const &_open);
      //! Gets Q matrix. 
      t_rmatrix const & get_Q() const { return Q_; }
      //! Sets open_ matrix. 
      t_partition const & get_open_states() const { return open_; }
  
    protected:
      //! Transition matrix Q
      t_rmatrix Q_;
      //! Partition vector indicating open states.
      t_partition open_;
  };
}

#endif
