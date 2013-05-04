#ifndef DCPROGS_LIKELIHOOD_IDEALG
#define DCPROGS_LIKELIHOOD_IDEALG

#include <DCProgsConfig.h>

#include <unsupported/Eigen/MatrixFunctions>

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
      auto aa(t_real t) const -> decltype( t_rmatrix::Zero(3, 3) ) { return std::move(t_rmatrix::Zero(3, 3)); }
      //! Shut to open transitions.
      auto fa(t_real t) const -> decltype( (t*StateMatrix::ff()).exp()*StateMatrix::fa() )
        { return (t*StateMatrix::ff()).exp()*StateMatrix::fa(); }
      //! Open to shut transitions.
      auto af(t_real t) const -> decltype( (t*StateMatrix::aa()).exp()*StateMatrix::af() )
        { return (t*StateMatrix::aa()).exp()*StateMatrix::af(); }
      //! Shut to shut transitions.
      auto ff(t_real t) const -> decltype( t_rmatrix::Zero(3, 3) ) { return t_rmatrix::Zero(3, 3); }
  };
}

#endif
