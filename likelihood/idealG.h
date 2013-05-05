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
      //! \brief Constructor with parameters.
      //! \details Calls set method with input parameters.
      //! \param[in] _matrix: Any matrix or matrix expression from Eigen. Will become the transition
      //!                     matrix. Diagonal elements are transformed as explain in set(). Open
      //!                     states should be in the top rows.
      //! \param[in] _nopen: Number of open states. 
      //! \throws errors::Domain if input has incorrect values or size.
      template<class T>
        IdealG(Eigen::DenseBase<T> const &_matrix, t_int _nopen);
      //! Destructor 
      virtual ~IdealG() {}; 
  
      //! \brief Sets Q matrix and the number of open states.
      //! \details Enforces \f[Q{ii} = -\sum_{j\neqi} Q{ij}]\f.
      //!          It is expected that open states are the top rows [0, _nopen].
      void set(t_rmatrix const &_Q, t_int const &_nopen);
      //! Sets state matrix on which to act.
      void set(StateMatrix const &_in) { set(_in.matrix, _in.nopen); }
      //! Gets Q matrix. 
      t_rmatrix const & get_Q() const { return this->matrix; }
      //! Gets the number of open states
      t_int const & get_nopen() const { return this->nopen; }

      //! Open to open transitions.
      auto aa(t_real t) const -> decltype( t_rmatrix::Zero(nopen, nopen) )
        { return std::move(t_rmatrix::Zero(nopen, nopen)); }
      //! Shut to shut transitions.
      t_rmatrix ff(t_real t) const;
      //! Shut to open transitions.
      auto fa(t_real t) const -> decltype( (t*StateMatrix::ff()).exp()*StateMatrix::fa() )
        { return (t*StateMatrix::ff()).exp()*StateMatrix::fa(); }
      //! Open to shut transitions.
      auto af(t_real t) const -> decltype( (t*StateMatrix::aa()).exp()*StateMatrix::af() )
        { return (t*StateMatrix::aa()).exp()*StateMatrix::af(); }

      //! Laplace transform of open to open transitions.
      auto laplace_aa(t_real s) const -> decltype(aa(0)) { return this->aa(0); }
      //! Laplace transform of shut to shut transitions.
      auto laplace_ff(t_real s) const -> decltype(ff(0)) { return this->ff(0); }
      //! Laplace transform of shut to open transitions.
      t_rmatrix laplace_fa(t_real s) const;
      //! Open to shut transitions.
      t_rmatrix laplace_af(t_real t) const;
  };

  template<class T>
    IdealG :: IdealG   (Eigen::DenseBase<T> const &_matrix, t_int _nopen)
                     : StateMatrix() {
      try { this->set(_matrix, _nopen); }
      catch(...) {
        this->matrix.resize(0, 0);
        this->nopen = 0;
        throw;
      }
    }
}

#endif
