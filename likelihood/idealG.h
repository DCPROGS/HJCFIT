#ifndef DCPROGS_LIKELIHOOD_IDEALG
#define DCPROGS_LIKELIHOOD_IDEALG

#include <DCProgsConfig.h>

#include <utility>

#include <unsupported/Eigen/MatrixFunctions>

#include "state_matrix.h"
#include "errors.h"

//! General namespace for all things DCProgs.
namespace DCProgs {

  //! \brief Ideal transition matrix of open and shut intervals
  //! \details Given a transition matrix $Q$ it is possible to figure out the evolution of any given
  //! system. 
  class MSWINDOBE IdealG : protected StateMatrix {

    //! Just trying to figure out a complex return type...
    typedef decltype( (t_real(0) * std::declval<const StateMatrix>().ff()).exp()
                      * std::declval<const StateMatrix>().fa() ) t_laplace_result;
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

      //! Shut to open transitions.
      t_laplace_result fa(t_real t) const 
        { return (t*StateMatrix::ff()).exp()*StateMatrix::fa(); }
      //! Open to shut transitions.
      t_laplace_result af(t_real t) const 
        { return (t*StateMatrix::aa()).exp()*StateMatrix::af(); }

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
