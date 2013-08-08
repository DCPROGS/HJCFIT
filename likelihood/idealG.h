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

#ifndef DCPROGS_LIKELIHOOD_IDEALG
#define DCPROGS_LIKELIHOOD_IDEALG

#include <DCProgsConfig.h>

#include <utility>

#include <unsupported/Eigen/MatrixFunctions>

#include "qmatrix.h"
#include "errors.h"

//! General namespace for all things DCProgs.
namespace DCProgs {

  //! \brief Ideal transition matrix of open and shut intervals
  //! \details Given a transition matrix $Q$ it is possible to figure out the evolution of any given
  //! system. 
  class MSWINDOBE IdealG : protected QMatrix {

    //! Just trying to figure out a complex return type...
    typedef decltype( (t_real(0) * std::declval<const QMatrix>().ff()).exp()
                      * std::declval<const QMatrix>().fa() ) t_time_result;
    public:
      //! Constructor
      IdealG() : QMatrix() {}
      //! \brief Constructor with parameters.
      //! \details Calls set method with input parameters.
      //! \param[in] _matrix Any matrix or matrix expression from Eigen. Will become the transition
      //!                     matrix. Diagonal elements are transformed as explain in set(). Open
      //!                     states should be in the top rows.
      //! \param[in] _nopen Number of open states. 
      //! \throws errors::Domain if input has incorrect values or size.
      template<class T>
        IdealG(Eigen::DenseBase<T> const &_matrix, t_uint _nopen);
      //! Constructor
      IdealG(QMatrix const &_c) : QMatrix(_c) {};
      //! Destructor 
      virtual ~IdealG() {}; 
  
      //! \brief Sets Q matrix and the number of open states.
      //! \details Enforces \f$Q_{ii} = -\sum_{j\neq i} Q_{ij}\f$.
      //!          It is expected that open states are the top rows [0, _nopen].
      void set(t_rmatrix const &_Q, t_uint const &_nopen);
      //! Sets state matrix on which to act.
      void set(QMatrix const &_in) { set(_in.matrix, _in.nopen); }
      //! Gets Q matrix. 
      t_rmatrix const & get_matrix() const { return this->matrix; }
      //! Gets the number of open states
      t_uint get_nopen() const { return this->nopen; }
      //! Gets the number of open states
      t_uint get_nshut() const { return this->nshut(); }

      //! Shut to open transitions.
      t_time_result fa(t_real t) const 
        { return (t*QMatrix::ff()).exp()*QMatrix::fa(); }
      //! Open to shut transitions.
      t_time_result af(t_real t) const 
        { return (t*QMatrix::aa()).exp()*QMatrix::af(); }

      //! Laplace transform of shut to open transitions.
      t_rmatrix laplace_fa(t_real s) const;
      //! Open to shut transitions.
      t_rmatrix laplace_af(t_real t) const;
  };

  template<class T>
    IdealG :: IdealG(Eigen::DenseBase<T> const &_matrix, t_uint _nopen) : QMatrix() {
      try { this->set(_matrix, _nopen); }
      catch(...) {
        this->matrix.resize(0, 0);
        this->nopen = 0;
        throw;
      }
    }

  //! Dumps object to stream.
  MSWINDOBE std::ostream & operator<< (std::ostream &_stream, IdealG const &_mat);
}

#endif
