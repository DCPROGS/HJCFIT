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

#ifndef DCPROGS_STATE_MATRIX_H
#define DCPROGS_STATE_MATRIX_H

#include <DCProgsConfig.h>
#include <ostream>
#include <tuple>
#include "errors.h"

//! General namespace for all things DCProgs.
namespace DCProgs {

  //! \brief State matrix that can  be partitioned into open/shut states.
  //! \details In practice, this is a two tuple with some helper functions to get corners.
  struct MSWINDOBE QMatrix {
 
    //! Number of open states.
    t_uint nopen; 
    //! The matrix itself.
    t_rmatrix matrix; 
 
    //! Constructor
    QMatrix() : nopen(0), matrix(0,0) {}
    //! Constructor
    template<class T>
      QMatrix(Eigen::DenseBase<T> const &_c, t_uint _nopen = 0) : nopen(_nopen), matrix(_c) {}

    //! Returns rate at location
    t_rmatrix::Scalar const & operator()( t_rmatrix::Index const &_i,
                                          t_rmatrix::Index const &_j ) const {
      return matrix(_i, _j);
    }
    
    //! Returns rate at location
    t_rmatrix::Scalar & operator()( t_rmatrix::Index const &_i, t_rmatrix::Index const &_j ) {
      return matrix(_i, _j);
    }
  
    //! Open to open transitions.
    Eigen::Block<t_rmatrix> aa() { return matrix.topLeftCorner(nopen, nopen); }
    //! Open to shut transitions.
    Eigen::Block<t_rmatrix> af() { return matrix.topRightCorner(nopen, nshut()); }
    //! Shut to open transitions.
    Eigen::Block<t_rmatrix> fa() { return matrix.bottomLeftCorner(nshut(), nopen); }
    //! Shut to shut transitions.
    Eigen::Block<t_rmatrix> ff() 
      { return matrix.bottomRightCorner(nshut(), nshut()); }
    //! Open to open transitions.
    Eigen::Block<t_rmatrix const> aa() const 
      { return matrix.topLeftCorner(nopen, nopen); }
    //! Open to shut transitions.
    Eigen::Block<t_rmatrix const> af() const 
      { return matrix.topRightCorner(nopen, nshut()); }
    //! Shut to open transitions.
    Eigen::Block<t_rmatrix const> fa() const 
      { return matrix.bottomLeftCorner(nshut(), nopen); }
    //! Shut to shut transitions.
    Eigen::Block<t_rmatrix const> ff() const 
      { return matrix.bottomRightCorner(nshut(), nshut()); }

    //! Number of shut states
    t_uint nshut() const { return static_cast<t_uint>(matrix.cols()) - nopen; }

    //! \brief Returns transpose of this \f$Q\f$-matrix.
    //! \details A states become F states, and F states become A states, and the matrix is
    //! transposed such that the new AA block is in the top left corner.
    QMatrix transpose() const;
 
    //! \brief Computes eigenvalues and eigenvectors
    //! \details Solves the *transpose* eigenproblem \f$\phi = \phi\cdot\mathcal{Q}\f$.
    std::tuple<t_cvector, t_cmatrix> eigenstuff() const;
  };

  //! Dumps object to stream.
  MSWINDOBE std::ostream & operator<< (std::ostream &_stream, QMatrix const &_mat);

  //! Verify that a QMatrix is not to large to be used in a stack allocated
  //! likelihood, exact_survivor or asymptotes calculation.
  void verify_qmatrix(QMatrix const &_qmatrix);
}

#endif
