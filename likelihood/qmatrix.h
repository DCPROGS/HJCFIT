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
    t_int nopen; 
    //! The matrix itself.
    t_rmatrix matrix; 
 
    //! Constructor
    QMatrix() : matrix(0,0), nopen(0) {}
    //! Constructor
    template<class T>
      QMatrix(Eigen::DenseBase<T> const &_c, t_int _nopen = 0) : matrix(_c), nopen(_nopen) {}
  
    //! Open to open transitions.
    Eigen::Block<t_rmatrix> aa() { return matrix.topLeftCorner(nopen, nopen); }
    //! Open to shut transitions.
    Eigen::Block<t_rmatrix> af() { return matrix.topRightCorner(nopen, matrix.rows() - nopen); }
    //! Shut to open transitions.
    Eigen::Block<t_rmatrix> fa() { return matrix.bottomLeftCorner(matrix.rows() - nopen, nopen); }
    //! Shut to shut transitions.
    Eigen::Block<t_rmatrix> ff() 
      { return matrix.bottomRightCorner(matrix.rows() - nopen, matrix.rows() - nopen); }
    //! Open to open transitions.
    Eigen::Block<t_rmatrix const> aa() const 
      { return matrix.topLeftCorner(nopen, nopen); }
    //! Open to shut transitions.
    Eigen::Block<t_rmatrix const> af() const 
      { return matrix.topRightCorner(nopen, matrix.rows() - nopen); }
    //! Shut to open transitions.
    Eigen::Block<t_rmatrix const> fa() const 
      { return matrix.bottomLeftCorner(matrix.rows() - nopen, nopen); }
    //! Shut to shut transitions.
    Eigen::Block<t_rmatrix const> ff() const 
      { return matrix.bottomRightCorner(matrix.rows() - nopen, matrix.rows() - nopen); }

    t_int nshut() const { return matrix.cols() - nopen; }

    //! \brief Returns transpose of state matrix.
    //! \details Means A states become F states, and F states become A states, and the partitionned
    //! matrix is transposed such that the new AA block is top left corner.
    QMatrix transpose() const;
 
    //! \brief Computes eigenvalues and eigenvectors
    //! \details Solves the *transpose* eigenproblem \f$\phi = \phi\cdot\mathcal{Q}\f$.
    std::tuple<t_cvector, t_cmatrix> eigenstuff() const;
  };

  //! Dumps object to stream.
  MSWINDOBE std::ostream & operator<< (std::ostream &_stream, QMatrix const &_mat);
}

#endif
