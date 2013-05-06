#ifndef DCPROGS_STATE_MATRIX_H
#define DCPROGS_STATE_MATRIX_H

#include <DCProgsConfig.h>
#include <tuple>
#include "errors.h"

//! General namespace for all things DCProgs.
namespace DCProgs {

  //! \brief State matrix that can  be partitioned into open/shut states.
  //! \details In practice, this is a two tuple with some helper functions to get corners.
  struct StateMatrix {


    //! Number of open states.
    t_int nopen; 
    //! The matrix itself.
    t_rmatrix matrix; 

    //! Constructor
    StateMatrix() : matrix(0,0), nopen(0) {}
    //! Constructor
    template<class T>
      StateMatrix   (Eigen::DenseBase<T> const &_c, t_int _nopen = 0)
                  : matrix(_c), nopen(_nopen) {}
  
    //! Open to open transitions.
    auto aa() const -> decltype( matrix.topLeftCorner(nopen, nopen) ) 
      { return matrix.topLeftCorner(nopen, nopen); }
    //! Open to shut transitions.
    auto af() const -> decltype( matrix.topRightCorner(nopen, nopen) ) 
      { return matrix.topRightCorner(nopen, matrix.rows() - nopen); }
    //! Shut to open transitions.
    auto fa() const -> decltype( matrix.bottomLeftCorner(nopen, nopen) ) 
      { return matrix.bottomLeftCorner(matrix.rows() - nopen, nopen); }
    //! Shut to shut transitions.
    auto ff() const -> decltype( matrix.bottomRightCorner(nopen, nopen) ) 
      { return matrix.bottomRightCorner(matrix.rows() - nopen, matrix.rows() - nopen); }

    //! \brief Computes eigenvalues and eigenvectors
    //! \details Solves the *transpose* eigenproblem \f$\phi = \phi\cdot\mathcal{Q}\f$.
    std::tuple<t_cvector, t_cmatrix> eigenstuff() {
      Eigen::EigenSolver<t_rmatrix> eigsolver(matrix.transpose());
      if(eigsolver.info() != Eigen::Success) 
        throw errors::Mass("Could not solve eigenvalue problem.");
      t_cvector const eigs = eigsolver.eigenvalues();
      t_cmatrix const vecs = eigsolver.eigenvectors();
      return std::make_tuple(eigs, vecs.transpose());
    }
  };
}

#endif
