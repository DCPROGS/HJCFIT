#ifndef DCPROGS_STATE_MATRIX_H
#define DCPROGS_STATE_MATRIX_H

#include <DCProgsConfig.h>
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
      { return matrix.topRightCorner(nopen, nopen); }
    //! Shut to open transitions.
    auto fa() const -> decltype( matrix.bottomLeftCorner(nopen, nopen) ) 
      { return matrix.bottomLeftCorner(nopen, nopen); }
    //! Shut to shut transitions.
    auto ff() const -> decltype( matrix.bottomRightCorner(nopen, nopen) ) 
      { return matrix.bottomRightCorner(nopen, nopen); }
  };
}

#endif
