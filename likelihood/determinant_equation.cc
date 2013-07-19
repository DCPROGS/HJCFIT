#include <DCProgsConfig.h>

#include <sstream>

#include "determinant_equation.h"

namespace DCProgs {

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, DeterminantEq const & _self) {
    
    return _stream << "Determinant equation:\n"
                   << "=====================\n\n" 
                   << "  * Transition Rate matrix:\n" << numpy_io(_self.get_qmatrix().matrix) << "\n"
                   << "  * Number of 'A' states: " << _self.get_nopen() << "\n"
                   << "  * Resolution time tau: " << _self.get_tau() << "\n"
                   << "  * FF eigenvalues: " << _self.get_ff_eigenvalues().transpose() << "\n";
  }
}
