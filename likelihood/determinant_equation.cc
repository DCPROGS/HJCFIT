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
