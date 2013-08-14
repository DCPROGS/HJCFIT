#include <iostream>

#include <likelihood/determinant_equation.h>
#include <likelihood/root_finder.h>
 
int main() {

  // Define parameters.
  DCProgs::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;
  DCProgs::QMatrix qmatrix(matrix, /*nopen=*/2);

  // Create determinant using a QMatrix and a matrix+nopen.
  DCProgs::DeterminantEq det0(qmatrix, 1e-4);
  DCProgs::DeterminantEq det1(matrix, 2, 1e-4);

  std::cout << det0 << "\n\n" << det1 << "\n";

  if( std::abs(det0(0) - det1(0)) > 1e-6 
      or std::abs(det0(-1) - det1(-1)) > 1e-6
      or std::abs(det0(-1e2) - det1(-1e2)) > 1e-6)
    throw DCProgs::errors::Runtime("instanciations differ.");

  if( std::abs( det0(-3045.285776037674) ) > 1e-6 * 3e3
      or std::abs( det0(-162.92946543451328) ) > 1e-6 * 2e2 )
    throw DCProgs::errors::Runtime("Roots are not roots.");

  DCProgs::DeterminantEq transpose = det0.transpose();
  if( std::abs( transpose(-17090.192769236815) ) > 1e-6 * 2e5
      or std::abs( transpose(-2058.0812921673496) ) > 1e-6 * 2e3
      or std::abs( transpose(-0.24356535498785126) ) > 1e-6  )
    throw DCProgs::errors::Runtime("Roots are not roots.");

  std::cout << "  * H(0):\n" << DCProgs::numpy_io(det0.H(0)) << "\n\n"
            << "  * H(-1e2):\n" << DCProgs::numpy_io(det0.H(-1e2)) << "\n\n";

  std::cout << "  * d[sI-H(s)]/ds for s=0:\n" << DCProgs::numpy_io(det0.s_derivative(0)) << "\n\n"
            << "  * d[sI-H(s)]/ds for s=-1e2:\n" << DCProgs::numpy_io(det0.s_derivative(-1e2)) << "\n\n";

  return 0;
}
