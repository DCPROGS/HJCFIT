#include <iostream>

#include <likelihood/approx_survivor.h>
#include <likelihood/errors.h>
 
int main() {

  // Define parameters.
  DCProgs::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;
  DCProgs::QMatrix qmatrix(matrix, /*nopen=*/2);

  DCProgs::ApproxSurvivor survivor(qmatrix, 1e-4);

  std::cout << survivor << std::endl;

  std::cout << "AF values\n"
               "---------\n\n";
  std::cout << "  * at time t=" << 1e-4 <<":\n    "
            << DCProgs::numpy_io(survivor.af(1e-4), "    ") << "\n"
            << "  * at time t=" << 1.5e-4 <<":\n    " 
            << DCProgs::numpy_io(survivor.af(1.5e-4), "    ") << "\n"
            << "  * at time t=" << 2.0e-4 <<":\n    " 
            << DCProgs::numpy_io(survivor.af(2.0e-4), "    ") << "\n"
            << "  * at time t=" << 2.5e-4 <<":\n    "
            << DCProgs::numpy_io(survivor.af(2.5e-4), "    ") << "\n\n";

  std::cout << " * Exponents: ";
  for(DCProgs::t_uint i(0); i < survivor.nb_af_components(); ++i) 
    std::cout << std::get<1>(survivor.get_af_components(i)) << " ";
  std::cout << std::endl;

  return 0;
}


