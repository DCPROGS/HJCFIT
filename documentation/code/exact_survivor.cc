#include <iostream>

#include <likelihood/exact_survivor.h>
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

  DCProgs::ExactSurvivor survivor(qmatrix, 1e-4);

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

  std::cout << "AF recusion matrices\n"
               "--------------------\n\n";
  for(DCProgs::t_uint i(0); i < 5; ++i)
    for(DCProgs::t_uint m(1); m < 3; ++m)
      for(DCProgs::t_uint l(0); l <= m; ++l) 
        std::cout << "  * C_{" << i << m << l << "}:\n    "
                  << DCProgs::numpy_io(survivor.recursion_af(i, m, l), "    ") << "\n\n";

  return 0;
}


