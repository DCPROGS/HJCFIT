#include <iostream>
#include <likelihood/likelihood.h>

int main() {

  DCProgs::t_Bursts bursts{
     {0.1, 0.2, 0.1},                  /* 1st burst */
     {0.2},                            /* 2nd burst */
     {0.15, 0.16, 0.18, 0.05, 0.1}     /* 3rd burst */
  };
  
  DCProgs::Log10Likelihood likelihood(bursts, /*nopen=*/2, /*tau=*/1e-2, /*tcrit=*/-1);
 
  std::cout << likelihood << std::endl;

  DCProgs::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;

  DCProgs::t_real const result = likelihood(matrix);

  std::cout << "Computation: " << result << std::endl;
 
  return 0;
}
