#include <iostream>
#include <exception>

#include <likelihood/idealG.h>
#include <likelihood/errors.h>
 
int main() {

  // Define parameters.
  HJCFIT::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;
  HJCFIT::QMatrix qmatrix(matrix, /*nopen=*/2);

  HJCFIT::IdealG idealG(qmatrix);

  std::cout << idealG << std::endl;

  HJCFIT::t_rmatrix const idealG_fa = (2e-4*qmatrix.ff()).exp()*qmatrix.fa();
  if( ((idealG.fa(2e-4) - idealG_fa).array().abs() > 1e-8).any() )
    throw HJCFIT::errors::Runtime("Not so ideal idealG");

  HJCFIT::t_rmatrix const inversion = -0.5 * HJCFIT::t_rmatrix::Identity(2, 2) - qmatrix.aa();
  if( ((inversion * idealG.laplace_af(-0.5) - qmatrix.af()).array().abs() > 1e-8).any() )
    throw HJCFIT::errors::Runtime("Not so ideal idealG");

  return 0;
}


