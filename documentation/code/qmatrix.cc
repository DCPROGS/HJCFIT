#include <iostream>
#include <typeinfo>

#include <exception>
#include <likelihood/qmatrix.h>
 

int main() {

  DCProgs::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;

  DCProgs::QMatrix qmatrix(matrix, /*nopen=*/2);
 
  std::cout << qmatrix << std::endl;

  if(qmatrix.nopen != 2) return 1;
  if(qmatrix.nshut() != 3) return 1;
  if( ((qmatrix.matrix - matrix).array().abs() > 1e-12).any() ) return 1;
  
  // Compute sum over rows, row by row.
  for(DCProgs::t_int i(0); i < qmatrix.matrix.rows(); ++i) 
    std::cout << "sum(row[" << i << "]): " << qmatrix.matrix.row(i).sum() << std::endl;

  // Compute sum over rows, but let eigen do it.
  std::cout << "Sum over rows: " << qmatrix.matrix.rowwise().sum().transpose() << std::endl;
  
  // Check that sum over rows returns zero.
  if( (qmatrix.matrix.rowwise().sum().array().abs() > 1e-10).any() ) throw std::exception();
  
  // Gets AA block.
  std::cout << "AA block:\n---------\n" << qmatrix.aa() << std::endl;
                
  // aa returns a wrapper around the matrix block. Hence memory is the same.
  if(&(qmatrix.aa()(0,0)) != &(qmatrix(0, 0))) throw std::exception();
 
  return 0;
}
