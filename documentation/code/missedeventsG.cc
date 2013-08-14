#include <iostream>
#include <exception>

#include <likelihood/missed_eventsG.h>
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
  DCProgs::t_real const tau(1e-4); // in seconds

  // Create eG from prior knowledge of roots
  DCProgs::DeterminantEq determinant_eq(qmatrix, tau);
  std::vector<DCProgs::Root> af_roots{
    { /*root=*/ -3045.285776037674,   /*multiplicity=*/ 1}, 
    { /*root=*/ -162.92946543451328,  /*multiplicity=*/ 1}
  };
  std::vector<DCProgs::Root> fa_roots{
    { /*root=*/ -17090.192769236815,      /*multiplicity=*/ 1},
    { /*root=*/  -2058.0812921673496,     /*multiplicity=*/ 1},
    { /*root=*/     -0.24356535498785126, /*multiplicity=*/ 1}
  };
  DCProgs::MissedEventsG eG_from_roots( determinant_eq, af_roots, 
                                        determinant_eq.transpose(), fa_roots );


  // Create eG by giving home-made root-finding function.
  auto find_roots = [](DCProgs::DeterminantEq const &_det) {
    return DCProgs::find_roots(_det, 1e-12, 1e-12, 100, DCProgs::quiet_nan, DCProgs::quiet_nan);
  };
  DCProgs::MissedEventsG eG_from_func(qmatrix, tau, find_roots);
 
  // Create eG automaticallye
  DCProgs::MissedEventsG eG_automatic(qmatrix, tau);
 

  // Checks the three initialization are equivalent
  for(DCProgs::t_real t(tau);  t < 10*tau; t += tau * 0.1) {

    if(    ((eG_from_roots.af(t) - eG_from_func.af(t)).array().abs() > 1e-8).any()  
        or ((eG_from_roots.fa(t) - eG_from_func.fa(t)).array().abs() > 1e-8).any() )
      throw DCProgs::errors::Runtime("root != func");

    if(    ((eG_from_roots.af(t) - eG_automatic.af(t)).array().abs() > 1e-8).any() 
        or ((eG_from_roots.fa(t) - eG_automatic.fa(t)).array().abs() > 1e-8).any() )
      throw DCProgs::errors::Runtime("root != automatic");
  }

  return 0;
}
