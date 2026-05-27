#include <iostream>
#include <exception>

#include <likelihood/missed_eventsG.h>
#include <likelihood/root_finder.h>
 
int main() {

  // Define parameters.
  HJCFIT::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;
  HJCFIT::QMatrix qmatrix(matrix, /*nopen=*/2);
  HJCFIT::t_real const tau(1e-4); // in seconds

  // Create eG from prior knowledge of roots
  HJCFIT::DeterminantEq determinant_eq(qmatrix, tau);
  std::vector<HJCFIT::Root> af_roots{
    { /*root=*/ -3045.285776037674,   /*multiplicity=*/ 1}, 
    { /*root=*/ -162.92946543451328,  /*multiplicity=*/ 1}
  };
  std::vector<HJCFIT::Root> fa_roots{
    { /*root=*/ -17090.192769236815,      /*multiplicity=*/ 1},
    { /*root=*/  -2058.0812921673496,     /*multiplicity=*/ 1},
    { /*root=*/     -0.24356535498785126, /*multiplicity=*/ 1}
  };
  HJCFIT::MissedEventsG eG_from_roots( determinant_eq, af_roots, 
                                        determinant_eq.transpose(), fa_roots );


  // Create eG by giving home-made root-finding function.
  auto find_roots = [](HJCFIT::DeterminantEq const &_det) {
    return HJCFIT::find_roots(_det, 1e-12, 1e-12, 100, HJCFIT::quiet_nan, HJCFIT::quiet_nan);
  };
  HJCFIT::MissedEventsG eG_from_func(qmatrix, tau, find_roots);
 
  // Create eG automaticallye
  HJCFIT::MissedEventsG eG_automatic(qmatrix, tau);
 

  // Checks the three initialization are equivalent
  for(HJCFIT::t_real t(tau);  t < 10*tau; t += tau * 0.1) {

    if(    ((eG_from_roots.af(t) - eG_from_func.af(t)).array().abs() > 1e-8).any()  
        or ((eG_from_roots.fa(t) - eG_from_func.fa(t)).array().abs() > 1e-8).any() )
      throw HJCFIT::errors::Runtime("root != func");

    if(    ((eG_from_roots.af(t) - eG_automatic.af(t)).array().abs() > 1e-8).any() 
        or ((eG_from_roots.fa(t) - eG_automatic.fa(t)).array().abs() > 1e-8).any() )
      throw HJCFIT::errors::Runtime("root != automatic");
  }

  return 0;
}
