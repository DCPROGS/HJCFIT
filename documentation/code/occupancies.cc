#include <iostream>

#include <likelihood/missed_eventsG.h>
#include <likelihood/idealG.h>
#include <likelihood/occupancies.h>
 
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

  // Create missed-events G
  DCProgs::MissedEventsG eG(qmatrix, tau);
  // Create ideal G
  DCProgs::IdealG idealG(qmatrix);
  
  DCProgs::t_real const tcritical(5e-3);

  std::cout << "Equilibrium Occupancies\n"
            << "=======================\n\n"
            << "Ideal Likelihood\n"
            << "----------------\n\n"
            << "  * initial: " << DCProgs::occupancies(idealG) << "\n"
            << "  * final: "   << DCProgs::occupancies(idealG, false) << "\n\n\n"
            << "Missed-events Likelihood\n"
            << "------------------------\n\n"
            << "  * initial: " << DCProgs::occupancies(eG) << "\n"
            << "  * final: "   << DCProgs::occupancies(eG, false) << "\n\n\n\n"
            << "CHS Occupancies\n"
            << "===============\n\n"
            << "Missed-events Likelihood\n"
            << "------------------------\n\n"
            << "  * tcritical: " << tcritical << "\n"
            << "  * initial: " << DCProgs::CHS_occupancies(eG, tcritical) << "\n"
            << "  * final: "   << DCProgs::CHS_occupancies(eG, tcritical, false) << "\n";

  return 0;
}
