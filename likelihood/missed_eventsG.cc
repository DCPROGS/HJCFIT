#include <iostream>
#include "missed_eventsG.h"

namespace DCProgs {

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, MissedEventsG const & _self) {
    
    return _stream << "Missed-Event Likelihood:\n"
                   << "=======================\n\n" 
                   << "  * Transition Rate matrix:\n" << numpy_io(_self.get_qmatrix().matrix) << "\n"
                   << "  * Number of 'A' states: " << _self.get_qmatrix().nopen << "\n"
                   << "  * Resolution time tau: " << _self.get_tau() << "\n"
                   << "  * Exact events computed for: t < " << _self.get_nmax() << " tau\n";
  }
  
  // CHS matrices \f$H_{FA}\f$
  t_rmatrix CHS_matrix_Hfa(MissedEventsG const &_g, t_real _tcrit) {
    return partial_CHS_matrix(*_g.asymptotes_fa_, _g.get_tau(), _tcrit) * _g.get_fa_factor();
  }
  // CHS matrices \f$H_{FA}\f$
  t_rmatrix CHS_matrix_Haf(MissedEventsG const &_g, t_real _tcrit) {
    return partial_CHS_matrix(*_g.asymptotes_af_, _g.get_tau(), _tcrit) * _g.get_af_factor();
  }
}
