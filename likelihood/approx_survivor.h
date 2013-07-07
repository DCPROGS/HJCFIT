#ifndef DCPROGS_LIKELIHOOD_APPROX_SURVIVOR_H
#define DCPROGS_LIKELIHOOD_APPROX_SURVIVOR_H

#include <DCProgsConfig.h>

#include <functional>
#include <memory>

#include "asymptotes.h"
#include "root_finder.h"

namespace DCProgs {

  //! \brief Implementation of asymptotic missed events.
  //! \details This object merely puts together two Asymptotes objects.
  class MSWINDOBE ApproxSurvivor {
    public:
      //! Type of the function used to minimize roots.
      typedef std::function<std::vector<Root>(DeterminantEq const &)> t_RootFinder;
      //! Initializes approx survivor functor.
      //! \param[in] _matrix: Partitionned matrix with open states in top left corner.
      //! \param[in] _tau: Missed event cutoff time.
      ApproxSurvivor(DeterminantEq const &_af, std::vector<Root> const &_roots_af, 
                     DeterminantEq const &_fa, std::vector<Root> const &_roots_fa );
      //! Initializes approx survivor functor.
      ApproxSurvivor(StateMatrix const &_matrix, t_real _tau, t_RootFinder const &_findroots);

      //! Open to close transitions 
      t_rmatrix af(t_real _t) const { return asymptotes_af_->operator()(_t); }
      //! Close to open transitions
      t_rmatrix fa(t_real _t) const { return asymptotes_fa_->operator()(_t); }
    protected:
#     ifndef HAS_CXX11_UNIQUE_PTR
        //! Type of the pointers holding recursion interfaces.
        typedef std::auto_ptr<Asymptotes> t_AsymptotesPtr;
#     else
        //! Type of the pointers holding recursion interfaces.
        typedef std::unique_ptr<Asymptotes> t_AsymptotesPtr;
#     endif
      //! Pointer to AF recursion interface
      t_AsymptotesPtr asymptotes_af_;
      //! Pointer to FA recursion interface
      t_AsymptotesPtr asymptotes_fa_;
  };
}
#endif 

