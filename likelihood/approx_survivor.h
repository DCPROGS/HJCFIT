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
      //! Initializes approximate survivor functor.
      //! \param[in] _af: Determinantal equation for open->shut transitions
      //! \param[in] _roots_af: Roots of _af equation
      //! \param[in] _fa: Determinantal equation for shut->open transitions
      //! \param[in] _roots_fa: Roots of _fa equation
      ApproxSurvivor(DeterminantEq const &_af, std::vector<Root> const &_roots_af, 
                     DeterminantEq const &_fa, std::vector<Root> const &_roots_fa );
      //! Initializes approx survivor functor.
      //! \param[in] _matrix: Transition matrix
      //! \param[in] _tau: resolution/max length missed events
      //! \param[in] _findroots: A functor with which to find all roots.
      //!                        This function should take a DeterminantEq as its sole argument and
      //!                        return a std::vector<RootIntervals>
      ApproxSurvivor(QMatrix const &_matrix, t_real _tau, t_RootFinder const &_findroots);

      //! Open to close transitions 
      t_rmatrix af(t_real _t) const { return asymptotes_af_->operator()(_t); }
      //! Close to open transitions
      t_rmatrix fa(t_real _t) const { return asymptotes_fa_->operator()(_t); }
      //! Number of exponential components for af
      t_int nb_af_components() const { return asymptotes_af_->size(); }
      //! Number of exponential components for fa
      t_int nb_fa_components() const { return asymptotes_fa_->size(); }
      //! AF exponential components
      Asymptotes::t_MatrixAndRoot const & get_af_components(t_int i) const {
        return (*asymptotes_af_)[i]; 
      }
      //! FA exponential components
      Asymptotes::t_MatrixAndRoot const & get_fa_components(t_int i) const {
        return (*asymptotes_fa_)[i]; 
      }
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

  //! Dumps Survivor function equation to stream
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, ApproxSurvivor const &_self);
}
#endif 

