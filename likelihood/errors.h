#ifndef DCPROGS_ERRORS_H
#define DCPROGS_ERRORS_H
#include "DCProgsConfig.h"
#include <exception>
#include <stdexcept>

namespace DCProgs {

  //! Exceptions of DCprogs
  namespace errors {
    
    //! All explicit DCProgs exception derive from this, for easy catching,
    class Root : public std::exception {};
    //! Input size error
    class Domain : public Root, virtual public std::domain_error {
      public:
        explicit Domain(char const *_message) DC_NOEXCEPT : Root(), std::domain_error(_message) {};
        virtual char const* what() const DC_NOEXCEPT { return this->std::domain_error::what(); }
    };
  }
}
#endif
