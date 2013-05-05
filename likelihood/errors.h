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
        explicit Domain(char const *_message) noexcept : Root(), std::domain_error(_message) {};
        explicit Domain(std::string const &_message) noexcept : Root(), std::domain_error(_message) {};
        virtual char const* what() const noexcept { return this->std::domain_error::what(); }
    };
    //! Matrix is not invertible
    class NotInvertible : public Domain {
      public:
        explicit NotInvertible   (char const *_message) noexcept
                               : Domain(_message), std::domain_error(_message) {};
        explicit NotInvertible   (std::string const &_message) noexcept 
                               : Domain(_message), std::domain_error(_message) {};
    };
  }
}
#endif
