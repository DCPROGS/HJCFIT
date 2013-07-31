/***********************
    DCProgs computes missed-events likelihood as described in
    Hawkes, Jalali and Colquhoun (1990, 1992)

    Copyright (C) 2013  University College London

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
************************/

#ifndef DCPROGS_ERRORS_H
#define DCPROGS_ERRORS_H
#include "DCProgsConfig.h"
#include <string>
#include <exception>
#include <stdexcept>

namespace DCProgs {

  //! Exceptions of DCprogs
  namespace errors {
    
    //! All explicit DCProgs exception derive from this, for easy catching,
    class Root : public std::exception {};
    //! Math (convergence, domain, etc) error
    class Math : public Root { };
    //! Math error which carries a message.
    class Mass : public Math {
      public:
        Mass(std::string const &_message) noexcept : Math() {
          try { message_ = _message; }
          catch(...) { try { message_ = ""; } catch(...) {} }
        }
        virtual char const * what() const noexcept {
          try { return message_.c_str(); } catch(...) { return ""; }
        }
        virtual ~Mass() noexcept {};
      protected:
        std::string message_;
    };

    //! Found unexpected complex eigenvalues. 
    class ComplexEigenvalues : public Mass {
      public:
        ComplexEigenvalues(std::string const &_message) noexcept : Mass(_message) {
          try { message_ = _message; }
          catch(...) { try { message_ = ""; } catch(...) {} }
        }
        virtual char const * what() const noexcept {
          try {
            return ("Found complex eigenvalues: " + message_).c_str(); 
          } catch(...) { return ""; }
        }
        virtual ~ComplexEigenvalues() noexcept {};
    };

    //! Found a Not a Number 
    class NaN : public Mass { 
      public:
        NaN(std::string const &_message) noexcept : Mass(_message) {
          try { message_ = _message; }
          catch(...) { try { message_ = ""; } catch(...) {} }
        }
        virtual char const * what() const noexcept {
          try {
            return ("Found Not a Number: " + message_).c_str(); 
          } catch(...) { return ""; }
        }
        virtual ~NaN() noexcept {};
    };

    //! Input error to a math problem
    class Domain : public Math, virtual public std::domain_error {
      public:
        explicit Domain(char const *_message) noexcept : std::domain_error(_message), Math()  {};
        explicit Domain(std::string const &_message) noexcept : std::domain_error(_message), Math() {};
        virtual char const* what() const noexcept { return this->std::domain_error::what(); }
    };
    //! Index error
    class Index : public Root, public virtual std::out_of_range {
      public:
        explicit Index(char const *_message) noexcept : std::out_of_range(_message), Root() {};
        explicit Index(std::string const &_message) noexcept : std::out_of_range(_message), Root() {};
        virtual char const* what() const noexcept { return this->std::out_of_range::what(); }
    };


    //! Matrix is not invertible
    class NotInvertible : public Domain {
      public:
        NotInvertible   (char const *_message) noexcept
                      : std::domain_error(_message), Domain(_message) {};
        NotInvertible   (std::string const &_message) noexcept 
                      : std::domain_error(_message), Domain(_message) {};
    };

    //! Runtime error which carries a message.
    class Runtime : public Root {
      public:
        Runtime(std::string const &_message) noexcept : Root() {
          try { message_ = _message; }
          catch(...) { try { message_ = ""; } catch(...) {} }
        }
        virtual char const * what() const noexcept {
          try { return message_.c_str(); } catch(...) { return ""; }
        }
        virtual ~Runtime() noexcept {};
      private:
        std::string message_;
    };

    //! NotImplemented error which carries a message.
    class NotImplemented : public Root {
      public:
        NotImplemented(std::string const &_message) noexcept : Root() {
          try { message_ = _message; }
          catch(...) { try { message_ = ""; } catch(...) {} }
        }
        virtual char const * what() const noexcept {
          try { return message_.c_str(); } catch(...) { return ""; }
        }
        virtual ~NotImplemented() noexcept {};
      private:
        std::string message_;
    };


#   ifdef DCPROGS_PYTHON_BINDINGS
      //! Exception thrown in python modules 
      class Python : public Root {
        public:
          Python(std::string const &_message) noexcept : Root() {
            try { message_ = _message; }
            catch(...) { try { message_ = ""; } catch(...) {} }
          }
          virtual char const * what() const noexcept {
            try { return message_.c_str(); } catch(...) { return ""; }
          }
          virtual ~Python() noexcept {};
        private:
          std::string message_;
      };
      //! Exception was thrown by python API.
      class PythonErrorAlreadyThrown : public Python {
        public:
          PythonErrorAlreadyThrown() noexcept: Python("") {};
      };
      //! Exception thrown in python module when converting to C types.
      class PythonTypeError : public Python { 
        public:
          PythonTypeError(std::string const &_message) noexcept: Python(_message) {};
      };
      //! Exception thrown in python module when converting to C types.
      class PythonValueError : public Python {
        public:
          PythonValueError(std::string const &_message) noexcept: Python(_message) {};
      };
#   endif
  }
}
#endif
