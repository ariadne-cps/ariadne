/***************************************************************************
 *            exceptions.h
 *
 *  Copyright  2005-6  Pieter Collins, Alberto Casagrande
 *  Email  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file exceptions.h
 *  \brief Exceptions, error handling and assertions.
 */

#ifndef _ARIADNE_EXCEPTION_H
#define _ARIADNE_EXCEPTION_H

#include <stdexcept>
#include <iosfwd>

namespace Ariadne {


/*! \brief The method or function has not been implemented yet. */
class NotImplemented : public std::logic_error {
 public:
  NotImplemented(const std::string& str) : std::logic_error(str) { }
};

/*! \brief The implementation of a virtual method has been deferred to a derived class. 
 *  However, the method is not critical for the functioning of the class, so it is not
 *  declared pure virtual. 
 */
class DeferredImplementation : public std::logic_error {
 public:
  DeferredImplementation(const std::string& str) : std::logic_error(str) { }
};

/*! \brief The method or function has been deprecated. */
class Deprecated : public std::logic_error {
 public:
  Deprecated(const std::string& str) : std::logic_error(str) { }
};


class invalid_input : public std::runtime_error {
 public:
  invalid_input(const std::string& str) : std::runtime_error(str) { }
};

struct IncompatibleSizes : public std::runtime_error {
  IncompatibleSizes(const std::string& str) : std::runtime_error(str) { }
};

struct IncompatibleDimensions : public std::runtime_error {
  IncompatibleDimensions(const std::string& str) : std::runtime_error(str) { }
};

struct InvalidIndex : public std::runtime_error {
  InvalidIndex(const std::string& str) : std::runtime_error(str) { }
};

struct InvalidGenerators : public std::runtime_error {
  InvalidGenerators(const std::string& str) : std::runtime_error(str) { }
};

struct InvalidNumberOfGenerators : public std::runtime_error {
  InvalidNumberOfGenerators(const std::string& str) : std::runtime_error(str) { }
};

struct IncompatibleGrid : public std::runtime_error {
  IncompatibleGrid(const std::string& str) : std::runtime_error(str) { }
};
  
struct UnboundedSet : public std::runtime_error {
  UnboundedSet(const std::string& str) : std::runtime_error(str) { }
};
  
template<class S1, class S2>
void check_grid(const S1& s1, const S2& s2, const char* where="") {
  if(s1.grid()!=s2.grid()) {
    throw IncompatibleGrid(where);
  }
}
    
template<class S1, class S2>
void check_dimension(const S1& s1, const S2& s2, const char* where="") {
  if(s1.dimension()!=s2.dimension()) {
    throw IncompatibleDimensions(where);
  }
}
    
template<class S1, class V2>
void check_dimension_size(const S1& s1, const V2& v2, const char* where="") {
  if(s1.dimension()!=v2.size()) {
    throw IncompatibleDimensions(where);
  }
}
    
template<class R1, class R2>
void check_dimension(const Geometry::Point<R1>& pt1, const LinearAlgebra::Vector<R2>& v2, const char* where="") {
  if(pt1.dimension()!=v2.size()) {
    throw IncompatibleDimensions(where);
  }
}
    
template<class S1>
void check_index(const S1& s1, const dimension_type& i2, const char* where="") {
  if(i2>=s1.dimension()) {
    throw InvalidIndex(where);
  }
}

template<class S>
void check_bounded(const S& s, const char* where="") {
  if(!s.bounded()) {
    throw UnboundedSet(where);
  }
}

template<class M, class S>
void check_argument_dimension(const M& m1, const S& s2, const char* where="") {
  if(m1.argument_dimension()!=s2.dimension()) {
    throw IncompatibleDimensions(where);
  }
}

template<class M, class S>
void check_result_dimension(const M& m1, const S& s2, const char* where="") {
  if(m1.result_dimension()!=s2.dimension()) {
    throw IncompatibleDimensions(where);
  }
}


namespace Evaluation {

/*! \brief %Base class for exceptions in the Evaluation module. */
class EvaluationException
  : public std::exception 
{
 public:
  EvaluationException(const std::string& s) : _what(s) { }
  ~EvaluationException() throw () { }
  const char* what() const throw () { return this->_what.c_str(); }
 private:
  std::string _what;
};

} // namespace Evaluation

} // namespace Ariadne

#endif /* _ARIADNE_EXCEPTION_H */
