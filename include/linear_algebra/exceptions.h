/***************************************************************************
 *            linear_algebra.except.h
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
 
/*! \file linear_algebra/exceptions.h
 *  \brief Exceptions, error handling and assertions for the Linear Algebra module.
 */

#ifndef ARIADNE_LINEAR_ALGEBRA_EXCEPT_H
#define ARIADNE_LINEAR_ALGEBRA_EXCEPT_H

#include <stdexcept>
#include <iosfwd>

#include "base/types.h"

namespace Ariadne {
  namespace LinearAlgebra {
    //@{ \name Exceptions
    /*! \brief The sizes of the operands are incompatible. */
    struct IncompatibleSizes : public std::runtime_error {
      IncompatibleSizes(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief The index to an arrayed object was invalid. */
    struct InvalidIndex : public std::runtime_error {
      InvalidIndex(const std::string& str) : std::runtime_error(str) { }
    };
    //@}

    template<class V1> inline
    void check_size(const V1& v1, const size_type& n2, const char* where="") {
      if(v1.size()!=n2) { throw IncompatibleSizes(where); }
    }

    template<class V1, class V2> inline
    void check_equal_sizes(const V1& v1, const V2& v2, const char* where="") {
      if(v1.size()!=v2.size()) { throw IncompatibleSizes(where); }
    }

    template<class Mx> inline
    void check_square(const Mx& A, const char* where="") {
      if(A.number_of_rows()!=A.number_of_columns()) { throw IncompatibleSizes(where); }
    }

    template<class S1> inline
    void check_index(const S1& s1, const size_type& i2, const char* where="") {
      if(i2>=s1.size()) { throw InvalidIndex(where); }
    }

  }
}

#endif /* ARIADNE_LINEAR_ALGEBRA_EXCEPT_H */
