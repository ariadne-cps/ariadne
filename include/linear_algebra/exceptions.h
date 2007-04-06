/***************************************************************************
 *            linear_algebra/exceptions.h
 *
 *  Copyright  2005-7  Pieter Collins, Alberto Casagrande
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

#ifndef ARIADNE_LINEAR_ALGEBRA_EXCEPTIONS_H
#define ARIADNE_LINEAR_ALGEBRA_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

#include "../throw.h"
#include "../base/types.h"
#include "../base/exceptions.h"

namespace Ariadne {
  namespace LinearAlgebra {
    //@{ \name Exceptions
    /*! \brief The sizes of the operands are incompatible. */
    struct IncompatibleSizes : public std::runtime_error {
      IncompatibleSizes(const std::string& str) : std::runtime_error(str) { }
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

#define ARIADNE_CHECK_SIZE(vec,sz,func) \
  { if((vec).size()!=sz) { ARIADNE_THROW(IncompatibleSizes,func,#vec"="<<vec<<", "#sz"="<<sz); } }

#define ARIADNE_CHECK_EQUAL_SIZES(vec1,vec2,func) \
  { if((vec1).size()!=vec2.size()) { ARIADNE_THROW(IncompatibleSizes,func,#vec1"="<<vec1<<", "#vec2"="<<vec2); } }

#define ARIADNE_CHECK_MATRIX_EQUAL_SIZES(mx1,mx2,func) \
  { if((mx1).number_of_rows()!=mx2.number_of_rows() || (mx1).number_of_columns()!=mx2.number_of_columns()) { \
      ARIADNE_THROW(IncompatibleSizes,func,#mx1"="<<mx1<<", "#mx2"="<<mx2); } }

#define ARIADNE_CHECK_SQUARE(mx,func) \
  { if((mx).number_of_rows()!=(mx).number_of_columns()) { ARIADNE_THROW(IncompatibleSizes,func,#mx"="<<mx<<" is not square"); } }

#define ARIADNE_CHECK_INDEX(vec,ind,func)                                  \
  { if((vec).size()>=ind) { ARIADNE_THROW(InvalidIndex,func,#vec"="<<vec<<", "#ind"="<<ind); } }

#endif /* ARIADNE_LINEAR_ALGEBRA_EXCEPTIONS_H */
