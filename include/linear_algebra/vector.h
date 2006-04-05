/***************************************************************************
 *            vector.h
 *
 *  Mon May  3 12:31:15 2004
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file vector.h
 *  \brief Vectors and vector operations.
  */

#ifndef _ARIADNE_VECTOR_H
#define _ARIADNE_VECTOR_H 

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include "../base/basic_type.h"
#include "../base/numerical_type.h"
#include "../base/interval.h"

namespace boost { 
  namespace numeric { 
    namespace ublas {

      template <typename R>
      std::ostream&
      operator<<(std::ostream& os, const vector<R>& v)
      {
        os << "[";
        if(v.size()>0) {
          os << v[0];
        }
        for(uint i=1; i!=v.size(); ++i) {
          os << "," << v[i];
        }
        os << "]";
        return os;
      }

    }
  }
}


namespace Ariadne {
  namespace LinearAlgebra {

    using boost::numeric::ublas::vector;
    using boost::numeric::ublas::vector_expression;
        
/*
    template<typename R>
    class Vector
      : public boost::numeric::ublas::vector<R> 
    {
      typedef boost::numeric::ublas::vector<R> Base;
     public:
      template<typename E>
      Vector(const boost::numeric::ublas::vector_expression<E>& v) : Base(v().size()) {
        for(dimension_type i=0; i!=this->size(); ++i) {
          (*this)[i]=v()[i];
        }
      }
    };
*/
    
    template <typename R>
    inline vector<R> zero_vector(dimension_type dim) {
      vector<R> v(dim);
      for (dimension_type i=0; i<dim; ++i) {
        v(i)=0;
      }
      return v;
    }
    
    template <typename R>
    inline bool linear_multiple(const vector<R> u, const vector<R>& v) {
      assert(u.dimension()==v.dimension());
      R multiple=0;
      for (dimension_type i=0; i<u.dimension(); ++i) {
        if(u(i)==0 && multiple!=0 && v(i)!=0) {
          return false;
        }
        if(multiple==0) {
          multiple=v(i)/u(i);
        }
        else {
          if(multiple*u(i)!=v(i)) {
            return false;
          }
        }
      }
      return true;
    }
    
    template<typename R>
    inline
    R
    norm(const vector<R>& v) {
      R norm=R(0);
      for (size_type i=0; i<v.size(); ++i) {
        norm=std::max(norm,R(abs(v(i))));
      }
      return norm;
    }
    
    template <typename R>
    inline Integer common_denominator(const vector<R>& b) 
    {
      Integer denom=1;
      for (dimension_type i=0; i< b.size(); ++i) {
        denom=lcm( denom, denominator(b(i)) );
      }
      return denom;
    }

    template <typename R>
    inline R norm_infinite(const vector<R>& b) 
    {
      R norm=0.0;
      for (size_t i=0; i< b.size(); i++) {
        if (abs(b[i])>norm) norm=abs(b[i]);
      }
      return norm;
    }

  }
}




#endif /* _ARIADNE_VECTOR_H */
