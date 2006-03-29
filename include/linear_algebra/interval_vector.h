/***************************************************************************
 *            interval_vector.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file interval_vector.h
 *  \brief Vectors of intervals.
  */

#ifndef _ARIADNE_INTERVAL_VECTOR_H
#define _ARIADNE_INTERVAL_VECTOR_H 

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "../base/basic_type.h"
#include "../base/numerical_type.h"
#include "../base/interval.h"

namespace boost {
  namespace numeric {
    namespace ublas {
      
      using Ariadne::Interval;
      using Ariadne::size_type;
      
      template <typename Real>
      vector< Interval<Real> > 
      operator+(const vector< Interval<Real> >& iv, const vector<Real> v)
      {
        vector< Interval<Real> > result(v.size());
        for(size_type i=0; i!=result.size(); ++i) {
          result(i)=v(i)+iv(i);
        }
        return result;
      }

      template <typename Real>
      vector< Interval<Real> > 
      operator+(const vector<Real>& v, const vector< Interval<Real> > iv)
      {
        vector< Interval<Real> > result(v.size());
        for(size_type i=0; i!=result.size(); ++i) {
          result(i)=v(i)+iv(i);
        }
        return result;
      }

      template<typename Real>
      inline
      vector< Interval<Real> >
      operator*(const Real& s, const vector< Interval<Real> >& v)
      {
        return Interval<Real>(s)*v;
      }
      
    }
  }
}
  
namespace Ariadne {
  namespace LinearAlgebra {
    template<typename Real>
    inline
    Real
    norm(const vector< Interval<Real> >& v)
    {
      Real norm=Real(0);
      for (size_type i=0; i<v.size(); i++) {
        norm=std::max(norm,Real(abs(v(i).lower())));
        norm=std::max(norm,Real(abs(v(i).upper())));
      }
      return norm;
    }
    
  }
}  

#endif /* _ARIADNE_INTERVAL_VECTOR_H */
