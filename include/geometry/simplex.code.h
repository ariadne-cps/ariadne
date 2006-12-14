/***************************************************************************
 *            simplex.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "simplex.h"

#include <iostream>
#include <vector>

#include "../base/stlio.h"

namespace Ariadne {
  namespace Geometry {

    template<class R>
    Simplex<R>::Simplex()
      : Polytope<R>(LinearAlgebra::Matrix<R>(0,1))
    {
    }
  
    template<class R>
    Simplex<R>::Simplex(const LinearAlgebra::Matrix<R>& A)
      : Polytope<R>(A)
    {
      check_dimension(*this,A.number_of_columns()-1,__PRETTY_FUNCTION__);
    }
    
    template<class R>
    Simplex<R>::Simplex(const PointList<R>& v)
      : Polytope<R>(v)
    {
      check_dimension(*this,v.size()-1,__PRETTY_FUNCTION__);
    }
    
    template<class R>
    tribool
    Simplex<R>::contains(const Point<R>& pt) const
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      LinearAlgebra::Vector<F> coords=this->coordinates(pt);
      //std::cerr << "coordinates" << pt << "=" << coords << std::endl;
      tribool result=true;
      for(dimension_type i=0; i!=this->dimension()+1; ++i) {
        result = result && (coords(i)>=0);
        //if(!result) { return result; }
      }
      return result;
    }

    template<class R>
    LinearAlgebra::Vector<typename Numeric::traits<R>::arithmetic_type>
    Simplex<R>::coordinates(const Point<R>& pt) const
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      LinearAlgebra::Vector<F> v(pt.dimension()+1);
      for(dimension_type i=0; i!=pt.dimension(); ++i) { 
        v(i)=pt[i];
      }
      v(pt.dimension())=1;
      LinearAlgebra::Matrix<F> Ginv=this->generators().inverse();
      return Ginv*v;
    }
    
    template<class R>
    std::ostream&
    Simplex<R>::write(std::ostream& os) const
    {
      const Simplex<R>& s=*this;
      if(s.dimension() > 0) {
        os << "Simplex( vertices=" << s.vertices() << " )";
      }
      return os;
    }
    
    template<class R>
    std::istream& 
    Simplex<R>::read(std::istream& is)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
      
  }
}
