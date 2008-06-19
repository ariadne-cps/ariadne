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

#include "base/stlio.h"

namespace Ariadne {


template<class X>
Simplex<X>::Simplex()
  : Polytope<X>(Matrix<X>(0,1))
{
}

template<class X>
Simplex<X>::Simplex(const Matrix<X>& A)
  : Polytope<X>(A)
{
  ARIADNE_CHECK_DIMENSION(*this,A.number_of_columns()-1,"Simplex::Simplex(Matrix A)");
}

template<class X>
Simplex<X>::Simplex(const PointList<X>& v)
  : Polytope<X>(v)
{
  ARIADNE_CHECK_DIMENSION(*this,v.size()-1,"Simplex::Simplex(PointList v)");
}

template<class X>
tribool
Simplex<X>::contains(const Point<X>& pt) const
{
  typedef typename traits<X>::arithmetic_type F;
  Vector<F> coords=this->coordinates(pt);
  //std::cerr << "coordinates" << pt << "=" << coords << std::endl;
  tribool result=true;
  for(dimension_type i=0; i!=this->dimension()+1; ++i) {
    result = result && (coords(i)>=0);
    //if(!result) { return result; }
  }
  return result;
}

template<class X>
Vector<typename traits<X>::arithmetic_type>
Simplex<X>::coordinates(const Point<X>& pt) const
{
  typedef typename traits<X>::arithmetic_type F;
  Vector<F> v(pt.dimension()+1);
  for(dimension_type i=0; i!=pt.dimension(); ++i) { 
    v(i)=pt[i];
  }
  v(pt.dimension())=1;
  Matrix<F> Ginv=inverse(this->generators());
  return Ginv*v;
}

template<class X>
std::ostream&
Simplex<X>::write(std::ostream& os) const
{
  const Simplex<X>& s=*this;
  if(s.dimension() > 0) {
    os << "Simplex( vertices=" << s.vertices() << " )";
  }
  return os;
}

template<class X>
std::istream& 
Simplex<X>::read(std::istream& is)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


} // namespace Ariadne
