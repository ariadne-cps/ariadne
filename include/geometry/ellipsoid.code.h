/***************************************************************************
 *            ellispoid.code.h
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
 
#include <stdexcept>

#include "numeric/rational.h"
#include "ellipsoid.h"

#include "geometry/exceptions.h"
#include "geometry/sphere.h"


namespace Ariadne {

 
      
template<class R>
Ellipsoid<R>::Ellipsoid(size_type n)
  : _centre(n), _bilinear_form(Matrix<R>::identity(n))
{
}

template<class R>
Ellipsoid<R>::Ellipsoid(const Point<R>& c, const Matrix<R>& A)
  : _centre(c), _bilinear_form(A)
{
  if(c.dimension()!=A.number_of_rows() && A.number_of_rows()!=A.number_of_columns()) {
    ARIADNE_THROW(IncompatibleDimensions,"Ellipsoid::Ellipsoid(Point c, Matrix A)","c="<<c<<", A="<<A);
  }
}

template<class R>
Ellipsoid<R>::Ellipsoid(const std::string& s)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Ellipsoid<R>::Ellipsoid(const Sphere<R>& s)
  : _centre(s.centre()), _bilinear_form(s.dimension(),s.dimension())
{ 
  for(size_type i=0; i!=this->dimension(); ++i) {
    this->_bilinear_form(i,i) = div_up(R(1),mul_down(s.radius(),s.radius()));
  }
}


template<class R>
tribool 
Ellipsoid<R>::contains(const Point<R>& point) const 
{
  Vector<Rational> p=point.position_vector();
  Vector<Rational> c=this->centre().position_vector();
  Vector<Rational> d=p-c;
  Matrix<Rational> A=this->bilinear_form();
  Rational r=inner_product(d,Vector<Rational>(A*d));
  if(r<1) { return true; }
  if(r>1) { return false; }
  return indeterminate;
}


template<class R>
void
Ellipsoid<R>::_instantiate_geometry_operators()
{
  R sf=1.0;
  Ellipsoid<R>* e=0;
  *e=scale(*e,sf);
}



template<class R>
Ellipsoid<R> 
scale(const Ellipsoid<R>& s, const R& scale_factor) 
{
  const Point<R>& centre=s.centre();
  const Matrix<R>& bilinear_form=s.bilinear_form();
  
  Point<R> new_centre(s.dimension());
  Matrix<R> new_bilinear_form(s.dimension(),s.dimension());
  
  for(size_type i=0; i!=s.dimension(); ++i) {
    new_centre[i]=mul_approx(scale_factor,centre[i]);
  }
  
  for(size_type i=0; i!=s.dimension(); ++i) {
    for(size_type j=0; j!=s.dimension(); ++j) {
      new_bilinear_form(i,j)=mul_up(scale_factor,bilinear_form(i,j));
    }
  }
  
  return Ellipsoid<R>(new_centre, new_bilinear_form);
}





template<class R>
std::ostream&
Ellipsoid<R>::write(std::ostream& os) const
{
  if(this->empty()) {
    os << "Empty";
  }
  else if(this->dimension() > 0) {
    os << "Ellipsoid( centre=" << this->centre() << ", axes=" << this->bilinear_form() << " )";
  }
  return os;
}

template<class R>
std::istream& 
Ellipsoid<R>::read(std::istream& is)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  
}

} // namespace Ariadne
