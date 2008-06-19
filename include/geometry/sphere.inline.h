/***************************************************************************
 *            sphere.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
     
namespace Ariadne {


template<class R> inline
Sphere<R>::Sphere(const Sphere<R>& original)
  : _centre(original._centre), _radius(original._radius)
{ 
}

template<class R> inline
Sphere<R>& 
Sphere<R>::operator=(const Sphere<R>& original) 
{
  if(this != &original) {
    this->_centre = original._centre;
    this->_radius = original._radius;
  }
  return *this;
}


template<class R> inline
bool 
Sphere<R>::operator==(const Sphere<R>& other) const
{
  return this->_centre==other._centre && this->_radius==other._radius;
}

template<class R> inline
bool 
Sphere<R>::operator!=(const Sphere<R>& other) const 
{
  return !(*this == other);
}


template<class R> inline
const Point<R>& 
Sphere<R>::centre() const 
{
  return this->_centre;
}

template<class R> inline
const R& 
Sphere<R>::radius() const 
{
  return this->_radius;
}

template<class R> inline
size_type 
Sphere<R>::dimension() const
{
  return this->_centre.dimension();
}

template<class R> inline
bool 
Sphere<R>::empty() const 
{
  return this->_radius<R(0);
}


template<class R> inline
bool
Sphere<R>::empty_interior() const 
{
  return this->_radius <= R(0);
}


template<class R> inline
tribool
Sphere<R>::contains(const Point<R>& point) const 
{
  if(euclidean_distance_square_down(point,this->_centre) < mul_approx(this->_radius , this->_radius)) {
    return true; 
  } else if(euclidean_distance_square_up(point,this->_centre) > mul_approx(this->_radius , this->_radius)) {
    return false;
  } else {
    return indeterminate;
  }
}


template<class R> inline
R
euclidean_distance_square_down(const Point<R>& pt1, const Point<R>& pt2) 
{
  R result=R(0);
  R tmp;
  for(dimension_type i=0; i!=pt1.dimension(); ++i) {
    tmp = (pt1[i]>=pt2[i]) ? sub_down(pt1[i],pt2[i]) : sub_down(pt2[i],pt1[i]);
    result=add_down(result,mul_down(tmp,tmp));
  }
  return result;
}

template<class R> inline
R
euclidean_distance_square_up(const Point<R>& pt1, const Point<R>& pt2) 
{
  R result=R(0);
  R tmp;
  for(dimension_type i=0; i!=pt1.dimension(); ++i) {
    tmp = (pt1[i]>=pt2[i]) ? sub_up(pt1[i],pt2[i]) : sub_up(pt2[i],pt1[i]);
    result=add_up(result,mul_up(tmp,tmp));
  }
  return result;
}

template<class R> inline 
tribool 
disjoint(const Sphere<R>& A, const Sphere<R>& B) 
{
  if(euclidean_distance_down(A.centre(),B.centre()) > 
     pow_up(add_up(A.radius(),B.radius()),2)) {
    return true; 
  } else if(euclidean_distance_up(A.centre(),B.centre()) < 
            pow_down(add_down(A.radius(),B.radius()),2)) {
    return false;
  } else { 
    return indeterminate;
  }
}

template<class R> inline 
tribool 
disjoint(const Sphere<R>& A, const Box<R>& B) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R> inline 
tribool 
disjoint(const Box<R>& A, const Sphere<R>& B) 
{
  return disjoint(B,A);
}



template<class R> inline 
tribool 
subset(const Sphere<R>& A, const Sphere<R>& B) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  //return A.radius()<=B.radius && euclidean_distance_square(A.centre()-B.centre()) <= square(B.centre()-A.centre());
}

template<class R> inline 
tribool 
subset(const Sphere<R>& A, const Box<R>& B) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  //return subset(A.bounding_box(),B);
}

template<class R> inline 
tribool 
subset(const Box<R>& A, const Sphere<R>& B) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  //array< Point<R> > vertices=A.vertices();
  //for(class Box<R>::vertex_iterator vertex_iter=vertices.begin(); vertex_iter!=vertices.end(); ++vertex_iter) {
  //  if(! B.contains(*vertex_iter) ) {
  //    return false;
  //  }
  //}
  return true;
}


template<class R> inline
Sphere<R> 
scale(const Sphere<R>& s, const R& scale_factor) {
  
  const Point<R>& centre=s.centre();
  Point<R> new_centre(s.dimension());
  
  for(size_type i=0; i!=s.dimension(); ++i) {
    new_centre[i]=scale_factor*centre[i];
  }
  
  return Sphere<R>(new_centre, scale_factor*s.radius());
}

template<class R> inline
std::ostream& 
operator<<(std::ostream& os, const Sphere<R>& s) {
  return s.write(os);
}


template<class R> inline
std::istream& 
operator>>(std::istream& is, Sphere<R>& s) {
  return s.read(is);
}



} // namespace Ariadne
