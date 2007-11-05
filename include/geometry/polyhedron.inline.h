/***************************************************************************
 *            polyhedron.inline.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
  namespace Geometry {
    template<class X>
    class PolyhedronConstraintsIterator
      : public boost::iterator_facade<PolyhedronConstraintsIterator<X>,
                                      PolyhedralConstraint<X>,
                                      boost::forward_traversal_tag,
                                      const PolyhedralConstraint<X>&,
                                      const PolyhedralConstraint<X>*
                                      >
    {
     public:
      PolyhedronConstraintsIterator(const Polyhedron<X>& ply, const size_type& n)
        : _c(ply.dimension(),ply.constraints().begin()+n*(ply.dimension()+1)) { }
      bool equal(const PolyhedronConstraintsIterator<X>& other) const { 
        return this->_c._a==other._c._a; }
      const PolyhedralConstraint<X>& dereference() const { return _c; }
      void increment() { _c._a+=_c._d+1u; ; }
     private:
      PolyhedralConstraint<X> _c;
    };
  
  } // namespace Geometry
} // namespace Ariadne 




namespace Ariadne {

template<class X> template<class XX> inline
Geometry::Polyhedron<X>::Polyhedron(const Rectangle<XX>& r)
  : _dimension(r.dimension()), 
    _number_of_constraints(r.dimension()*2u),
    _data((r.dimension()+1u)*r.dimension()*2u,static_cast<X>(0))
{
  dimension_type d=r.dimension();
  LinearAlgebra::MatrixSlice<X> constraints=this->_constraints();
  for(size_type i=0; i!=d; ++i) {
    constraints(i,i)=static_cast<X>(1);
    constraints(i,d)=-r.lower_bound(i);
    constraints(i+d,i)=static_cast<X>(-1); 
    constraints(i+d,d)=r.upper_bound(i);
  }
}

template<class X> template<class XX> inline 
Geometry::Polyhedron<X>::Polyhedron(const Polytope<XX>& p)
{ 
  (*this)=polyhedron(p);
}

template<class X> template<class XX> inline 
Geometry::Polyhedron<X>::Polyhedron(const Polyhedron<XX>& p)
  : _dimension(p.dimension()), 
    _number_of_constraints(p.number_of_constraints()), 
    _data(p.data())
{ 
}

template<class X> template<class XX> inline
Geometry::Polyhedron<X>&
Geometry::Polyhedron<X>::operator=(const Polyhedron<XX>& plhd)
  
{
  if(this!=&plhd) {
    this->_dimension=plhd.dimension();
    this->_number_of_constraints=plhd.number_of_constraints();
    this->_data=plhd.data();
  }
  return *this;
}

template<class X> template<class XX> inline
tribool 
Geometry::Polyhedron<X>::contains(const Point<XX>& pt) const
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,pt,"tribool Polyhedron::contains(Point pt)");
  tribool result=true;
  for(constraints_const_iterator i=this->constraints_begin(); 
      i!=this->constraints_end(); ++i) {
    result=result && i->satisfied_by(pt); 
    if(!result) { return result; }
  }
  return result;
}


template<class X> inline 
array<X>& 
Geometry::Polyhedron<X>::data() 
{ 
  return this->_data; 
}

template<class X> inline 
const array<X>& 
Geometry::Polyhedron<X>::data() const
{ 
  return this->_data; 
}


template<class X> inline 
X* 
Geometry::Polyhedron<X>::begin() 
{
  return this->_data.begin(); 
}


template<class X> inline 
const X* 
Geometry::Polyhedron<X>::begin() const 
{
  return this->_data.begin(); 
}



template<class X> inline 
const LinearAlgebra::MatrixSlice<X> 
Geometry::Polyhedron<X>::constraints() const 
{ 
  return LinearAlgebra::MatrixSlice<X>(this->_number_of_constraints,
                                       this->_dimension+1u,
                                       const_cast<X*>(this->begin()),
                                       this->_dimension+1u,
                                       1u);
}


template<class X> inline 
LinearAlgebra::Matrix<X> 
Geometry::Polyhedron<X>::A() const 
{ 
  return -LinearAlgebra::MatrixSlice<X>(this->number_of_constraints(),
                                        this->dimension(),
                                        const_cast<X*>(this->begin()),
                                        this->dimension()+1,1u);
}

template<class X> inline 
LinearAlgebra::Vector<X> 
Geometry::Polyhedron<X>::b() const 
{ 
  return LinearAlgebra::VectorSlice<X>(this->number_of_constraints(),
                                       const_cast<X*>(this->begin()+this->dimension()),
                                       this->dimension()+1u);
}


template<class X> inline 
size_type 
Geometry::Polyhedron<X>::number_of_constraints() const 
{ 
  return this->_number_of_constraints;
}





template<class X> inline
Geometry::PolyhedralConstraint<X>::PolyhedralConstraint(const dimension_type d, const X* a)
  : _d(d), _a(a) 
{ 
}


template<class X> inline
dimension_type
Geometry::PolyhedralConstraint<X>::dimension() const
{ 
  return this->_d;
}


template<class X1> template<class X2> inline 
tribool 
Geometry::PolyhedralConstraint<X1>::satisfied_by(const Point<X2>& pt) const
{
  typedef typename Numeric::traits<X1,X2>::arithmetic_type F;
  F prod=0;
  for(dimension_type i=0; i!=pt.dimension(); ++i) {
    prod+=this->_a[i]*pt[i]; 
  }
  prod+=this->_a[pt.dimension()]; 
  return prod>=0;
}











template<class X> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const PolyhedralConstraint<X>& c) 
{
  return c.write(os); 
}


template<class X> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const Polyhedron<X>& p) 
{
  return p.write(os); 
}


template<class X> inline
std::istream& 
Geometry::operator>>(std::istream& os, Polyhedron<X>& p) 
{
  return p.read(os); 
}




} // namespace Ariadne
