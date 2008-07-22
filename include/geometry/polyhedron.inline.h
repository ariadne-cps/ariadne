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

#include "geometry/halfspace.h"

namespace Ariadne {  
  

    template<class X>
    class PolyhedronConstraintsIterator
      : public boost::iterator_facade<PolyhedronConstraintsIterator<X>,
                                      Halfspace<X>,
                                      boost::forward_traversal_tag,
                                      const Halfspace<X>&,
                                      const Halfspace<X>*
                                      >
    {
     public:
      PolyhedronConstraintsIterator(const Polyhedron<X>& ply, const size_type& i)
        : _p(&ply), _i(i), _c(ply.dimension()) { }
      bool equal(const PolyhedronConstraintsIterator<X>& other) const { 
        return this->_i==other._i && this->_p ==other._p; }
      void increment() { 
        ++this->_i; }
      const Halfspace<X>& dereference() const { 
        this->_c=this->_p->constraint(this->_i); return this->_c; }
     private:
      const Polyhedron<X>* _p; size_type _i; mutable Halfspace<X> _c;
    };
  



template<class X> template<class XX>
Polyhedron<X>::Polyhedron(const Box<XX>& r)
  : _dimension(r.dimension()), 
    _number_of_constraints(r.dimension()*2u),
    _data((r.dimension()+1u)*r.dimension()*2u,static_cast<X>(0))
{
  dimension_type d=r.dimension();
  MatrixSlice<X> constraints=this->_constraints();
  for(size_type i=0; i!=d; ++i) {
    constraints(i,i)=static_cast<X>(1);
    constraints(i,d)=-r.lower_bound(i);
    constraints(i+d,i)=static_cast<X>(-1); 
    constraints(i+d,d)=r.upper_bound(i);
  }
}


template<class X> template<class XX>  
Polyhedron<X>::Polyhedron(const Polytope<XX>& p)
{ 
  (*this)=polyhedron(p);
}

template<class X> template<class XX>  
Polyhedron<X>::Polyhedron(const Polyhedron<XX>& p)
  : _dimension(p.dimension()), 
    _number_of_constraints(p.number_of_constraints()), 
    _data(p.data())
{ 
}

template<class X> template<class XX> inline
Polyhedron<X>&
Polyhedron<X>::operator=(const Polyhedron<XX>& plhd)
  
{
  if(this!=(void*)&plhd) {
    this->_dimension=plhd.dimension();
    this->_number_of_constraints=plhd.number_of_constraints();
    this->_data=plhd.data();
  }
  return *this;
}

template<class X> template<class XX> inline
tribool 
Polyhedron<X>::contains(const Point<XX>& pt) const
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,pt,"tribool Polyhedron::contains(Point pt)");
  tribool result=true;
  for(constraints_const_iterator i=this->constraints_begin(); 
      i!=this->constraints_end(); ++i) {
    result=result && satisfies(pt,*i); 
    if(!result) { return result; }
  }
  return result;
}


template<class X> inline 
array<X>& 
Polyhedron<X>::data() 
{ 
  return this->_data; 
}

template<class X> inline 
const array<X>& 
Polyhedron<X>::data() const
{ 
  return this->_data; 
}





template<class X> inline 
const MatrixSlice<X> 
Polyhedron<X>::constraints() const 
{ 
  return MatrixSlice<X>(this->_number_of_constraints,
                                       this->_dimension+1u,
                                       const_cast<X*>(this->data().begin()),
                                       this->_dimension+1u,
                                       1u);
}


template<class X> inline 
Matrix<X> 
Polyhedron<X>::A() const 
{ 
  return -MatrixSlice<X>(this->number_of_constraints(),
                                        this->dimension(),
                                        const_cast<X*>(this->data().begin()),
                                        this->dimension()+1,1u);
}

template<class X> inline 
Vector<X> 
Polyhedron<X>::b() const 
{ 
  return VectorSlice<X>(this->number_of_constraints(),
                                       const_cast<X*>(this->data().begin()+this->dimension()),
                                       this->dimension()+1u);
}


template<class X> inline 
size_type 
Polyhedron<X>::number_of_constraints() const 
{ 
  return this->_number_of_constraints;
}


template<class X> inline
Box<typename Polyhedron<X>::real_type> 
bounding_box(const Polyhedron<X>& plhd)
{
  return plhd.bounding_box();
}



template<class X> inline
std::ostream& 
operator<<(std::ostream& os, const Polyhedron<X>& p) 
{
  return p.write(os); 
}


template<class X> inline
std::istream& 
operator>>(std::istream& os, Polyhedron<X>& p) 
{
  return p.read(os); 
}




} // namespace Ariadne
