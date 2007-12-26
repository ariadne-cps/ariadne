/***************************************************************************
 *            grid.code.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "grid.h"

#include <ostream>
#include <iomanip>
#include <algorithm>

#include "exceptions.h"
#include "base/stlio.h"

#include "numeric/arithmetic.h"

#include "geometry/box.h"
#include "geometry/list_set.h"

#include "geometry/grid_block.h" // for GridBlock<R>
#include "geometry/grid_approximation.h"
#include "output/logging.h"

namespace Ariadne {


template<class Arg1, class Arg2, class Op> inline 
std::ostream& operator<<(std::ostream& os, const Numeric::Expression< Numeric::Binary<Op,Arg1,Arg2> >& e) { 
  return os << e.op << '(' << e.arg1 <<',' << e.arg2 << ')';
}


template<class R> 
Geometry::Grid<R>::~Grid()
{
}



template<class R> 
Geometry::Grid<R>::Grid(const dimension_type& d, const R& l)
  : _origin(d,static_cast<R>(0)),
    _lengths(d,l)
{
  this->create();
}


template<class R> 
Geometry::Grid<R>::Grid(const LinearAlgebra::Vector<R>& v)
  : _origin(v.size(),static_cast<R>(0)),
    _lengths(v.data())
{
  this->create();
}


template<class R> 
Geometry::Grid<R>::Grid(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
  : _origin(pt.data()),
    _lengths(v.data())
{
  if(pt.dimension() != v.size()) {
    throw IncompatibleDimensions(__PRETTY_FUNCTION__);
  }
  this->create();
}


template<class R> 
Geometry::Grid<R>::Grid(const Box<R>& r, const Combinatoric::LatticeBlock& lb)
  : _origin(r.dimension()),
    _lengths(r.dimension())
{
  if(r.dimension() != lb.dimension()) {
    throw IncompatibleDimensions(__PRETTY_FUNCTION__);
  }
  for(dimension_type i=0; i!=r.dimension(); ++i) {
    const R& l=r.lower_bound(i);
    const R& u=r.upper_bound(i);
    const index_type& a=lb.lower_bound(i);
    const index_type& b=lb.upper_bound(i);
    R& o=this->_origin[i];
    R& s=this->_lengths[i];
    
    //std::cerr << "a=" << a << ", b=" << b << ";  l=" << l << ", u=" << u << ";  "  << std::flush;
    
    if(a==b) {
      throw std::runtime_error("Invalid grid block");
    }
    if(a>0) { // o < l
      s=div_up(sub_up(u,l),b-a);
      o=sub_down(l,mul_up(s,a));
      while(add_approx(o,mul_approx(s,b))<u) {
        s=next_up(s);
        o=sub_down(l,mul_up(s,a));
      }
    } else if(b<0) { // o > u
      s=div_up(sub_up(u,l),b-a);
      o=sub_up(u,mul_down(s,b));
      while(add_approx(o,mul_approx(s,a))>l) {
        s=next_up(s);
        o=sub_up(u,mul_down(s,b));
      }
    } else { // l <= o <= u
      s=div_up(sub_up(u,l),b-a);
      o=div_approx(sub_approx(mul_approx(l,b),mul_approx(u,a)),b-a);
      while(add_approx(o,mul_approx(s,a))>l || add_approx(o,mul_approx(s,b))<u ) {
        s=next_up(s);
        o=div_approx(sub_approx(mul_approx(l,b),mul_approx(u,a)),b-a);
      }
    }
    
    //std::cerr << "o=" << o << "; s=" << s << std::endl;
    
    // check  o+sa <= l <= u <= o+sb
    //std::cerr << "add_approx(o,mul_approx(s,a))=" << add_approx(o,mul_approx(s,a)) << std::endl;
    assert(add_approx(o,mul_approx(s,a)) <= l);
    //std::cerr << "add_approx(o,mul_approx(s,b))=" << add_approx(o,mul_approx(s,b)) << std::endl;
    assert(add_approx(o,mul_approx(s,b)) >= u);
  }
  
  this->create();
}







template<class R>
void
Geometry::Grid<R>::create() 
{
}



template<class R>
dimension_type
Geometry::Grid<R>::dimension() const
{
  return this->_lengths.size();
}



template<class R>
Geometry::Point<R>
Geometry::Grid<R>::origin() const
{
  return Point<R>(LinearAlgebra::Vector<R>(this->_origin));
}


template<class R>
LinearAlgebra::Vector<R>
Geometry::Grid<R>::lengths() const
{
  return LinearAlgebra::Vector<R>(this->_lengths);
}




template<class R>
R
Geometry::Grid<R>::subdivision_coordinate(dimension_type d, dyadic_type x) const 
{
  return add_approx(this->_origin[d],mul_approx(this->_lengths[d],x));
}

template<class R>
R
Geometry::Grid<R>::subdivision_coordinate(dimension_type d, index_type n) const 
{
  return add_approx(this->_origin[d],mul_approx(this->_lengths[d],n));
}


template<class R> 
index_type 
Geometry::Grid<R>::subdivision_index(dimension_type d, const real_type& x) const 
{
  using namespace Numeric;
  
  R half=0.5;
  index_type n=floor(add_approx(div_approx(sub_approx(x,this->_origin[d]),this->_lengths[d]),half));
  R sc=add_approx(this->_origin[d],mul_approx(this->_lengths[d],n));
  ARIADNE_LOG(9,std::setprecision(20) << std::boolalpha << "sc=" << sc << " x=" << x << " sc-x=" << Interval<R>(sc-x) << "\n")
    if(sc == x) { 
      return n; 
    } else {
      std::cerr << std::setprecision(20) << std::boolalpha
                << "sc=" << sc << " x=" << x << " sc-x=" << Interval<R>(sc-x) << "\n"
                << "sc==x=" << (sc==x) << " sc!=x=" << (sc!=x)
                << " sc<x=" << (sc<x) << " sc>x=" << (sc>x) << " sc<=x=" << (sc<=x) << " sc>=x=" << (sc>=x) << std::endl; 
      ARIADNE_THROW(InvalidGridPosition,std::setprecision(20)<<"Grid::subdivision_index(dimension_type d,real_type x)","d="<<d<<", x="<<x<<", this->origin[d]="<<this->_origin[d]<<", this->lengths[d]="<<this->_lengths[d]<<" (closest value is "<<sc<<")");
    }
}


template<class R> 
index_type 
Geometry::Grid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const 
{
  using namespace Numeric;
  
  index_type n=floor(div_down(sub_down(x,this->_origin[d]),this->_lengths[d]));
  if(x>=add_approx(this->_origin[d],mul_approx(this->_lengths[d],(n+1)))) {
    return n+1;
  } else {
    return n;
  }
}


template<class R> 
index_type 
Geometry::Grid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
{
  using namespace Numeric;
  
  index_type n=ceil(div_up(sub_up(x,this->_origin[d]),this->_lengths[d]));
  if(x<=add_approx(this->_origin[d],mul_approx(this->_lengths[d],(n-1)))) {
    return n-1;
  } else {
    return n;
  }
}



template<class R> 
bool 
Geometry::Grid<R>::operator==(const Grid<R>& g) const
{
  if(this==&g) { 
    return true; 
  } else {
    return this->_origin==g._origin && this->_lengths==g._lengths;
  }
}


template<class R> 
bool 
Geometry::Grid<R>::operator!=(const Grid<R>& g) const
{
  return !(*this==g);
}


template<class R>
IndexArray 
Geometry::Grid<R>::index(const Point<R>& s) const
{
  IndexArray res(s.dimension());
  for(size_type i=0; i!=res.size(); ++i) {
    res[i]=subdivision_index(i,s[i]);
  }
  return res;
}


template<class R>
IndexArray  
Geometry::Grid<R>::lower_index(const Box<R>& r) const {
  IndexArray res(r.dimension());
  for(size_type i=0; i!=res.size(); ++i) {
    res[i]=subdivision_lower_index(i,r.lower_bound(i));
  }
  return res;
}


template<class R>
IndexArray  
Geometry::Grid<R>::upper_index(const Box<R>& r) const {
  IndexArray res(r.dimension());
  for(size_type i=0; i!=res.size(); ++i) {
    res[i]=subdivision_upper_index(i,r.upper_bound(i));
  }
  return res;
}


template<class R>
Combinatoric::LatticeBlock  
Geometry::Grid<R>::index_block(const Box<R>& r) const {
  Combinatoric::LatticeBlock res(r.dimension());
  for(size_type i=0; i!=res.dimension(); ++i) {
    res.set_lower_bound(i,this->subdivision_lower_index(i,r.lower_bound(i)));
    res.set_upper_bound(i,this->subdivision_upper_index(i,r.upper_bound(i)));
  }
  return res;
}


template<class R>
Geometry::Point<R> 
Geometry::Grid<R>::point(const IndexArray& a) const
{
  Point<R> res(a.size());
  for(size_type i=0; i!=res.dimension(); ++i) {
    res[i]=subdivision_coordinate(i,a[i]);
  }
  return res;
}


template<class R>
Geometry::Box<R> 
Geometry::Grid<R>::rectangle(const Combinatoric::LatticeBlock& lb) const
{
  Box<R> res(lb.dimension());
  for(size_type i=0; i!=res.dimension(); ++i) {
    res.set_lower_bound(i,this->subdivision_coordinate(i,lb.lower_bound(i)));
    res.set_upper_bound(i,this->subdivision_coordinate(i,lb.upper_bound(i)));
  }
  return res;
}


template<class R>
std::ostream&
Geometry::Grid<R>::write(std::ostream& os) const
{
  os << "Grid( ";
  if(this->_origin!=array<R>(this->dimension(),0)) {
    os << "origin=" << this->_origin << ", ";
  }
  os  << "lengths=" << this->_lengths << " )";
  return os;
}



template<class R>
std::istream&
Geometry::Grid<R>::read(std::istream& is) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}







template<class R>
Geometry::FiniteGrid<R>::~FiniteGrid() 
{
  if(this->_own_ptr) {
    delete this->_grid_ptr;
  }
}


template<class R> inline
Geometry::FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Combinatoric::LatticeBlock& b) 
  : _own_ptr(0), _grid_ptr(&g), _lattice_block(b)
{ 
}

template<class R> inline
Geometry::FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Box<R>& bb) 
  : _own_ptr(0), _grid_ptr(&g), _lattice_block(over_approximation(bb,g).lattice_set())
{ 
}

template<class R>
Geometry::FiniteGrid<R>::FiniteGrid(const Box<R>& bb, const size_type& s)
  : _own_ptr(1), _grid_ptr(0), _lattice_block(bb.dimension())
{
  dimension_type d=bb.dimension();
  array<R> subdivision_lengths(bb.dimension());
  for(dimension_type i=0; i!=bb.dimension(); ++i) {
    subdivision_lengths[i]=div_up(bb[i].width(),R(s));
  }
  this->_grid_ptr=new Grid<R>(LinearAlgebra::Vector<R>(d,subdivision_lengths.begin()));
  this->_lattice_block=this->grid().index_block(bb);
}


template<class R>
dimension_type 
Geometry::FiniteGrid<R>::dimension() const 
{
  return this->_grid_ptr->dimension();
}


template<class R>
Geometry::Box<R>
Geometry::FiniteGrid<R>::extent() const
{
  return this->_grid_ptr->rectangle(this->_lattice_block);
}


template<class R>
std::ostream&
Geometry::FiniteGrid<R>::write(std::ostream& os) const
{
  return os << "FiniteGrid(grid=" << this->grid() << ", lattice_block=" << this->lattice_block() << ")";
}



} // namespace Ariadne
