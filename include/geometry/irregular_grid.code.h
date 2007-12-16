/***************************************************************************
 *            irregular_grid.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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


#include "base/stlio.h"

#include "rectangle.h"
#include "list_set.h"

namespace Ariadne {



template<class R>
Geometry::IrregularGrid<R>::IrregularGrid(const ListSet< Rectangle<R> >& ls)
  : _subdivision_coordinates(ls.dimension()),
    _centre_positions(ls.dimension(),0u)
{
  dimension_type d=ls.dimension();
  R one=static_cast<R>(1);
  
  std::vector<R>* rectangle_bounds=new std::vector<R>[d];
  
  for(typename ListSet< Rectangle<R> >::const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
    for(dimension_type i=0; i!=d; ++i) {
      rectangle_bounds[i].push_back(riter->lower_bound(i));
      rectangle_bounds[i].push_back(riter->upper_bound(i));
    }
  }
  
  for(dimension_type i=0; i!=d; ++i) {
    std::vector<R>& bounds=rectangle_bounds[i];
    bounds.push_back(static_cast<R>(0));
    std::sort(bounds.begin(),bounds.end());
    typename std::vector<R>::const_iterator bounds_end=std::unique(bounds.begin(),bounds.end());
    bounds.resize(bounds_end-bounds.begin());
    
    array<R>& coordinates=this->_subdivision_coordinates[i];
    coordinates.resize(bounds.size());
    std::copy(bounds.begin(),bounds.end(),coordinates.begin());
    
    this->_centre_positions[i] = std::find(coordinates.begin(),coordinates.end(),static_cast<R>(0))-coordinates.begin();
  }
  delete[] rectangle_bounds;
  this->create();
}


template<class R>
Geometry::IrregularGrid<R>::IrregularGrid(const IrregularGrid<R>& g1, IrregularGrid<R>& g2)
  : _subdivision_coordinates(g1.dimension()),
    _centre_positions(g1.dimension(),0)
{
  for(dimension_type d=0; d!=this->dimension(); ++d) {
    array<R>& sc(this->_subdivision_coordinates[d]);
    const array<R>& sc1(g1._subdivision_coordinates[d]);
    const array<R>& sc2(g2._subdivision_coordinates[d]);
    sc.resize(sc1.size()+sc2.size());
    std::merge(sc1.begin(),sc1.end(),sc2.begin(),sc2.end(),sc.begin());
    index_type n=std::unique(sc.begin(),sc.end())-sc.begin();
    sc.resize(n);
    R x=sc1[g1._centre_positions[d]];
    this->_centre_positions[d]=std::find(sc.begin(),sc.end(),x)-sc.begin();
  }
  
  this->create();
}



template<class R>
dimension_type
Geometry::IrregularGrid<R>::dimension() const
{
  return this->_subdivision_coordinates.size();
}



template<class R>
R
Geometry::IrregularGrid<R>::subdivision_coordinate(dimension_type d, index_type n) const 
{
  const array<R>& subdiv=this->_subdivision_coordinates[d];
  if( (n < -int(_centre_positions[d])) || 
      n > int(subdiv.size()-_centre_positions[d]) ) {
    throw std::runtime_error("Index lies outside grid block");
  }
  return subdiv[n+_centre_positions[d]];
}



template<class R> 
index_type 
Geometry::IrregularGrid<R>::subdivision_index(dimension_type d, const real_type& x) const 
{
  const array<R>& subdiv=this->_subdivision_coordinates[d];
  if(x < subdiv.front() || x > subdiv.back()) {
    std::cerr << __FUNCTION__ << *this << " " << d << "  " << x << std::endl;
    throw std::runtime_error("Value does not lie in extent of irregular grid");
  }
  typename array<R>::const_iterator pos=std::upper_bound(subdiv.begin(),subdiv.end(),x);
  if(*pos!=x) {
    throw std::runtime_error("Value is not a coordinate of the irregular grid");
  }
  return (pos-subdiv.begin())+this->_centre_positions[d];
}


template<class R> 
index_type 
Geometry::IrregularGrid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const 
{
  const array<R>& subdiv=this->_subdivision_coordinates[d];
  if(x < subdiv.front() || x > subdiv.back()) {
    std::cerr << __FUNCTION__ << *this << " " << d << "  " << x << std::endl;
    throw std::runtime_error("Value does not lie in extent of irregular grid");
  }
  typename array<R>::const_iterator pos = std::upper_bound(subdiv.begin(),subdiv.end(),x);
  return (pos - subdiv.begin()) + this->_centre_positions[d] - 1;
}


template<class R> 
index_type 
Geometry::IrregularGrid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
{
  const array<R>& subdiv=this->_subdivision_coordinates[d];
  if(x < subdiv.front() || x > subdiv.back()) {
    std::cerr << __FUNCTION__ << *this << " " << d << "  " << x << std::endl;
    throw std::runtime_error("Value does not lie in extent of irregular grid");
  }
  typename array<R>::const_iterator pos = std::lower_bound(subdiv.begin(),subdiv.end(),x);
  return (pos - subdiv.begin()) + this->_centre_positions[d] + 1;
}


template<class R> inline 
Combinatoric::LatticeBlock
Geometry::IrregularGrid<R>::index_block(const Rectangle<R>& r) const 
{ 
  Combinatoric::LatticeBlock result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result.set_lower_bound(i,this->subdivision_lower_index(i,r.lower_bound(i)));
    result.set_upper_bound(i,this->subdivision_upper_index(i,r.upper_bound(i)));
  }
  return result;
}


template<class R> inline 
Geometry::Rectangle<R>
Geometry::IrregularGrid<R>::rectangle(const Combinatoric::LatticeBlock& lb) const 
{ 
  Rectangle<R> result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result.set_lower_bound(i,this->subdivision_coordinate(i,lb.lower_bound(i)));
    result.set_upper_bound(i,this->subdivision_coordinate(i,lb.upper_bound(i)));
  }
  return result;
}


template<class R> inline 
Combinatoric::LatticeBlock
Geometry::IrregularGrid<R>::block() const 
{ 
  Combinatoric::LatticeBlock result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result.set_lower_bound(i,-this->_centre_positions[i]);
    result.set_upper_bound(i,this->_subdivision_coordinates[i].size()-this->_centre_positions[i]-1);
  }
  return result;
}


template<class R> inline 
SizeArray 
Geometry::IrregularGrid<R>::sizes() const 
{ 
  return this->block().sizes(); 
}

template<class R> inline 
size_type 
Geometry::IrregularGrid<R>::capacity() const 
{ 
  return this->block().size(); 
}

template<class R> inline 
size_type 
Geometry::IrregularGrid<R>::size(dimension_type i) const 
{ 
  return this->_subdivision_coordinates[i].size()-1; 
}


template<class R>
void
Geometry::IrregularGrid<R>::create() 
{
}


template<class R>
std::ostream&
Geometry::IrregularGrid<R>::write(std::ostream& os) const
{
  os << "IrregularGrid( subdivision_coordinates=" << this->_subdivision_coordinates
     << ", centre_positions=" << this->_centre_positions << " )";
  return os;
}


} // namespace Ariadne



