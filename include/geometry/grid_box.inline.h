/***************************************************************************
 *            grid_box.inline.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 * 
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
 
#include "geometry/rectangle.h"

namespace Ariadne {

template<class R> inline 
const Geometry::Grid<R>& 
Geometry::GridBox<R>::grid() const 
{ 
  return *this->_grid_ptr; 
}

template<class R> inline
dimension_type 
Geometry::GridBox<R>::dimension() const 
{
  return this->_coordinates.size()/2; 
}



template<class R> inline
R
Geometry::GridBox<R>::lower_bound(dimension_type i) const 
{
  return _grid_ptr->subdivision_coordinate(i,_coordinates[2*i]);
}


template<class R> inline
R
Geometry::GridBox<R>::upper_bound(dimension_type i) const 
{
  return _grid_ptr->subdivision_coordinate(i,_coordinates[2*i+1]);
}

template<class R> inline
Geometry::GridBox<R>
Geometry::GridBox<R>::subdivide(dimension_type i, bool lr) const 
{
  GridBox<R> result(*this);
  mc=(this->_coordinates[2*i+1]-this->_coordinates[2*i])/2
  if(lr==left) {
    result[2*i+1]=mc;
  } else {
    result[2*i]=mc;
  }
  return result;
}

template<class R> inline
Geometry::GridBox<R>
Geometry::GridBox<R>::subdivide(bool lr) const 
{
  dyadic_type mr=0;
  dimension_type mi=0;
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    dyadic_type r=this->_coordinates[2*i+1]-this->_coordinates[2*i];
    if(r>mr) {
      mr=r;
      mi=i;
    }
  }
  this->subdivide(mi,lr);
}


template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const GridBox<R>& gc) {
  return gc.write(os);
}
    

}
