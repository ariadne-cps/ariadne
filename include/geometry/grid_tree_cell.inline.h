/***************************************************************************
 *            grid_tree_cell.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


namespace Ariadne {


template<class R> inline
Geometry::GridTreeCell<R>::GridTreeCell(const Grid<R>& g, 
                                        const Combinatoric::BinaryWord& nw)
  : _grid(r), _word(bw)
{
}


template<class R> inline
const Geometry::Grid<R>& 
Geometry::GridTreeCell<R>::grid() const 
{
  return *this->_grid
}



template<class R> inline
dimension_type 
Geometry::GridTreeCell<R>::dimension() const 
{
  return this->_grid.dimension(); 
}

template<class R>
R
Geometry::GridTreeCell<R>::lower_bound(dimension_type i) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
R
Geometry::GridTreeCell<R>::upper_bound(dimension_type i) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  //  return add_approx(_unit_box.lower_bound(i),
  //                    mul_approx(_subdivision_cell.upper_bound(i),
  //                               sub_approx(_unit_box.upper_bound(i),_unit_box.lower_bound(i))));
}

template<class R> inline
GridTreeCell<R>
Geometry::PartitionTreeCell<R>::subdivide(bool lr) const 
{
  GridTreeCell result(*this);
  result._word.push_back(lr);
  return result;
}



template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const GridTreeCell<R>& gtc)
{ 
  return gtc.write(os);
}


} // namespace Ariadne

