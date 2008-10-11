/***************************************************************************
 *            grid.inline.h
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
 

namespace Ariadne {


template<class R> inline
const Geometry::Grid<R>& 
Geometry::FiniteGrid<R>::grid() const
{ 
  return this->_grid; 
}

template<class R> inline
const Combinatoric::LatticeBlock& 
Geometry::FiniteGrid<R>::lattice_block() const 
{
  return this->_lattice_block; 
}


template<class R> inline
R 
Geometry::FiniteGrid<R>::subdivision_coordinate(dimension_type d, integer_type n) const
{
  return this->_grid.subdivision_coordinate(d,n);
}

template<class R> inline
R 
Geometry::FiniteGrid<R>::subdivision_coordinate(dimension_type d, long_integer_type n) const
{
  return this->_grid.subdivision_coordinate(d,n);
}

template<class R> inline
typename Geometry::FiniteGrid<R>::integer_type 
Geometry::FiniteGrid<R>::subdivision_index(dimension_type d, const real_type& x) const 
{
  return this->_grid.subdivision_index(d,x);
}

template<class R> inline
typename Geometry::FiniteGrid<R>::integer_type 
Geometry::FiniteGrid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const 
{
  return this->_grid.subdivision_lower_index(d,x);
}

template<class R> inline
typename Geometry::FiniteGrid<R>::integer_type 
Geometry::FiniteGrid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
{
  return this->_grid.subdivision_upper_index(d,x);
}





template<class R> inline
std::istream& 
Geometry::operator>>(std::istream& is, Grid<R>& g)
{
  return g.read(is);
}

template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const Grid<R>& g)
{
  return g.write(os);
}

template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const FiniteGrid<R>& fg)
{
  return fg.write(os);
}

template<class R> template<class A> 
void 
Geometry::Grid<R>::serialize(A& archive, const unsigned int version) {
  archive & this->_data->_origin;
  archive & this->_data->_lengths;
}


} // namespace Ariadne
