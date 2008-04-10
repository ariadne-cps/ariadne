/***************************************************************************
 *            fast_approximator.code.h
 *
 *  Copyright  2007  Pieter Collins
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
 
#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/zonotope.h"
#include "geometry/grid_cell.h"
#include "geometry/grid_block.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "geometry/grid_approximation.h"
#include "evaluation/fast_approximator.h"

namespace Ariadne {


template<class BS>
Evaluation::FastApproximator<BS>::FastApproximator()
{
}

template<class BS>
Evaluation::FastApproximator<BS>::FastApproximator(const FastApproximator<BS>& approx)
{
}

template<class BS>
Evaluation::FastApproximator<BS>* 
Evaluation::FastApproximator<BS>::clone() const
{
  return new FastApproximator<BS>(*this);
}

template<class BS>
BS
Evaluation::FastApproximator<BS>::enclosure_set(const Geometry::Box<R>& r) const
{
  return BS(r);
}

template<class BS>
typename BS::real_type
Evaluation::FastApproximator<BS>::radius(const BS& bs) const
{
  return Geometry::bounding_box(bs).radius();
}

template<class BS>
Geometry::Box<typename BS::real_type>
Evaluation::FastApproximator<BS>::bounding_box(const BS& bs) const
{
  return Geometry::bounding_box(bs);
}

template<class BS>
Geometry::BoxListSet<typename BS::real_type>
Evaluation::FastApproximator<BS>::lower_approximation(const BS& bs) const
{
  Geometry::BoxListSet<R> result;
  result.adjoin(bs.bounding_box());
  return result;
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::FastApproximator<BS>::inner_approximation(const BS& bs, const Geometry::Grid<R>& g) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::FastApproximator<BS>::outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const
{
  return Geometry::fuzzy_outer_approximation(bs,g);
}


} // namespace Ariadne
