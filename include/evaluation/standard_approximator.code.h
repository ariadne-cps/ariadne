/***************************************************************************
 *            standard_approximator.code.h
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
#include "geometry/grid_cell.h"
#include "geometry/grid_block.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "geometry/grid_approximation.h"
#include "evaluation/standard_approximator.h"

namespace Ariadne {



template<class BS>
Evaluation::StandardApproximator<BS>::StandardApproximator()
{
}

template<class BS>
Evaluation::StandardApproximator<BS>::StandardApproximator(const StandardApproximator<BS>& approx)
{
}


template<class BS>
Evaluation::StandardApproximator<BS>* 
Evaluation::StandardApproximator<BS>::clone() const
{
  return new StandardApproximator<BS>(*this);
}


template<class BS>
BS
Evaluation::StandardApproximator<BS>::over_approximation(const Geometry::Box<R>& r) const
{
  return BS(r);
}

template<class BS>
Geometry::Box<typename BS::real_type>
Evaluation::StandardApproximator<BS>::bounding_box(const BS& bs) const
{
  return Geometry::bounding_box(bs);
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::StandardApproximator<BS>::outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const
{
  return Geometry::outer_approximation(bs,g);
}

template<class BS>
std::pair<BS,BS>
Evaluation::StandardApproximator<BS>::subdivide(const BS& bs) const
{
  return Geometry::subdivide(bs);
  //Geometry::ListSet<BS> ls=Geometry::subdivide(bs);
  //assert(ls.size()==2);
  //return std::make_pair(ls[0],ls[1]);
}



} // namespace Ariadne
