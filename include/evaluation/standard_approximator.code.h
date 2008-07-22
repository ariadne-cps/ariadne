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
#include "geometry/partition_tree_set.h"
#include "evaluation/standard_approximator.h"

namespace Ariadne {






template<class ES>
ES
StandardApproximator<ES>::enclosure_set(const Box<R>& r) const
{
  return ES(r);
}

template<class ES>
typename ES::real_type
StandardApproximator<ES>::radius(const ES& es) const
{
  return Ariadne::bounding_box(es).radius();
}

template<class ES>
Box<typename ES::real_type>
StandardApproximator<ES>::bounding_box(const ES& es) const
{
  return Ariadne::bounding_box(es);
}

template<class ES>
tribool
StandardApproximator<ES>::disjoint(const ES& es, const Box<R>& bx) const
{
  return Ariadne::disjoint(es,bx);
}

template<class ES>
tribool
StandardApproximator<ES>::superset(const ES& es, const Box<R>& bx) const
{
  return Ariadne::subset(bx,es);
}

template<class ES>
void
StandardApproximator<ES>::adjoin_over_approximation(BoxListSet<R>& bxls, const ES& es) const
{
  bxls.adjoin(Ariadne::bounding_box(es));
}

template<class ES>
void
StandardApproximator<ES>::adjoin_outer_approximation(GridCellListSet<R>& gcls, const ES& es) const
{
  gcls.adjoin(Ariadne::outer_approximation(es,gcls.grid()));
}

template<class ES>
void
StandardApproximator<ES>::adjoin_outer_approximation(GridMaskSet<R>& gms, const ES& es) const
{
  gms.adjoin(Ariadne::outer_approximation(es,gms.grid()));
  //gms.adjoin_outer_approximation(es);
}



template<class ES>
std::ostream&
StandardApproximator<ES>::write(std::ostream& os) const
{
  return os << "StandardApproximator( )\n";
}


} // namespace Ariadne
