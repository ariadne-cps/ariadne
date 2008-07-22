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


template<class ES>
FastApproximator<ES>* 
FastApproximator<ES>::clone() const
{
  return new FastApproximator<ES>(*this);
}

template<class ES>
ES
FastApproximator<ES>::enclosure_set(const Box<R>& r) const
{
  return ES(r);
}

template<class ES>
typename ES::real_type
FastApproximator<ES>::radius(const ES& es) const
{
  return Ariadne::bounding_box(es).radius();
}

template<class ES>
tribool
FastApproximator<ES>::disjoint(const ES& es, const Box<R>& bx) const
{
  return Ariadne::disjoint(es,bx);
}

template<class ES>
tribool
FastApproximator<ES>::superset(const ES& es, const Box<R>& bx) const
{
  return Ariadne::superset(es,bx);
}

template<class ES>
Box<typename ES::real_type>
FastApproximator<ES>::bounding_box(const ES& es) const
{
  return Ariadne::bounding_box(es);
}

template<class ES>
BoxListSet<typename ES::real_type>
FastApproximator<ES>::over_approximation(const ES& es) const
{
  BoxListSet<R> result;
  result.adjoin(Ariadne::bounding_box(es));
  return result;
}

template<class ES>
void
FastApproximator<ES>::adjoin_outer_approximation(GridCellListSet<R>& gcls, const ES& es) const
{
  gcls.adjoin(Ariadne::fuzzy_outer_approximation(es,gcls.grid()));
}

template<class ES>
void
FastApproximator<ES>::adjoin_outer_approximation(GridMaskSet<R>& gms, const ES& es) const
{
  gms.adjoin(Ariadne::fuzzy_outer_approximation(es,gms.grid()));
}

template<class ES>
std::ostream&
FastApproximator<ES>::write(std::ostream& os) const
{
  return os<<"FastApproximator()\n";
}

} // namespace Ariadne
