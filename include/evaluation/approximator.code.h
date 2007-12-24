/***************************************************************************
 *            approximator.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/grid_cell.h"
#include "geometry/grid_block.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "geometry/grid_approximation.h"
#include "evaluation/approximator.h"

namespace Ariadne {



template<class BS>
Evaluation::Approximator<BS>::~Approximator()
{
}

template<class BS>
Evaluation::Approximator<BS>::Approximator()
{
}

template<class BS>
Evaluation::Approximator<BS>::Approximator(const Approximator<BS>& approx)
{
}

template<class BS>
Evaluation::Approximator<BS>* 
Evaluation::Approximator<BS>::clone() const
{
  return new Approximator<BS>(*this);
}

template<class BS>
BS
Evaluation::Approximator<BS>::over_approximation(const Geometry::Box<R>& r) const
{
  return BS(r);
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::Approximator<BS>::outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const
{
  return Geometry::outer_approximation(bs,g);
}




template<class R>
Evaluation::Approximator< Geometry::Rectangle<R> >* 
Evaluation::Approximator< Geometry::Rectangle<R> >::clone() const
{
  return new Approximator< Geometry::Rectangle<R> >(*this);
}

template<class R>
Geometry::Rectangle<R> 
Evaluation::Approximator< Geometry::Rectangle<R> >::over_approximation(const Geometry::Box<R>& r) const
{
  return r;
}

template<class R>
Geometry::GridCellListSet<R>
Evaluation::Approximator< Geometry::Rectangle<R> >::outer_approximation(const Geometry::Rectangle<R>& bs, const Geometry::Grid<R>& g) const
{
  Geometry::GridCellListSet<R> gcls(g);
  gcls.adjoin_outer_approximation(bs);
  return gcls;
}





template<class BS>
Evaluation::FastApproximator<BS>::~FastApproximator()
{
}

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
Evaluation::FastApproximator<BS>::over_approximation(const Geometry::Box<R>& r) const
{
  return BS(r);
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::FastApproximator<BS>::outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const
{
  return Geometry::fuzzy_outer_approximation(bs,g);
}


} // namespace Ariadne
