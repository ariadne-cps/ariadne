/***************************************************************************
 *            affine_applicator.code.h
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
 
#include "affine_applicator.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/affine_function.h"

#include "geometry/zonotope.h"

#include "system/affine_map.h"
#include "evaluation/standard_applicator.h"

#include "output/logging.h"

namespace Ariadne {



template<class R>
const AffineMap<R>*
AffineApplicator< Zonotope<R> >::
cast(const Map<R>* f) const
{
  if(dynamic_cast<const AffineFunction<R>*>(&f->function())) {
    return static_cast<const AffineMap<R>*>(f);
  } else {
    return 0;
  }
}

template<class R>
Zonotope<R> 
AffineApplicator< Zonotope<R> >::
apply(const Map<R>& map, 
      const Zonotope<R>& initial_set) const
{
  const AffineMap<R>* affine_map=this->cast(&map);
  if(!affine_map) {
    ARIADNE_THROW(std::runtime_error,"AffineApplicator::apply(...)","map is not affine.");
  }
  return this->apply(*affine_map,initial_set);
}




template<class R>
Zonotope<R>
AffineApplicator< Zonotope<R> >::
apply(const AffineMap<R>& f, const Zonotope<R>& z) const
{
  return StandardApplicator< Zonotope<R> >().apply(f,z);
}


} // namespace Ariadne
