/***************************************************************************
 *            kuhn_applicator.code.h
 *
 *  Copyright  2006-7  Pieter Collins
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
 
#include "standard_applicator.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/affine_model.h"
#include "geometry/zonotope.h"

#include "system/map.h"

#include "base/stlio.h"
#include "output/logging.h"

namespace Ariadne {




template<class R>
Geometry::Zonotope<R> 
Evaluation::KuhnApplicator<R>::apply(const System::Map<R>& f, const Geometry::Zonotope<R>& z) const
{
  Function::AffineModel<R> model(z.bounding_box(),z.centre(),f.function());
  Geometry::Zonotope<R> fz=Geometry::apply(model,z);
  return Geometry::cascade_over_approximation(fz,this->_cascade_size);
}


}
