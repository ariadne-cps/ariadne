/***************************************************************************
 *            constraint.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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


#include "numeric/float.h"

#include "geometry/constraint_interface.h"

#include "geometry/constraint.h"
#include "geometry/constraint.code.h"

#include "geometry/linear_constraint.h"
#include "geometry/linear_constraint.code.h"

#include "geometry/set_constraint.h"
#include "geometry/set_constraint.code.h"

namespace Ariadne {
  

    
    
#ifdef ENABLE_FLOAT64
    template class ConstraintInterface<Float64>;
    template class Constraint<Float64>;
    template class LinearConstraint<Float64>;
    template class SetConstraint<Float64>;
#endif
  
#ifdef ENABLE_FLOATMP
    template class ConstraintInterface<FloatMP>;
    template class Constraint<FloatMP>;
    template class LinearConstraint<FloatMP>;
    template class SetConstraint<FloatMP>;
#endif

  
} // namespace Ariadne
