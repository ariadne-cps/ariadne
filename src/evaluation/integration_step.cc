/***************************************************************************
 *            integration_step.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "base/numerical_type.h"
#include "evaluation/integration_step.h"
#include "evaluation/integration_step.tpl"

#include "real_typedef.h"

namespace Ariadne {
  namespace Evaluation {

    template 
    Geometry::Rectangle<Real> 
    integration_step(const VectorField<Real>&, const Geometry::Rectangle<Real>&, Real& step_size);
  
    template 
    Geometry::Parallelotope<Real> 
    integration_step(const VectorField<Real>&, const Geometry::Parallelotope<Real>&, Real& step_size);
  
    template 
    Geometry::Parallelotope<Real> 
    integration_step(const AffineVectorField<Real>&, const Geometry::Parallelotope<Real>&, Real& step_size);
  
    
    template 
    Geometry::Rectangle<Real> 
    reach_step(const VectorField<Real>&, const Geometry::Rectangle<Real>&, Real& step_size);
  
    template 
    Geometry::Parallelotope<Real> 
    reach_step(const VectorField<Real>&, const Geometry::Parallelotope<Real>&, Real& step_size);
  
    template 
    Geometry::Zonotope<Real> 
    reach_step(const VectorField<Real>&, const Geometry::Zonotope<Real>&, Real& step_size);
  

  }
}
