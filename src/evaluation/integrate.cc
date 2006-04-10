/***************************************************************************
 *            integrate.cc
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


#include "evaluation/integrate.h"
#include "evaluation/integrate.tpl"

#include "real_typedef.h"

namespace Ariadne {
  namespace Evaluation {

    template 
    Geometry::Rectangle<Real> 
    integrate(const VectorField<Real>&, 
              const Geometry::Rectangle<Real>&, 
              const Real&,
              const Real&);

    template 
    Geometry::Parallelotope<Real> 
    integrate(const VectorField<Real>&, 
              const Geometry::Parallelotope<Real>&, 
              const Real&,
              const Real&);

    
    template 
    Geometry::ListSet<Real,Geometry::Rectangle> 
    integrate(const VectorField<Real>&, 
              const Geometry::ListSet<Real,Geometry::Rectangle>&, 
              const Real&,
              const Real&);
  
    template 
    Geometry::ListSet<Real,Geometry::Parallelotope> 
    integrate(const VectorField<Real>&, 
              const Geometry::ListSet<Real,Geometry::Parallelotope>&, 
              const Real&,
              const Real&);
  
    template 
    Geometry::GridMaskSet<Real> 
    integrate(const VectorField<Real>&, 
              const Geometry::GridMaskSet<Real>&, 
              const Geometry::GridMaskSet<Real>&,
              const Real&,
              const Real&);
  
  
    template 
    Geometry::ListSet<Real,Geometry::Parallelotope> 
    reach(const VectorField<Real>&, 
          const Geometry::ListSet<Real,Geometry::Parallelotope>&, 
          const Real&,
          const Real&);
  
    template 
    Geometry::GridMaskSet<Real> 
    reach(const VectorField<Real>&, 
          const Geometry::GridMaskSet<Real>&, 
          const Geometry::GridMaskSet<Real>&, 
          const Real&,
          const Real&);
  
  }
}
