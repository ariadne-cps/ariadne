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

namespace Ariadne {
  namespace Evaluation {

    template 
    Geometry::Rectangle<Dyadic> 
    integrate(const VectorField<Dyadic>&, const Geometry::Rectangle<Dyadic>&, const Dyadic& t);
  
    template 
    Geometry::Parallelotope<Dyadic> 
    integrate(const VectorField<Dyadic>&, const Geometry::Parallelotope<Dyadic>&, const Dyadic& t);
  
    template 
    Geometry::Parallelotope<Dyadic> 
    integrate(const VectorField<Dyadic>&, const Geometry::Parallelotope<Dyadic>&, const Interval<Dyadic>& t);
  
    template 
    Geometry::ListSet<Dyadic,Geometry::Parallelotope> 
    integrate(const VectorField<Dyadic>&, 
              const Geometry::ListSet<Dyadic,Geometry::Parallelotope>&, 
              const Dyadic& t);
  
  }
}
