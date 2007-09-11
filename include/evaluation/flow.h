/***************************************************************************
 *            vector_field_evolver.h
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
 
/*! \file flow.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_FLOW_H
#define ARIADNE_FLOW_H

#include "../base/types.h"
#include "../base/declarations.h"
#include "../numeric/declarations.h"
#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"

namespace Ariadne {

  namespace Evaluation {

   
    /*! \brief %Base class for integration schemes. 
     *  \ingroup Integrate
     */
    template<class R>
    class Flow {
      typedef Numeric::Interval<R> I;
     public:
      Flow(const VectorFieldEvolver<R>& i, const System::VectorFieldInterface<R>& vf);
      Geometry::Point<I> evaluate(const Geometry::Point<I>& pt);
      Geometry::Point<I> jacobian(const Geometry::Point<I>& pt);
      
      
    };
    
  }
}

#endif /* ARIADNE_FLOW_H */
