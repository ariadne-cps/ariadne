/***************************************************************************
 *            taylor_flow.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file taylor_flow.h
 *  \brief Interface for flows in Euclidean space. 
 */

#ifndef ARIADNE_TAYLOR_FLOW_H
#define ARIADNE_TAYLOR_FLOW_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "function/affine_variable.h"
#include "function/taylor_series.h"
#include "system/flow_interface.h"

namespace Ariadne {
  
   
    /*! \brief Flows defined by Taylor series 
     *  \ingroup Integrate
     */
    template<class R>
    class TaylorFlow
      : public FlowInterface<R>
    {
      typedef Interval<R> I;
     public:
      TaylorFlow(const TaylorSeries< AffineDerivative<I> >&);
      virtual Point<I> evaluate(const I& t, const Point<I>& x);
      virtual Vector<I> tangent(const I& t, const Point<I>& pt);
      virtual Matrix<I> jacobian(const I& t, const Point<I>& pt);
     private:
      TaylorSeries< AffineVariable<I> > _data;
    };
    
  }
}

#endif /* ARIADNE_TAYLOR_FLOW_H */
