/***************************************************************************
 *            flow_interface.h
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
 
/*! \file flow_interface.h
 *  \brief Interface for flows in Euclidean space. 
 */

#ifndef ARIADNE_FLOW_INTERFACE_H
#define ARIADNE_FLOW_INTERFACE_H

#include "../base/types.h"
#include "../base/declarations.h"
#include "../numeric/declarations.h"
#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"

namespace Ariadne {

  namespace System {

   
    /*! \brief %Base class for flows. 
     *  \ingroup Integrate
     */
    template<class R>
    class FlowInterface {
      typedef Numeric::Interval<R> I;
     public:
      virtual ~FlowInterface() { }
      virtual Geometry::Point<I> evaluate(const I& t, const Geometry::Point<I>& x) = 0;
      virtual LinearAlgebra::Vector<I> tangent(const I& t, const Geometry::Point<I>& pt) = 0;
      virtual LinearAlgebra::Matrix<I> jacobian(const I& t, const Geometry::Point<I>& pt) = 0;
    };
    
  }
}

#endif /* ARIADNE_FLOW_INTERFACE_H */
