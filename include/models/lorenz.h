/***************************************************************************
 *            lorenz_system.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.itm, Pieter.Collins@cwi.nl
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
 
/*! \file lorenz_system.h
 *  \brief The Lorenz system \f$(\dot{x},\dot{y},\dot{z}) = (\sigma(y-x),\rho x-y-xz,-\beta z+xy)\f$.
 */

#ifndef ARIADNE_LORENZ_SYSTEM_H
#define ARIADNE_LORENZ_SYSTEM_H

#include "system/build_vector_field.h"

namespace Ariadne {
  namespace System {

    template<class R, class A, class P>
    void
    lorenz_function(R& r, const A& x, const P& p) 
    {
      r[0]=p[2]*(x[1]-x[0]);
      r[1]=p[1]*x[0]-x[1]-x[0]*x[2];
      r[2]=-p[0]*x[2]+x[0]*x[1];
    }
      
    /*! \class LorenzSystem
     *  \brief The Lorenz system \f$(\dot{x},\dot{y},\dot{z}) = (\sigma(y-x),\rho x-y-xz,-\beta z+xy)\f$. 
     *     Variables: x, y, z
     *     Parameters: beta, rho, sigma
     *     System: dotx=sigma*(y-x)
     *             doty=rho*x-y-x*z
     *             dotz=-beta*z+x*y
     *  The standard parameters for the Lorenz attractos are \f$\sigma=10\f$, \f$\beta = 8/3\f$ and \f$\rho=28\f$.
     */
    ARIADNE_BUILD_VECTOR_FIELD(LorenzSystem,lorenz_function,3,3,255);

    template<class R>
    std::ostream& operator<<(std::ostream& os, const LorenzSystem<R>& ls) {
      os << "LorenzSystem( beta=" << ls.parameter(0) << ", rho=" << ls.parameter(1) << ", sigma=" << ls.parameter(2) << " )";
      return os;
    }

    
    
  }
}


#endif /* ARIADNE_LORENZ_SYSTEM_H */
