/***************************************************************************
 *            vanderpol.h
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

/*! \file vanderpol.h
 *  \brief The van der Pol equation \f$\ddot{y}-\mu(1-y^2)\dot{y}+y=0\f$.
 */

#ifndef ARIADNE_VANDERPOL_EQUATION_H
#define ARIADNE_VANDERPOL_EQUATION_H

#include "system/build_vector_field.h"

namespace Ariadne {
  namespace System {
  

    template<class R, class A, class P>
    void
    van_der_pol_function(R& r, const A& x, const P& p)
    {
      r[0]=x[1];
      r[1]=p[0]*(1-x[0]*x[0])*x[1]-x[0];
    }

    /*! \class VanDerPolEquation
     *  \brief The van der Pol equation \f$\ddot{x}=\mu*(1-x^2)*\dot{x}-x\f$.
     *     Variables:  x, v
     *     Parameters: mu
     *     System:     dotx=v
     *                 dotv=mu*(1-x*x)*v-x
     */

    ARIADNE_BUILD_VECTOR_FIELD(VanDerPolEquation,van_der_pol_function,3,1,255)

    template<class R>
    std::ostream& operator<<(std::ostream& os, const VanDerPolEquation<R>& vdp) {
      os << "VanDerPolEquation( mu=" << vdp.parameter(0) << " )";
      return os;
    }



    
    
  }
}


#endif /* ARIADNE_VANDERPOL_EQUATION_H */
