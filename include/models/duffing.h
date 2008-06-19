/***************************************************************************
 *            duffing.h
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

/*! \file duffing.h
 *  \brief The Duffing equation \f$(\ddot{x}+\delta\dot{x}+(\beta x^3\pm\omega_0^2)=\gamma\cos(\omega t+\phi)\f$.
 */

#ifndef ARIADNE_DUFFING_EQUATION_H
#define ARIADNE_DUFFING_EQUATION_H

#include "system/build_vector_field.h"

namespace Ariadne {
  namespace Models {  


    template<class R, class A, class P>
    void
    duffing_function(R& r, const A& x, const P& p)
    {
      r[0]=x[1];
      r[1]=-p[0]*x[1]-x[0]*(p[2]+p[1]*x[0]*x[0])+p[3]*cos(p[4]*x[2]+p[5]);
      r[2]=1.0;
    }

    /*! \class DuffingEquation
     *  \brief The Duffing equation. 
     *  Variables: x, v, t
     *  Parameters: delta, beta, alpha, gamma, omega, phi
     *  System: dotx=v; 
     *          dotv=-delta*v-x*(alpha+beta*x*x)+gamma*cos(omega*t+phi);
     */

    ARIADNE_BUILD_VECTOR_FIELD(DuffingEquation,duffing_function,3,6,255);

    template<class R>
    std::ostream& operator<<(std::ostream& os, const DuffingEquation<R>& de) {
      os << "DuffingEquation( delta=" << de.parameter(0) << ", beta=" << de.parameter(1) << ", alpha=" << de.parameter(2) 
         << ", gamma=" << de.parameter(3) << ", omega=" << de.parameter(4) << ", phi=" << de.parameter(5) << " )";
      return os;
    }
    

     
    
  } // namespace Models
} // namespace Ariadne


#endif /* ARIADNE_DUFFING_EQUATION_H */
