/***************************************************************************
 *            models.h
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
 *
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
 
/*! \file models.h
 *  \brief Simple models for testing.
 */

#ifndef ARIADNE_MODELS_H
#define ARIADNE_MODELS_H


#include "function.h"
#include "hybrid_automaton.h"

namespace Ariadne {

static const uint SMOOTH=255;

struct Henon : FunctionData<2,2,2> {
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=p[0]-x[0]*x[0]-p[1]*x[1];
    r[1]=x[0];
  }
};
                   
struct HenonInverse : FunctionData<2,2,2> {
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=x[1]; 
    r[1]=(p[0]-x[1]*x[1]-x[0])/p[1];
  }
};
                   
  
/*! \class Duffing
 *  \brief The Duffing equation. 
 *  Variables: x, v, t
 *  Parameters: delta, beta, alpha, gamma, omega, phi
 *  System: dotx=v; 
 *          dotv=-delta*v-x*(alpha+beta*x*x)+gamma*cos(omega*t+phi);
 */
struct Duffing : FunctionData<3,3,6> {
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=x[1];
    r[1]=-p[0]*x[1]-x[0]*(p[2]+p[1]*x[0]*x[0])+p[3]*cos(p[4]*x[2]+p[5]);
    r[2]=1.0;
  }
};
                   


/*! \class VanDerPol
 *  \brief The Van der Pol equation \f$\ddot{x}=\mu*(1-x^2)*\dot{x}-x\f$.
 *     Variables:  x, v
 *     Parameters: mu
 *     System:     dotx=v
 *                 dotv=mu*(1-x*x)*v-x
 */
struct VanDerPol : FunctionData<2,2,1> {
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=x[1];
    r[1]=p[0]*(1-x[0]*x[0])*x[1]-x[0];
  }
};


/*! \class VanDerPol
 *  \brief The Van der Pol equation \f$\ddot{x}=\mu*(1-x^2)*\dot{x}-x+a\sin(\omega t)\f$.
 *     Variables:  x, v, t
 *     Parameters: mu, a, omega
 *     System:     dotx=v
 *                 dotv=mu*(1-x*x)*v-x+a*sin(omega*t)
 *                 dott=1
 */
struct ForcedVanDerPol : FunctionData<3,3,3> {
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=x[1];
    r[1]=p[0]*(1.0-x[0]*x[0])*x[1]-x[0]+p[1]*sin(p[2]*x[2]);
    r[2]=1.0;
  }
};


/*! \class LorenzSystem
 *  \brief The Lorenz system \f$(\dot{x},\dot{y},\dot{z}) = (\sigma(y-x),\rho x-y-xz,-\beta z+xy)\f$. 
 *     Variables: x, y, z
 *     Parameters: beta, rho, sigma
 *     System: dotx=sigma*(y-x)
 *             doty=rho*x-y-x*z
 *             dotz=-beta*z+x*y
 *  The standard parameters for the Lorenz attractos are \f$\sigma=10\f$, \f$\beta = 8/3\f$ and \f$\rho=28\f$.
 */
struct Lorenz : FunctionData<3,3,3> {
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=p[2]*(x[1]-x[0]);
    r[1]=p[1]*x[0]-x[1]-x[0]*x[2];
    r[2]=-p[0]*x[2]+x[0]*x[1];
  }
};
  

//! The singularly forced Van der Pol oscillator.
class SingularVanDerPol
  : public HybridAutomaton
{
  SingularVanDerPol(const Vector<Interval>& p);
};


} // namespace Ariadne


#endif // ARIADNE_MODELS_H
