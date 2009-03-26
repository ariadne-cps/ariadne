/***************************************************************************
 *            models.cc
 *
 *  Copyright 2008  Pieter Collins
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
 
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "taylor_model.h"

#include "models.h"

namespace Ariadne {

namespace Models {

//\dot{x}      & = -x+a(x^2-1) \sin(2 pi t) 
//\dot{t} & =  \omega(x^2- 1) 
//
// r(x,t) = (2,t-0.5)
//
// g(x,t) = (x<=1)

struct SingularVanDerPolDynamic 
    : public FunctionData<2,2,2>
{
    template<class R, class A, class P> 
    void compute(R& r, const A& x, const P& p) const {
        typename R::value_type z=sqr(x[0])-1.0;
        r[0] = -x[0]+p[0]*z+sin(2*pi<Float>()*x[1]);
        r[1] = p[1] * z;
    }
};

struct SingularVanDerPolReset 
    : public FunctionData<2,2,0>
{
    template<class R, class A, class P> 
    void compute(R& r, const A& x, const P& p) const {
        r[0] = 2.0;
        r[1] = x[1]-0.5;
    }
};

struct SingularVanDerPolGuard 
    : public FunctionData<1,2,0>
{
    template<class R, class A, class P> 
    void compute(R& r, const A& x, const P& p) const {
        r[0] = 1.0-x[0];
    }
};


SingularVanDerPol::
SingularVanDerPol(const Vector<Interval>& p)
    : HybridAutomaton("SingularVanDerPolOscillator")
{
    Function<SingularVanDerPolDynamic> dynamic(p);
    Function<SingularVanDerPolReset> reset;
    Function<SingularVanDerPolGuard> guard;
    this->new_mode(1,dynamic);
    this->new_forced_transition(1,1,1,reset,guard);
} 


} // namespace Models

} // namespace Ariadne
