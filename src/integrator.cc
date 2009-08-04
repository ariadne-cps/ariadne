/***************************************************************************
 *            integrator.cc
 *
 *  Copyright  2006-9  Pieter Collins
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

#include <iomanip>

#include "integrator.h"

#include "logging.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "function_interface.h"
#include "taylor_function.h"

namespace Ariadne {

typedef Vector<Interval> IVector;

Pair<Float,IVector>
IntegratorBase::flow_bounds(const FunctionInterface& vf, const IVector& dp, const IVector& dx, const Float& hmax) const
{

    ARIADNE_ASSERT(vf.result_size()==dx.size());
    ARIADNE_ASSERT(vf.argument_size()==dp.size()+dx.size());
    ARIADNE_ASSERT(hmax>0);

    // Set up constants of the method.
    // TODO: Better estimates of constants
    const double INITIAL_MULTIPLIER=2;
    const double MULTIPLIER=1.125;
    const double BOX_RADIUS_MULTIPLIER=1.03125;
    const uint EXPANSION_STEPS=8;
    const uint REDUCTION_STEPS=8;
    const uint REFINEMENT_STEPS=4;

    Vector<Interval> delta=(dx-midpoint(dx))*(BOX_RADIUS_MULTIPLIER-1);

    Float h=hmax;
    Float hmin=hmax/(1<<REDUCTION_STEPS);
    bool success=false;
    Vector<Interval> bx,nbx,df;
    Interval ih(0,h);
    while(!success) {
        ARIADNE_ASSERT_MSG(h>hmin," h="<<h<<", hmin="<<hmin);
        bx=dx+INITIAL_MULTIPLIER*ih*vf.evaluate(join(dp,dx))+delta;
        for(uint i=0; i!=EXPANSION_STEPS; ++i) {
            df=vf.evaluate(join(dp,bx));
            nbx=dx+delta+ih*df;
            if(subset(nbx,bx)) {
                success=true;
                break;
            } else {
                bx=dx+delta+MULTIPLIER*ih*df;
            }
        }
        if(!success) {
            h/=2;
            ih=Interval(0,h);
        }
    }

    ARIADNE_ASSERT(subset(nbx,bx));

    Vector<Interval> vfbx;
    vfbx=vf.evaluate(bx);

    for(uint i=0; i!=REFINEMENT_STEPS; ++i) {
        bx=nbx;
        vfbx=vf.evaluate(join(dp,bx));
        nbx=dx+delta+ih*vfbx;
        ARIADNE_ASSERT_MSG(subset(nbx,bx),std::setprecision(20)<<"refinement "<<i<<": "<<nbx<<" is not a inside of "<<bx);
    }


    // Check result of operation
    // We use subset rather than inner subset here since the bound may touch
    ARIADNE_ASSERT(subset(nbx,bx));

    bx=nbx;


    ARIADNE_ASSERT(subset(dx,bx));

    ARIADNE_ASSERT_MSG(subset(dx+h*vf.evaluate(join(dp,bx)),bx),
        "d="<<dx<<"\nh="<<h<<"\nf(b)="<<vf.evaluate(join(dp,bx))<<"\nd+hf(b)="<<Vector<Interval>(dx+h*vf.evaluate(join(dp,bx)))<<"\nb="<<bx<<"\n");

    return std::make_pair(h,bx);
}



TaylorFunction
TaylorIntegrator::flow(const FunctionInterface& vf, const IVector& dp, const IVector& dx, Float hmax) const
{
    IVector bx;
    Float h;
    make_lpair(h,bx)=this->flow_bounds(vf,dp,dx,hmax);
    //return parameterised_flow(vf,dp,dx,h,this->temporal_order());
}


} // namespace Ariadne
