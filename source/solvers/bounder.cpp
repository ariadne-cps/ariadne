/***************************************************************************
 *            bounder.cpp
 *
 *  Copyright  2018  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "bounder.hpp"

namespace Ariadne {

Pair<PositiveFloatDPValue,UpperBoxType> BounderBase::flow_bounds(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPApproximation hsug) const {
    const PositiveFloatDPValue INITIAL_STARTING_WIDENING=cast_positive(2.0_exact);
    const PositiveFloatDPValue INITIAL_REFINING_WIDENING=cast_positive(1.125_exact);
    const PositiveFloatDPValue LIPSCHITZ_TOLERANCE=cast_positive(0.5_exact);
    const Nat EXPANSION_STEPS=4;
    const Nat REFINEMENT_STEPS=4;

    PositiveFloatDPValue h=cast_exact(hsug);

    FloatDPUpperBound lipschitz = norm(f.jacobian(Vector<FloatDPBounds>(cast_singleton(dom)))).upper();
    PositiveFloatDPValue hlip = cast_positive(cast_exact(LIPSCHITZ_TOLERANCE/lipschitz));
    h=cast_positive(min(hlip,h));

    UpperBoxType B;
    UpperBoxType V(project(dom,range(f.result_size(),f.argument_size())));
    Bool success=false;
    while(!success) {
        B=this->_initial(dom,f,UpperBoxType(dom),h,INITIAL_STARTING_WIDENING);
        for(Nat i=0; i<EXPANSION_STEPS; ++i) {
            UpperBoxType Br=this->_refinement(dom,f,B,h);
            if(not definitely(is_bounded(Br))) {
                success=false;
                break;
            } else if(refines(Br,B)) {
                B=Br;
                success=true;
                break;
            } else {
                UpperBoxType BV=product(B,V);
                B=this->_initial(dom,f,BV,h,INITIAL_REFINING_WIDENING);
            }
        }
        if(!success) {
            h=hlf(h);
        }
    }

    for(Nat i=0; i<REFINEMENT_STEPS; ++i) {
        B = this->_refinement(dom,f,B,h);
    }

    return std::make_pair(h,B);
}

UpperBoxType BounderBase::_initial(BoxDomainType dom, ValidatedVectorFunction f, UpperBoxType arg, PositiveFloatDPValue h, PositiveFloatDPValue FORMULA_WIDENING) const {
    const PositiveFloatDPValue BOX_RADIUS_WIDENING=cast_positive(0.25_exact);
    SizeType n = f.result_size();
    SizeType p = f.argument_size();
    BoxDomainType D = project(dom,range(0,n));
    BoxDomainType V = project(dom,range(n,p));
    UpperBoxType wD = D + BOX_RADIUS_WIDENING*(D-D.midpoint());
    return wD + FORMULA_WIDENING*formula(D,V,f,arg,h);
}

UpperBoxType BounderBase::_refinement(BoxDomainType dom, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const {
    SizeType n = f.result_size();
    SizeType p = f.argument_size();
    BoxDomainType D = project(dom,range(0,n));
    BoxDomainType V = project(dom,range(n,p));
    UpperBoxType BV = product(B,UpperBoxType(V));
    return D + formula(D,V,f,BV,h);
}

UpperBoxType EulerBounder::formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType arg, PositiveFloatDPValue h) const {
    return IntervalDomainType(0,h)*apply(f,arg);
}

UpperBoxType HeunBounder::formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const {
    UpperBoxType k1 = IntervalDomainType(0,h)*apply(f,B);
    UpperBoxType B2 = D + k1;
    UpperBoxType BV2 = product(B2,UpperBoxType(V));
    UpperBoxType k2 = IntervalDomainType(0,h)*apply(f,BV2);
    return (k1+k2)/2;
}

UpperBoxType RalstonBounder::formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const {
    UpperBoxType k1 = IntervalDomainType(0,h)*apply(f,B);
    UpperBoxType B2 = D+2*k1/3;
    UpperBoxType BV2 = product(B2,UpperBoxType(V));
    UpperBoxType k2 = IntervalDomainType(0,h)*apply(f,BV2);
    return (k1+3*k2)/4;
}

UpperBoxType RungeKutta4Bounder::formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const {
    UpperBoxType k1 = IntervalDomainType(0,h)*apply(f,B);
    UpperBoxType B2 = D+k1/2;
    UpperBoxType BV2 = product(B2,UpperBoxType(V));
    UpperBoxType k2 = IntervalDomainType(0,h)*apply(f,BV2);
    UpperBoxType B3 = D+k2/2;
    UpperBoxType BV3 = product(B3,UpperBoxType(V));
    UpperBoxType k3 = IntervalDomainType(0,h)*apply(f,BV3);
    UpperBoxType B4 = D+k3;
    UpperBoxType BV4 = product(B4,UpperBoxType(V));
    UpperBoxType k4 = IntervalDomainType(0,h)*apply(f,BV4);
    return (k1+2*k2+2*k3+k4)/6;
}

BounderHandler BounderFactory::create(Bounder method) {
    switch(method) {
    case Bounder::EULER : return BounderHandler(SharedPointer<BounderInterface>(new EulerBounder()));
    case Bounder::HEUN : return BounderHandler(SharedPointer<BounderInterface>(new HeunBounder()));
    case Bounder::RALSTON : return BounderHandler(SharedPointer<BounderInterface>(new RalstonBounder()));
    case Bounder::RUNGEKUTTA4 : return BounderHandler(SharedPointer<BounderInterface>(new RungeKutta4Bounder()));
    default:
        ARIADNE_FAIL_MSG("Unexpected flow bounds method "<<method<<"\n");
    }
}

} // namespace Ariadne;
