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
    PositiveFloatDPValue h=cast_exact(hsug);

    UpperBoxType B = this->_initial(f,dom,h);

    while(not refines(this->_refinement(B,f,dom,h),B)) {
        h=hlf(h);
    }

    for(Nat i=0; i<4; ++i) {
        B = this->_refinement(B,f,dom,h);
    }

    return std::make_pair(h,B);
}

UpperBoxType BounderBase::_initial(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const {
    SizeType n = f.result_size();
    SizeType p = f.argument_size();
    BoxDomainType D = project(dom,range(0,n));
    BoxDomainType V = project(dom,range(n,p));
    UpperBoxType wD = D + (D-D.midpoint());
    return wD + 2*formula(D,V,f,UpperBoxType(dom),h);
}

UpperBoxType BounderBase::_refinement(UpperBoxType B, ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const {
    SizeType n = f.result_size();
    SizeType p = f.argument_size();
    BoxDomainType D = project(dom,range(0,n));
    BoxDomainType V = project(dom,range(n,p));
    UpperBoxType BV = product(B,UpperBoxType(V));
    return D + formula(D,V,f,BV,h);
}

UpperBoxType EulerBounder::formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const {
    return IntervalDomainType(0,h)*apply(f,B);
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
