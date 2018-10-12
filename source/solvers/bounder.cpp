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
#include "../function/formula.hpp"
#include "../function/taylor_model.hpp"

namespace Ariadne {

Pair<StepSizeType,UpperBoxType> BounderBase::compute(ValidatedVectorMultivariateFunction f, BoxDomainType dom, StepSizeType hsug) const {
    const PositiveFloatDPValue INITIAL_STARTING_WIDENING=cast_positive(2.0_exact);
    const PositiveFloatDPValue INITIAL_REFINING_WIDENING=cast_positive(1.125_exact);
    const PositiveFloatDPValue LIPSCHITZ_TOLERANCE=cast_positive(0.5_exact);
    const Nat EXPANSION_STEPS=4;
    const Nat REFINEMENT_STEPS=4;

    StepSizeType h=hsug;

    FloatDPUpperBound lipschitz = norm(f.jacobian(Vector<FloatDPBounds>(cast_singleton(dom)))).upper();
    StepSizeType hlip = static_cast<StepSizeType>(cast_exact(LIPSCHITZ_TOLERANCE/lipschitz));
    h=min(hlip,h);

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

UpperBoxType BounderBase::_initial(BoxDomainType dom, ValidatedVectorMultivariateFunction f, UpperBoxType arg, StepSizeType h, PositiveFloatDPValue FORMULA_WIDENING) const {
    const PositiveFloatDPValue BOX_RADIUS_WIDENING=cast_positive(0.25_exact);
    SizeType n = f.result_size();
    SizeType p = f.argument_size();
    BoxDomainType D = project(dom,range(0,n));
    BoxDomainType V = project(dom,range(n,p));
    UpperBoxType wD = D + BOX_RADIUS_WIDENING*(D-D.midpoint());
    return wD + FORMULA_WIDENING*formula(D,V,f,arg,h);
}

UpperBoxType BounderBase::_refinement(BoxDomainType dom, ValidatedVectorMultivariateFunction f, UpperBoxType B, StepSizeType h) const {
    SizeType n = f.result_size();
    SizeType p = f.argument_size();
    BoxDomainType D = project(dom,range(0,n));
    BoxDomainType V = project(dom,range(n,p));
    UpperBoxType BV = product(B,UpperBoxType(V));
    return D + formula(D,V,f,BV,h);
}

UpperBoxType EulerBounder::formula(BoxDomainType D, BoxDomainType V, ValidatedVectorMultivariateFunction f, UpperBoxType arg, StepSizeType h) const {
    return IntervalDomainType(0,h)*apply(f,arg);
}

} // namespace Ariadne;
