/***************************************************************************
 *            taylor_multifunction.cpp
 *
 *  Copyright 2008-12  Pieter Collins
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

#include "../config.hpp"

#include "taylor_multifunction.hpp"
#include "multifunction.hpp"

#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/algebra.hpp"
#include "../algebra/multi_index.hpp"
#include "../function/polynomial.hpp"
#include "../algebra/differential.hpp"
#include "../algebra/evaluate.hpp"
#include "../geometry/interval.hpp"
#include "../geometry/box.hpp"
#include "../geometry/set_wrapper.hpp"

#include "taylor_model.tpl.hpp"
#include "scaled_function_patch.tpl.hpp"
#include "geometry/box.tpl.hpp"



namespace Ariadne {

template<class F> auto ScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<F>>::
_call(Argument<Number<P>> const& x) const -> CompactSet<P,RES>
{
    Argument<Bounds<F>> cx(x,this->_upcast().precision());
    Interval<UpperBound<F>> ivl=this->_upcast().operator()(cx);
    IntervalSet<UpperBound<F>> ivls(std::move(ivl));
    return CompactSet<P,Real>(std::make_shared<CompactSetWrapper<IntervalSet<UpperBound<F>>,ValidatedTag,Real>>(ivls));
}

template<class F> auto VectorScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<F>>::
_call(Argument<Number<P>> const& x) const -> CompactSet<P,RES>
{
    typedef Box<Interval<UpperBound<F>>> BoxType;
    Vector<Bounds<F>> cx(x,this->_upcast().precision());
    BoxType rbx=this->_upcast().operator()(cx);
    return CompactSet<P,RealVector>(std::make_shared<CompactSetWrapper<BoxType,ValidatedTag,RealVector>>(rbx));
}

template class ScaledFunctionPatch<ValidatedIntervalTaylorModelDP>;
template class VectorScaledFunctionPatch<ValidatedIntervalTaylorModelDP>;
template class ScaledFunctionPatch<ValidatedIntervalTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedIntervalTaylorModelMP>;

template auto ScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<FloatDP>>::
_call(Argument<Number<P>> const& x) const -> CompactSet<P,RES>;
template auto ScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<FloatMP>>::
_call(Argument<Number<P>> const& x) const -> CompactSet<P,RES>;

template auto VectorScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<FloatDP>>::
_call(Argument<Number<P>> const& x) const -> CompactSet<P,RES>;
template auto VectorScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<FloatMP>>::
_call(Argument<Number<P>> const& x) const -> CompactSet<P,RES>;




template<class P, class FLT> auto
VectorIntervalTaylorParametrisedMultifunctionModel<P,FLT>::operator() (Vector<ValidatedNumber> const& x) const
    -> VectorIntervalTaylorFunctionModelImageSet<P,FLT>
{
    ARIADNE_PRECONDITION(this->argument_size()==x.size());
    Array<SizeType> a(x.size(),[](SizeType i){return i;});
    return VectorIntervalTaylorFunctionModelImageSet<P,FLT>(partial_evaluate(this->_ifm,a,x));
}

template<class P, class FLT> auto
VectorIntervalTaylorParametrisedMultifunctionModel<P,FLT>::_image(BoxDomainType const& s) const
    -> VectorIntervalTaylorFunctionModelImageSet<P,FLT>
{
    ARIADNE_PRECONDITION(this->argument_size()==s.dimension());
    BoxDomainType dom=this->_ifm.domain();
    for (SizeType i=0; i!=s.dimension(); ++i) { dom[i]=s[i]; }
    ARIADNE_PRECONDITION(subset(dom,this->_ifm.domain()));
    return VectorIntervalTaylorFunctionModelImageSet<P,FLT>(dom,this->_ifm);
}



template<class P, class FLT> auto
VectorIntervalTaylorFunctionModelImageSet<P,FLT>::separated(BasicSetType const& bs) const -> ValidatedLowerKleenean {
    VectorMultivariateTaylorFunctionModel<P,FLT> fm = explicitly_parametrise(this->_ifm, 0);
    using PR=PrecisionType<FLT>;
    return ValidatedConstrainedImageSet(fm.domain(),ValidatedVectorMultivariateFunctionModel<PR>(fm)).separated(bs);
}

template<class P, class FLT> auto
VectorIntervalTaylorFunctionModelImageSet<P,FLT>::inside(BasicSetType const& bs) const -> ValidatedLowerKleenean {
    return subset(this->_ifm.range(),bs);
}

template<class P, class FLT> auto
VectorIntervalTaylorFunctionModelImageSet<P,FLT>::bounding_box() const -> BoundingSetType {
    return static_cast<BoundingSetType>(this->_ifm.range());
}


template<class P, class FLT> auto
VectorTaylorFunctionModelImageSet<P,FLT>::separated(BasicSetType const& bs) const -> ValidatedLowerKleenean {
    std::cerr<<"separated..." <<std::endl;
    using PR=PrecisionType<FLT>;
    auto cis=ValidatedConstrainedImageSet(this->_fm.domain(),ValidatedVectorMultivariateFunctionModel<PR>(this->_fm));
    std::cerr<<"cis="<<cis<<"\n";
    ValidatedLowerKleenean r=cis.separated(bs);
    std::cerr<<"r="<<r<<"\n";
    return r;
}

template<class P, class FLT> auto
VectorTaylorFunctionModelImageSet<P,FLT>::inside(BasicSetType const& bs) const -> ValidatedLowerKleenean {
    return subset(this->_fm.range(),bs);
}

template<class P, class FLT> auto
VectorTaylorFunctionModelImageSet<P,FLT>::bounding_box() const -> BoundingSetType {
    return static_cast<BoundingSetType>(this->_fm.range());
}


template<class F, class X> F partial_evaluate(F f, Array<SizeType> const& a, Vector<X> const& x) {
    for (SizeType i=0; i!=a.size(); ++i) { f=partial_evaluate(f,a[i],x[i]); } return f; }

template<class FLT> auto
image(CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> const& s, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> const& f)
    -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FLT>
{
    Vector<ValidatedNumber> const& x=std::get<0u>(s);
    Array<SizeType> a(x.size(),[](SizeType i){return i;});
    ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> rf=partial_evaluate(f,a,x);
    return  ValidatedVectorIntervalTaylorFunctionModelImageSet<FLT>(std::get<1u>(s),rf);
}

template<class FLT> auto
image(BoxDomainType const& s, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> const& f)
    -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FLT>
{
    return  ValidatedVectorIntervalTaylorFunctionModelImageSet<FLT>(s,f);
}

template auto image(CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> const&, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatDP> const&) -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FloatDP>;
template auto image(CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> const&, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatMP> const&) -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FloatMP>;

template auto image(BoxDomainType const&, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatDP> const&) -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FloatDP>;
template auto image(BoxDomainType const&, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatMP> const&) -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FloatMP>;

template<class FLT>
ValidatedVectorMultivariateTaylorFunctionModel<FLT>
explicitly_parametrise(ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> const& mf, ExactDouble threshold_radius) {

    Vector<ValidatedIntervalTaylorModel<FLT>> const& vm=mf.models();

    std::vector<IntervalDomainType> extra_parameters;
    for (SizeType i=0; i!=vm.size(); ++i) {
        for (auto term : vm[i]) {
            if (possibly(term.coefficient().radius() > threshold_radius)) {
                extra_parameters.push_back(cast_exact_interval(static_cast<IntervalRangeType>(term.coefficient())));
            }
        }
    }

    SizeType rs=mf.size();
    SizeType as=mf.argument_size();
    SizeType ras=as + extra_parameters.size();
    BoxDomainType rdom = product(mf.domain(),BoxDomainType(Array<IntervalDomainType>(extra_parameters.begin(),extra_parameters.end())));

    Vector<ValidatedTaylorModel<FLT>> rvm(rs,ValidatedTaylorModel<FLT>(ras,vm[0].sweeper()));
    SizeType extra_parameter_number = mf.argument_size();

    MultiIndex ra(ras);
    for (SizeType i=0; i!=vm.size(); ++i) {
        auto const& m=vm[i];
        auto& rm = rvm[i];
        for (auto term : m) {
            auto const& a=term.index();
            for (SizeType j=0; j!=as; ++j) { ra[j]=a[j]; }
            Interval<UpperBound<FLT>> const& c=term.coefficient();
            FLT rc=cast_exact(c.midpoint());
            PositiveUpperBound<FLT> rce=c.radius();
            rm.expansion().append(ra,rc);
            if (possibly(rce > threshold_radius)) {
                assert(extra_parameter_number<ras);
                ra[extra_parameter_number]=1;
                rm.expansion().append(ra,cast_exact(rce));
                ra[extra_parameter_number]=0;
                ++extra_parameter_number;
            } else {
                rm.error()+=rce;
            }
        }
        rm.expansion().sort();
    }
    std::cerr<<"\nivm="<<vm<<"\nrvm="<<rvm<<"\n\n";
    return ValidatedVectorMultivariateTaylorFunctionModel<FLT>(rdom,rvm);
}

#warning
template ValidatedVectorMultivariateTaylorFunctionModel<FloatDP>
explicitly_parametrise(ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatDP> const& mf, ExactDouble threshold_radius);
template ValidatedVectorMultivariateTaylorFunctionModel<FloatMP>
explicitly_parametrise(ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatMP> const& mf, ExactDouble threshold_radius);

template class VectorIntervalTaylorFunctionModelImageSet<ValidatedTag,FloatDP>;
template class VectorIntervalTaylorFunctionModelImageSet<ValidatedTag,FloatMP>;

template class VectorTaylorFunctionModelImageSet<ValidatedTag,FloatDP>;
template class VectorTaylorFunctionModelImageSet<ValidatedTag,FloatMP>;

template class VectorIntervalTaylorParametrisedMultifunctionModel<ValidatedTag,FloatDP>;
template class VectorIntervalTaylorParametrisedMultifunctionModel<ValidatedTag,FloatMP>;

template<> String class_name<IntervalTaylorFunctionModel<ValidatedTag,RealVector(RealVector),FloatDP>>() {
    return "ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatDP>"; }
template<> String class_name<IntervalTaylorFunctionModel<ValidatedTag,RealVector(RealVector),FloatMP>>() {
    return "ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatMP>"; }


} // namespace Ariadne

