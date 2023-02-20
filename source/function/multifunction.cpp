/***************************************************************************
 *            multifunction.cpp
 *
 *  Copyright 2008-21  Pieter Collins
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

#include "multifunction.hpp"

#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "function/taylor_function.hpp"
#include "function/function.hpp"
#include "function/function_model.hpp"

namespace Ariadne {

#warning
CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> product(Vector<ValidatedNumber> x1, BoxDomainType bx2) {
    return {x1,bx2};
}

using ValidatedImageSet = ValidatedConstrainedImageSet;

auto
image(CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> const& s, ValidatedVectorMultivariateFunction const& f) -> ValidatedImageSet
{
    using SIG=RealVector(RealVector);

    auto x=std::get<0>(s);
    auto pd=std::get<1>(s);

    auto fx=ValidatedFunction<SIG>::constant(pd.dimension(),x);
    auto fid=ValidatedFunction<SIG>::identity(pd.dimension());
    auto fim=compose(f,join(fx,fid));
    return ValidatedImageSet(pd,fim);
}

template<class P, class PR> ParametrisedMultifunctionModel<P,RealVector(RealVector),PR>::
    ParametrisedMultifunctionModel(BoxDomainType dom, VectorMultivariateFunction<P> f, BoxDomainType params, Sweeper<FLT> swp)
        : _f(VectorMultivariateTaylorFunctionModel<P,FLT>(product(dom,params),f,swp)), _as(dom.dimension()) { }

template<class P, class PR> ParametrisedMultifunctionModel<P,RealVector(RealVector),PR>::
    ParametrisedMultifunctionModel(SizeType as, FunctionModel<P,SIG,PR> const& f)
        : _f(f), _as(as)
{
    ARIADNE_ASSERT(as<=f.argument_size());
}

template class Multifunction<ValidatedTag,RealVector(RealVector)>;
template class ParametrisedMultifunction<ValidatedVectorMultivariateFunction,BoxDomainType>;

template class ParametrisedMultifunctionModel<ValidatedTag,RealVector(RealVector),DoublePrecision>;
template class ParametrisedMultifunctionModel<ValidatedTag,RealVector(RealVector),MultiplePrecision>;


#warning
template<class T> LowerKleenean inside(EffectiveCompactSet<T> const&, EffectiveOpenSet<T> const&) { ARIADNE_NOT_IMPLEMENTED; }

//! An open set of the form \f$\{ x : \mathbb{X} \mid F(x) \subset V \}\f$ where \f$f:X\to\mathcal{K}(\mathbb{Y})\f$ and \f$V:\mathcal{O}(\mathbb{Y})\f$,
//! for \f$\mathbb{X}\f$ being \a ARG and \f$\mathbb{Y}\f$ being \a RES.
template<class P, class ARG, class RES> class StrongPreimage;

template<class ARG, class RES> class StrongPreimage<EffectiveTag,ARG,RES>
    : public virtual OpenSetInterface<EffectiveTag,ARG>
{
    using P=EffectiveTag;
    OpenSet<P,RES> _s; Multifunction<P,RES(ARG),CompactSet> _f;
  public:
    typedef typename SetTraits<ARG>::BasicSetType BasicSetType;
    StrongPreimage(OpenSet<P,RES> s, Multifunction<P,RES(ARG),CompactSet> f) : _s(s), _f(f) { }
    DimensionType dimension() const override { return _f.argument_size(); }
    LowerKleenean contains(ARG pt) const{ return inside(_f(pt),_s); }
    LowerKleenean covers(BasicSetType const& bx) const override { ARIADNE_NOT_IMPLEMENTED; }
    ValidatedLowerKleenean covers(BasicSetType const& pt, Effort eff) const override { ARIADNE_NOT_IMPLEMENTED; }
    friend OutputStream& operator<<(OutputStream& os, StrongPreimage<P,ARG,RES> const& self) {
        return os << "strong_preimage(" << self._s << "," << self._f << ")"; }
  private:
    OutputStream& _write(OutputStream& os) const override { return os << *this; }
};


static_assert(Same<typename SetTraits<Real>::DimensionType,typename Multifunction<EffectiveTag,Real(RealVector),CompactSet>::ResultSizeType>);

template<class ARG, class RES> class MultiImageSet<EffectiveTag,RES,ARG,CompactSet>
    : public virtual CompactSetInterface<EffectiveTag,RES>
{
    using P=EffectiveTag; using SIG=RES(ARG);
    CompactSet<P,ARG> _s; Multifunction<P,SIG,CompactSet> _f;
  public:
    typedef typename SetTraits<RES>::DimensionType DimensionType;
    typedef typename SetTraits<RES>::BoundingSetType BoundingSetType;
    typedef typename SetTraits<RES>::BasicSetType BasicSetType;
    MultiImageSet(CompactSet<P,ARG> s, Multifunction<P,SIG,CompactSet> f)
        : _s(s), _f(f) { }
    virtual DimensionType dimension() const final { return this->_f.result_size(); }
#warning
    virtual BoundingSetType bounding_box() const final { return image(_s,_f).bounding_box(); }
    virtual LowerKleenean separated(BasicSetType const&) const final { ARIADNE_NOT_IMPLEMENTED; }
    virtual ValidatedLowerKleenean separated(BasicSetType const&, Effort eff) const final { ARIADNE_NOT_IMPLEMENTED; }
    virtual LowerKleenean inside(BasicSetType const&) const final { ARIADNE_NOT_IMPLEMENTED; }
    virtual ValidatedLowerKleenean inside(BasicSetType const&, Effort eff) const final { ARIADNE_NOT_IMPLEMENTED; }
    LowerKleenean inside(EffectiveOpenSet<RES> const& ops) const {
        EffectiveOpenSet<ARG> preops=StrongPreimage<P,ARG,RES>(ops,_f);
        return inside(this->_s,preops);
    }
    friend OutputStream& operator<<(OutputStream& os, MultiImageSet<P,RES,ARG,CompactSet> const& self) {
        return os << "image(" << self._s <<"," << self._f << ")"; }
    operator CompactSet<P,RES> () const { return new MultiImageSet<P,RES,ARG,CompactSet>(*this); }
  private:
    MultiImageSet<P,RES,ARG,CompactSet>* clone() const { return new MultiImageSet<P,RES,ARG,CompactSet>(*this); }
    OutputStream& _write(OutputStream& os) const { return os << *this; }
};


template<class P, class RES, class ARG> CompactSet<P,RES> image(CompactSet<P,ARG> s, CompactMultifunction<P,RES(ARG)> f);

template<class ARG, class RES> class MultiImageSet<ValidatedTag,RES,ARG,CompactSet>
    : public virtual CompactSetInterface<ValidatedTag,RES>
{
    using P=ValidatedTag; using SIG=RES(ARG);
    CompactSet<P,ARG> _s; Multifunction<P,SIG,CompactSet> _f;
  public:
    typedef typename SetTraits<RES>::DimensionType DimensionType;
    typedef typename SetTraits<RES>::BoundingSetType BoundingSetType;
    typedef typename SetTraits<RES>::BasicSetType BasicSetType;
    MultiImageSet(CompactSet<P,ARG> s, Multifunction<P,SIG,CompactSet> f)
        : _s(s), _f(f) { }
    virtual DimensionType dimension() const final { return this->_f.result_size(); }
    virtual BoundingSetType bounding_box() const final { return image(_s,_f).bounding_box(); }
    virtual ValidatedLowerKleenean separated(BasicSetType const&) const final { ARIADNE_NOT_IMPLEMENTED; }
    virtual ValidatedLowerKleenean inside(BasicSetType const&) const final { ARIADNE_NOT_IMPLEMENTED; }
    ValidatedLowerKleenean inside(ValidatedOpenSet<RES> const& ops) const { ARIADNE_NOT_IMPLEMENTED; }
    friend OutputStream& operator<<(OutputStream& os, MultiImageSet<P,RES,ARG,CompactSet> const& self) {
        return os << "image(" << self._s <<"," << self._f << ")"; }
    operator CompactSet<P,RES> () const { return new MultiImageSet<P,RES,ARG,CompactSet>(*this); }
  private:
    MultiImageSet<P,RES,ARG,CompactSet>* clone() const { return new MultiImageSet<P,RES,ARG,CompactSet>(*this); }
    OutputStream& _write(OutputStream& os) const { return os << *this; }
};


template class StrongPreimage<EffectiveTag,RealVector,Real>;
template class StrongPreimage<EffectiveTag,RealVector,RealVector>;

template<class P, class RES, class ARG> CompactSet<P,RES> image(CompactSet<P,ARG> s, CompactMultifunction<P,RES(ARG)> f) {
    return CompactMultiImageSet<P,RES,ARG>(s,f);
}

template EffectiveCompactSet<Real> image(EffectiveCompactSet<RealVector>,EffectiveMultifunction<Real(RealVector),CompactSet> f);
template ValidatedCompactSet<Real> image(ValidatedCompactSet<RealVector>,ValidatedMultifunction<Real(RealVector),CompactSet> f);

ValidatedVectorMultivariateFunction embed(SizeType n, ValidatedVectorMultivariateFunction f);
ValidatedVectorMultivariateFunction embed(ValidatedVectorMultivariateFunction f, SizeType n);
ValidatedScalarMultivariateFunction embed(ValidatedScalarMultivariateFunction f, SizeType n) {
    List<ValidatedScalarMultivariateFunction> fcs=ValidatedScalarMultivariateFunction::coordinates(f.argument_size()+n);
    fcs.resize(f.argument_size());
    ValidatedVectorMultivariateFunction p(fcs);
    return compose(f,p);
}


ValidatedConstrainedImageSet product(ValidatedConstrainedImageSet s, BoxDomainType const& bx) {
    auto rdom = product(s.domain(),bx);
    auto rfn = combine(s.function(),ValidatedVectorMultivariateFunction::identity(bx.dimension()));
    List<ValidatedConstraint> rcs = {};
    for (SizeType i=0; i!=s.number_of_constraints(); ++i) {
        ValidatedConstraint c = s.constraint(i);
        rcs.append(ValidatedConstraint(embed(c.function(),bx.dimension()),c.bounds()));
    }
    return ValidatedConstrainedImageSet(rdom,rfn,rcs);
};

ValidatedConstrainedImageSet image(ValidatedConstrainedImageSet const& s, ValidatedVectorMultivariateParametrisedPatchMultifunction const& mf) {
    ARIADNE_ASSERT(s.dimension()==mf.argument_size());
    ValidatedConstrainedImageSet sp=product(s,mf.parameter_domain());
    return ValidatedConstrainedImageSet(sp.domain(),compose(mf.argument_parameter_function(),sp.function()),sp.constraints());
}


} // namespace Ariadne

