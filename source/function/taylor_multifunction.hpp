/***************************************************************************
 *            taylor_multifunction.hpp
 *
 *  Copyright 2012-21  Pieter Collins
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

/*! \file taylor_multifunction.hpp
 *  \brief Multivalued functions based on Taylor polynomials
 */

#ifndef ARIADNE_TAYLOR_MULTIFUNCTION_HPP
#define ARIADNE_TAYLOR_MULTIFUNCTION_HPP


#include "taylor_model.hpp"
#include "taylor_function.hpp"

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



namespace Ariadne {


inline ValidatedLowerKleenean operator>=(LowerBound<FloatMP> const& x1, FloatDP const& x2) {
    return cast_exact(x1)>=x2; }
inline ValidatedLowerKleenean operator<=(UpperBound<FloatMP> const& x1, FloatDP const& x2) {
    return cast_exact(x1)<=x2; }

inline decltype(auto) contains(IntervalDomainType const& ivl, Interval<UpperBound<FloatMP>> const& x) {
    return subset(x,ivl); }

template<class M> class ScaledFunctionPatchMixin;

template<class U> class IntervalSetBase : public Interval<U> {
  protected:
    IntervalSetBase(Interval<U>&& ivl) : Interval<U>(ivl) { }
    template<class UU> decltype(auto) _inside(Interval<UU> const& ivl) const { return inside(*this,ivl); }
    template<class UU> decltype(auto) _separated(Interval<UU> const& ivl) const { return separated(*this,ivl); }
};

template<class U> class IntervalSet : public IntervalSetBase<U> {
  public:
    IntervalSet(Interval<U>&& ivl) : IntervalSetBase<U>(std::move(ivl)) { }
    Interval<U> const& bounding_box() const { return *this; }
    template<class UU> decltype(auto) inside(Interval<UU> const& ivl) const { return this->_inside(ivl); }
    template<class UU> decltype(auto) separated(Interval<UU> const& ivl) const { return this->_separated(ivl); }
};



template<class FLT> class ScaledFunctionPatchFactory<ValidatedIntervalTaylorModel<FLT>> {
  public:
    ScaledFunctionPatchFactory(Sweeper<FLT> const&) { }
};

template<class FLT> class ScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<FLT>>
    : public MultifunctionInterface<ValidatedTag, Real(RealVector),CompactSet>
{
    using M=ValidatedIntervalTaylorModel<FLT>; using P=ValidatedTag; using RES=Real;
    ScaledFunctionPatch<M> const& _upcast() const { return static_cast<ScaledFunctionPatch<M>const&>(*this); }
  protected:
    virtual CompactSet<P,RES> _call(Argument<Number<P>> const& x) const;
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_upcast(); }
};

template<class FLT> class VectorScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<FLT>>
    : public MultifunctionInterface<ValidatedTag, RealVector(RealVector), CompactSet>
{
    using M=ValidatedIntervalTaylorModel<FLT>; using P=ValidatedTag; using RES=RealVector;
    VectorScaledFunctionPatch<M> const& _upcast() const { return static_cast<VectorScaledFunctionPatch<M>const&>(*this); }
  protected:
    virtual CompactSet<P,RES> _call(Argument<Number<P>> const& x) const;
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_upcast(); }
};


// Declare classes and type aliases
template<class P, class SIG, class FLT> class IntervalTaylorFunctionModel;
template<class P, class FLT> using ScalarMultivariateIntervalTaylorFunctionModel = IntervalTaylorFunctionModel<P,RealScalar(RealVector),FLT>;
template<class P, class FLT> using VectorMultivariateIntervalTaylorFunctionModel = IntervalTaylorFunctionModel<P,RealVector(RealVector),FLT>;
template<class SIG, class FLT> using ValidatedIntervalTaylorFunctionModel = IntervalTaylorFunctionModel<ValidatedTag,SIG,FLT>;
template<class FLT> using ValidatedScalarMultivariateIntervalTaylorFunctionModel = IntervalTaylorFunctionModel<ValidatedTag,RealScalar(RealVector),FLT>;
template<class FLT> using ValidatedVectorMultivariateIntervalTaylorFunctionModel = IntervalTaylorFunctionModel<ValidatedTag,RealVector(RealVector),FLT>;

template<class P, class FLT> class VectorTaylorFunctionModelImageSet;
template<class FLT> using ValidatedVectorTaylorFunctionModelImageSet = VectorTaylorFunctionModelImageSet<ValidatedTag,FLT>;

template<class P, class FLT> class VectorIntervalTaylorFunctionModelImageSet;
template<class FLT> using ValidatedVectorIntervalTaylorFunctionModelImageSet = VectorIntervalTaylorFunctionModelImageSet<ValidatedTag,FLT>;

template<class FN, class DOM> class ParametrisedMultifunction;
template<class P, class FLT> using VectorMultivariateTaylorParametrisedMultifunction
    = ParametrisedMultifunction<VectorMultivariateTaylorFunctionModel<P,FLT>,BoxDomainType>;
template<class FLT> using ValidatedVectorMultivariateTaylorParametrisedMultifunction
    = ParametrisedMultifunction<ValidatedVectorMultivariateTaylorFunctionModel<FLT>,BoxDomainType>;
template<class P, class FLT> using VectorMultivariateIntervalTaylorParametrisedMultifunction
    = ParametrisedMultifunction<VectorMultivariateIntervalTaylorFunctionModel<P,FLT>,BoxDomainType>;
template<class FLT> using ValidatedVectorMultivariateIntervalTaylorParametrisedMultifunction
    = ParametrisedMultifunction<ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT>,BoxDomainType>;





template<class FLT> auto
image(CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> const& s, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> const& f)
    -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FLT>;

template<class FLT> auto
image(BoxDomainType const& s, ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> const& f)
    -> ValidatedVectorIntervalTaylorFunctionModelImageSet<FLT>;




//! \ingroup Multifunction
//! \brief A multifunction formed as \f$F(x) = \bigcup_{p\in P}\widehat{F}(x,p)\f$ for \f$\widehat{F}\f$ a Taylor multifunction and \f$P\f$ a parameter domain.
template<class P, class FLT> class VectorIntervalTaylorParametrisedMultifunctionModel
{
    static_assert(SameAs<P,ValidatedTag>);

    SizeType _as;
    VectorMultivariateIntervalTaylorFunctionModel<P,FLT> _ifm;

  public:
    using BasicSetType = BoxDomainType;
    using BoundingSetType = BoxRangeType;

    explicit VectorIntervalTaylorParametrisedMultifunctionModel(SizeType as, VectorMultivariateIntervalTaylorFunctionModel<P,FLT> ifm)
        : _as(as), _ifm(ifm) { }
    explicit VectorIntervalTaylorParametrisedMultifunctionModel(VectorMultivariateIntervalTaylorFunctionModel<P,FLT> ifm, BoxDomainType pdom);

    SizeType result_size() const { return this->_ifm.result_size(); }
    SizeType argument_size() const { return this->_as; }

    VectorIntervalTaylorFunctionModelImageSet<P,FLT> operator() (Vector<ValidatedNumber> const& x) const;
#warning Only in ImageSet?
    //    friend VectorTaylorParametrisedMultifunctionModel<P,FLT> explicitly_parametrise(VectorIntervalTaylorParametrisedMultifunctionModel<P,FLT> const&, ExactDouble eps);
    friend VectorIntervalTaylorFunctionModelImageSet<P,FLT> image(VectorIntervalTaylorParametrisedMultifunctionModel<P,FLT> const& mf, BoxDomainType const& s) { return mf._image(s); }
  private:
    VectorIntervalTaylorFunctionModelImageSet<P,FLT> _image (BoxDomainType const& s) const;
};
template<class FLT> using ValidatedVectorIntervalTaylorParametrisedMultifunctionModel = VectorIntervalTaylorParametrisedMultifunctionModel<ValidatedTag,FLT>;


//! \ingroup Multifunction
//! \brief A set formed as the image of a domain under some interval Taylor multifunction
template<class P, class FLT> class VectorIntervalTaylorFunctionModelImageSet
    : public virtual CompactSetInterface<P,RealVector>
{
    static_assert(SameAs<P,ValidatedTag>);

    VectorMultivariateIntervalTaylorFunctionModel<P,FLT> _ifm;
  public:
    using BasicSetType = BoxDomainType;
    using BoundingSetType = BoxRangeType;

    explicit VectorIntervalTaylorFunctionModelImageSet(VectorMultivariateIntervalTaylorFunctionModel<P,FLT> ifm)
        : _ifm(ifm) { }
    explicit VectorIntervalTaylorFunctionModelImageSet(BoxDomainType dom, VectorMultivariateIntervalTaylorFunctionModel<P,FLT> fm)
        : _ifm(restriction(fm,dom)) { }

    virtual SizeType dimension() const override {
        return _ifm.result_size(); }
    virtual ValidatedLowerKleenean separated(BasicSetType const& bs) const override;
    virtual ValidatedLowerKleenean inside(BasicSetType const& bs) const override;
    virtual BoundingSetType bounding_box() const override;

    friend VectorTaylorFunctionModelImageSet<P,FLT> explicitly_parametrise(VectorIntervalTaylorFunctionModelImageSet<P,FLT> const& set, ExactDouble eps) {
        auto fm=explicitly_parametrise(set._ifm,eps); return VectorTaylorFunctionModelImageSet<P,FLT>(fm); }
    friend OutputStream& operator<<(OutputStream& os, VectorIntervalTaylorFunctionModelImageSet<P,FLT> const& imset) {
        return os << "ValidatedVectorIntervalTaylorFunctionModelImageSet<" << class_name<FLT>() << ">(" << imset._ifm << ")"; }
  private:
    virtual VectorIntervalTaylorFunctionModelImageSet<P,FLT>* clone() const override {
        return new VectorIntervalTaylorFunctionModelImageSet<P,FLT>(*this); }
    virtual OutputStream& _write(OutputStream& os) const override {
        return os << *this; }

};

//! \ingroup Multifunction
//! \brief A set formed as the image of a domain under some interval Taylor multifunction
template<class P, class FLT> class VectorTaylorFunctionModelImageSet
    : public virtual CompactSetInterface<P,RealVector>
{
    static_assert(SameAs<P,ValidatedTag>);

    VectorMultivariateTaylorFunctionModel<P,FLT> _fm;
  public:
    using BasicSetType = BoxDomainType;
    using BoundingSetType = BoxRangeType;

    explicit VectorTaylorFunctionModelImageSet(VectorMultivariateTaylorFunctionModel<P,FLT> ifm)
        : _fm(ifm) { }
    explicit VectorTaylorFunctionModelImageSet(BoxDomainType dom, VectorMultivariateTaylorFunctionModel<P,FLT> fm)
        : _fm(restriction(fm,dom)) { }

    virtual SizeType dimension() const override {
        return _fm.result_size(); }
    virtual ValidatedLowerKleenean separated(BasicSetType const& bs) const override;
    virtual ValidatedLowerKleenean inside(BasicSetType const& bs) const override;
    virtual BoundingSetType bounding_box() const override;

    friend OutputStream& operator<<(OutputStream& os, VectorTaylorFunctionModelImageSet<P,FLT> const& imset) {
        return os << "ValidatedVectorTaylorFunctionModelImageSet<" << class_name<FLT>() << ">(" << imset._fm << ")"; }
  private:
    virtual VectorTaylorFunctionModelImageSet<P,FLT>* clone() const override {
        return new VectorTaylorFunctionModelImageSet<P,FLT>(*this); }
    virtual OutputStream& _write(OutputStream& os) const override {
        return os << *this; }
};



template<class FLT> class IntervalTaylorFunctionModel<ValidatedTag,RealScalar(RealVector),FLT>
    : public ScaledFunctionPatch<ValidatedIntervalTaylorModel<FLT>>
{
    using P=ValidatedTag; using M=ValidatedIntervalTaylorModel<FLT>; using RES=Real;
  public:
    IntervalTaylorFunctionModel(ScaledFunctionPatch<M> const& f) : ScaledFunctionPatch<M>(f) { }
    using ScaledFunctionPatch<M>::ScaledFunctionPatch;
    Interval<UpperBound<FLT>> operator() (Vector<Bounds<FLT>> const& x) const {
        return evaluate(this->model(),Ariadne::unscale(x,this->domain())); }
    CompactSet<P,RES> operator() (Vector<Number<P>> const& x) const { return this->_call(x); }
};

template<class FLT> class IntervalTaylorFunctionModel<ValidatedTag,RealVector(RealVector),FLT>
    : public VectorScaledFunctionPatch<ValidatedIntervalTaylorModel<FLT>>
{
    using P=ValidatedTag; using M=ValidatedIntervalTaylorModel<FLT>; using RES=RealVector;
  public:
    IntervalTaylorFunctionModel(VectorScaledFunctionPatch<M> const& f) : VectorScaledFunctionPatch<M>(f) { }
    using VectorScaledFunctionPatch<M>::VectorScaledFunctionPatch;
    Box<Interval<UpperBound<FLT>>> operator() (Vector<Bounds<FLT>> const& x) const {
        return evaluate(this->models(),Ariadne::unscale(x,this->domain())); }
    CompactSet<P,RES> operator() (Vector<Number<P>> const& x) const { return this->_call(x); }
};


template<class T, class M> struct ScaledFunctionModelTypedef;
template<class M> struct ScaledFunctionModelTypedef<RealScalar,M> { typedef ScaledFunctionPatch<M> Type; };
template<class M> struct ScaledFunctionModelTypedef<RealVector,M> { typedef VectorScaledFunctionPatch<M> Type; };

//! \brief A compact set of functions over a bounded domain, represented as a polynomial with interval-valued coefficients.
template<class SIG, class FLT> class ValidatedIntervalTaylorFunctionModelSet
    : public CompactSet<ValidatedTag,SIG>::Interface
{
    using P=ValidatedTag; using M=ValidatedIntervalTaylorModel<FLT>;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    static_assert(Same<ARG,RealVector>);
    static_assert(Same<RES,RealVector> or Same<RES,RealScalar>);
    ValidatedIntervalTaylorFunctionModel<RES(ARG),FLT> _function_model;
    using ScaledFunctionModelType = typename ScaledFunctionModelTypedef<RES,M>::Type;
  public:
    ValidatedIntervalTaylorFunctionModelSet(ScaledFunctionModelType const& f) : _function_model(f) { }
    ValidatedIntervalTaylorFunctionModelSet(BoxDomainType dom, ScaledFunctionModelType const& f) : _function_model(restriction(f,dom)) { }
    ValidatedIntervalTaylorFunctionModel<SIG,FLT> multifunction() const { return this->_function_model; }
    BoxDomainType parameter_domain() const { return this->_function_model.domain(); }

    CompactSet<P,RES> operator() (Vector<Number<P>> const& x) const { return this->_function_model(x); }
    friend OutputStream& operator<<(OutputStream& os, ValidatedIntervalTaylorFunctionModelSet<SIG,FLT> const& fs) {
        return os << "FunctionSet(" << fs._function_model << ")"; }
  private:
    virtual CompactSet<P,RES> _call(Vector<Number<P>> const& x) const { return this->operator()(x); }
    virtual OutputStream& _write(OutputStream& os) const { return os << *this; }
};


template<class FLT> using ValidatedScalarMultivariateIntervalTaylorFunctionModelSet = ValidatedIntervalTaylorFunctionModelSet<RealScalar(RealVector),FLT>;
template<class FLT> using ValidatedVectorMultivariateIntervalTaylorFunctionModelSet = ValidatedIntervalTaylorFunctionModelSet<RealVector(RealVector),FLT>;
using ValidatedVectorMultivariateIntervalTaylorFunctionModelSetDP = ValidatedIntervalTaylorFunctionModelSet<RealVector(RealVector),FloatDP>;
using ValidatedVectorMultivariateIntervalTaylorFunctionModelSetMP = ValidatedIntervalTaylorFunctionModelSet<RealVector(RealVector),FloatMP>;


template<class FLT> using ValidatedVectorMultivariateTaylorFunctionModel = VectorMultivariateTaylorFunctionModel<ValidatedTag,FLT>;

template<class FLT>
ValidatedVectorMultivariateTaylorFunctionModel<FLT>
explicitly_parametrise(ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> const& mf, ExactDouble threshold_radius);



template<class FLT> class ValidatedIntervalTaylorCurriedTrajectorySet
    : public CompactSet<ValidatedTag,RealVector(Real)>::Interface
{
    using P=ValidatedTag; using M=ValidatedIntervalTaylorModel<FLT>; using RES=RealVector;
    ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> _flow_model;
    Vector<Number<P>> _initial_state;
  public:
    ValidatedIntervalTaylorCurriedTrajectorySet(ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> const& phi, Vector<Number<P>> x0)
        : _flow_model(phi), _initial_state(x0) { ARIADNE_PRECONDITION(phi.argument_size()==x0.size()+1u); }
    CompactSet<P,RES> operator() (Scalar<Number<P>> const& t) const { return this->_flow_model(join(this->_initial_state,t)); }
    friend OutputStream& operator<<(OutputStream& os, ValidatedIntervalTaylorCurriedTrajectorySet<FLT> const& trs) {
        return os << "TrajectorySet(" << trs._flow_model << ", " << trs._initial_state << ")"; }
  private:
    virtual CompactSet<P,RES> _call(Scalar<Number<P>> const& t) const override {
        return this->operator()(t); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << *this; }
};



template<class FLT> class ValidatedIntervalTaylorMultiflow
    : public ValidatedCompactMultiflow::Interface
{
    using P=ValidatedTag; using M=ValidatedIntervalTaylorModel<FLT>; using RES=Real;
    typedef ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> FlowModelType;
    typedef ValidatedIntervalTaylorCurriedTrajectorySet<FLT> TrajectorySetType;
    typedef ValidatedVectorMultivariateIntervalTaylorFunctionModel<FLT> ReachMapType;
    FlowModelType _flow_model;
  public:
    ValidatedIntervalTaylorMultiflow(VectorScaledFunctionPatch<M> const& f) : _flow_model(f) { }
    TrajectorySetType operator() (Vector<Number<P>> const& x0) const {
        return TrajectorySetType(this->_flow_model,x0); }
    ReachMapType curry() const {
        return ReachMapType(this->_flow_model); }
    friend OutputStream& operator<<(OutputStream& os, ValidatedIntervalTaylorMultiflow<FLT> const& phi) { return phi._write(os); }
  private: public:
    virtual CompactSet<P,RealVector(Real)> _call(Vector<Number<P>> const& x0) const override {
        return CompactSet<P,RealVector(Real)>(std::make_shared<TrajectorySetType>(this->_flow_model,x0)); }
    virtual Function<P,CompactSet<P,RealVector>(RealVector)> _curry() const override {
        return Function<P,CompactSet<P,RealVector>(RealVector)>(std::make_shared<ReachMapType>(this->curry())); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << *this; }
};

} // namespace Ariadne

#endif

