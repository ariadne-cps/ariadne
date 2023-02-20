/***************************************************************************
 *            multifunction.hpp
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

/*! \file multifunction.hpp
 *  \brief Multivalued functions
 */

#ifndef ARIADNE_MULTIFUNCTION_HPP
#define ARIADNE_MULTIFUNCTION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/declarations.hpp"
#include "utility/macros.hpp"
#include "utility/pointer.hpp"
#include "utility/container.hpp"
#include "utility/metaprogramming.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/differential.hpp"

#include "function/function_interface.hpp"
#include "function/function.hpp"
//#include "function/function_model.hpp"

#include "geometry/set.hpp"
#include "geometry/interval.hpp"
#include "geometry/box.hpp"
#include "geometry/function_set.hpp"



namespace Ariadne {

template<class P, class SIG, template<class,class>class SET> class Multifunction;


#warning
inline ValidatedLowerKleenean operator<(UpperBound<FloatMP> const& x1, LowerBound<FloatDP> const& x2) {
    if (x1.raw()<x2.raw()) { return true; } else { return ValidatedLowerKleenean(ValidatedKleenean(indeterminate)); } }
inline ValidatedLowerKleenean operator>(LowerBound<FloatMP> const& x1, UpperBound<FloatDP> const& x2) {
    if (x1.raw()>x2.raw()) { return true; } else { return ValidatedLowerKleenean(ValidatedKleenean(indeterminate)); } }
template<class F1, class F2> inline ValidatedKleenean operator<=(Bounds<F1> const& x1, Bounds<F2> const& x2) {
    if (x1.upper_raw()<=x2.lower_raw()) {  return true; } else if (x1.lower_raw()> x2.upper_raw()) { return false; } else { return indeterminate; } }
template<class F1, class F2> inline ValidatedKleenean operator>=(Bounds<F1> const& x1, Bounds<F2> const& x2) {
    if (x1.upper_raw()>=x2.lower_raw()) {  return true; } else if (x1.lower_raw()< x2.upper_raw()) { return false; } else { return indeterminate; } }
template<class F1, class F2> inline ValidatedKleenean operator<=(UpperBound<F1> const& x1, LowerBound<F2> const& x2) {
    if (x1.raw()<=x2.raw()) {  return true; } else { return indeterminate; } }
inline decltype(auto) operator<=(FloatDP const& x1, Bounds<FloatMP> const& x2) {
    return Bounds<FloatDP>(x1) <= x2;
}
inline decltype(auto) operator>=(FloatDP const& x1, Bounds<FloatMP> const& x2) {
    return Bounds<FloatDP>(x1) >= x2;
}



template<> struct ElementTraits<RealScalar> {
    typedef SizeOne SizeType;
    typedef IndexZero IndexType;
};
template<> struct ElementTraits<RealVector> {
    typedef Ariadne::SizeType SizeType;
    typedef Ariadne::SizeType IndexType;
};


template<class P> class GenericPoint {
    class Interface {
        virtual ~Interface() = 0;
        virtual Vector<Interval<Number<P>>> _bounds() const;
        GenericPoint<P> _apply(VectorFunction<P> const& f) const;
    };
    SharedPointer<const Interface> _ptr;
  public:
    GenericPoint(SharedPointer<const Interface> ptr) : _ptr(ptr) { }
    Vector<Interval<Number<P>>> bounds() { return this->_ptr->_bounds(); }
    Interval<Number<P>> operator[](SizeType i) { return this->bounds()[i]; }
    friend GenericPoint<P> image(GenericPoint<P> const& pt, VectorFunction<P> const& f) { return pt._ptr->_apply(f); }
};

template<class P> class ImageGenericPoint : public GenericPoint<P>::Interface {
    VectorFunction<P> _h;
    ImageGenericPoint(VectorFunction<P> const& h) : _h(h) { }
    friend ImageGenericPoint<P> image(GenericPoint<P> const& pt, VectorFunction<P> const& f) {
        return ImageGenericPoint<P>(f(pt._h)); }
};


template<class P, class SIG, template<class,class>class SET=LocatedSet> class Multifunction;
template<class P, class SIG> using OvertMultifunction = Multifunction<P,SIG,OvertSet>;
template<class P, class SIG> using CompactMultifunction = Multifunction<P,SIG,CompactSet>;
template<class P, class SIG> using LocatedMultifunction = Multifunction<P,SIG,LocatedSet>;
template<class SIG, template<class,class>class SET=LocatedSet> using EffectiveMultifunction = Multifunction<EffectiveTag,SIG,SET>;
template<class SIG, template<class,class>class SET=LocatedSet> using ValidatedMultifunction = Multifunction<ValidatedTag,SIG,SET>;
template<class SIG, template<class,class>class SET=LocatedSet> using ApproximateMultifunction = Multifunction<ApproximateTag,SIG,SET>;
template<class SIG> using EffectiveCompactMultifunction = Multifunction<EffectiveTag,SIG,CompactSet>;
template<class SIG> using ValidatedCompactMultifunction = Multifunction<ValidatedTag,SIG,CompactSet>;


template<class P, template<class,class>class SET=LocatedSet> using VectorMultivariateMultifunction = Multifunction<P,RealVector(RealVector),SET>;
using ValidatedVectorMultivariateMultifunction = VectorMultivariateMultifunction<ValidatedTag>;


template<class P, class SIG, template<class,class>class SET=LocatedSet> class MultifunctionInterface;

template<class P, class SIG, template<class,class>class SET> class MultifunctionInterface {
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    typedef typename ElementTraits<ARG>::SizeType ArgumentSizeType;
    typedef typename ElementTraits<RES>::SizeType ResultSizeType;
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    virtual ~MultifunctionInterface() = default;
    virtual ArgumentSizeType argument_size() const = 0;
    virtual ResultSizeType result_size() const = 0;
    virtual SET<P,RES> _call(Argument<Number<P>> const& x) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class F, class P, class SIG, template<class,class>class SET> struct IsMultifunction {
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;

    template<class FF, class=decltype(declval<SET<P,RES>>()=declval<FF>()(declval<Argument<Number<P>>>()))>
        static std::true_type test(int);
    template<class FF>
        static std::false_type test(...);
    static const bool value = decltype(test<F>(1))::value;
};

template<class F, class P, class SIG, template<class,class>class SET> concept AMultifunction = IsMultifunction<F,P,SIG,SET>::value;


//! \ingroup Multifunction
//! \brief Functions \f$\X\mvto\Y\f$ returning a set of type \ref SET at each point.
template<class P, class SIG, template<class,class>class SET> class Multifunction
    : public Handle<MultifunctionInterface<P,SIG,SET>>
{
  public:
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
    using ResultSetType = SET<P,RES>;
  public:
    typedef MultifunctionInterface<P,SIG,SET> Interface;
    typedef D DomainType;
    typedef typename ElementTraits<ARG>::SizeType ArgumentSizeType;
    typedef typename ElementTraits<RES>::SizeType ResultSizeType;
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;

    using Handle<Interface>::Handle;
    //! \brief <p/>
    template<AMultifunction<P,SIG,SET> MF> explicit Multifunction(MF const& mf);
    //! \brief <p/>
    ArgumentSizeType argument_size() const { return this->reference().argument_size(); }
    //! \brief <p/>
    ResultSizeType result_size() const { return this->reference().result_size(); }
    //! \brief <p/>
    SET<P,RES> operator() (Argument<Number<P>> const& x) const { return this->reference()._call(x); }
    //! \brief <p/>
    friend SET<P,RES> image(SET<P,ARG> const& x, Multifunction<P,SIG,SET> const& f) { ARIADNE_NOT_IMPLEMENTED; }
    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os,Multifunction<P,SIG,SET> const& mf) { return mf.reference()._write(os); }
};

static_assert(Same<typename Multifunction<EffectiveTag,Real(RealVector),CompactSet>::RES,RealScalar>);
static_assert(Same<typename Multifunction<EffectiveTag,Real(RealVector),CompactSet>::ARG,RealVector>);
static_assert(Same<typename Multifunction<EffectiveTag,Real(RealVector),CompactSet>::ResultSizeType,SizeOne>);
static_assert(Same<typename Multifunction<EffectiveTag,Real(RealVector),CompactSet>::ArgumentSizeType,SizeType>);

template<class MF, class P, class SIG, template<class,class>class SET> class MultifunctionWrapper
    : public virtual MultifunctionInterface<P,SIG,SET>
{
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
    using CD=typename SignatureTraits<SIG>::CodomainType;
  public:
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    using ArgumentSizeType = typename ElementTraits<D>::SizeType;
    using ResultSizeType = typename ElementTraits<CD>::SizeType;

    MultifunctionWrapper(MF mf) : _mf(mf) { }
    MF const& base() const { return this->_mf; }
    virtual ArgumentSizeType argument_size() const final override { return this->base().argument_size(); }
    virtual ResultSizeType result_size() const final override { return this->base().result_size(); }
    virtual SET<P,RES> _call(Argument<Number<P>> const& x) const final override {
        return static_cast<SET<P,RES>>(this->_mf(x)); }
    virtual OutputStream& _write(OutputStream& os) const  final override {
        return os << this->base(); }
  private:
    MF _mf;
};

template<class P, class SIG,template<class,class>class SET> template<AMultifunction<P,SIG,SET> MF>
Multifunction<P,SIG,SET>::Multifunction(MF const& mf)
    : Multifunction<P,SIG,SET>(std::make_shared<MultifunctionWrapper<MF,P,SIG,SET>>(mf)) {
}


#warning
template<class... TS> class CartesianProduct : public Tuple<TS...> { using Tuple<TS...>::tuple; };

CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> product(Vector<ValidatedNumber>, BoxDomainType);
ValidatedConstrainedImageSet image(CartesianProduct<Vector<ValidatedNumber>,BoxDomainType> const& set, ValidatedVectorMultivariateFunction const& f);

//! \ingroup Multifunction
//! A multifunction \f$F:\R^m \mvto \R^n\f$ defined as \f$F(x) = \{ f(x,p) \mid p \in P \}\f$ for some set \f$P\f$.
template<class F, class PDOM> class ParametrisedMultifunction {
    F _f; PDOM _pdom;
  public:
    typedef PDOM ParameterDomainType;

    ParametrisedMultifunction(F f, PDOM pdom) : _f(f), _pdom(pdom) { ARIADNE_PRECONDITION(f.argument_size()>=pdom.dimension());   }

    SizeType argument_size() const { return this->_f.argument_size()-this->_pdom.dimension(); }
    SizeType result_size() const { return this->_f.result_size(); }
    F argument_parameter_function() const { return this->_f; }
    ParameterDomainType parameter_domain() const { return this->_pdom; }

    template<class X> decltype(auto) operator() (X const& x) const {
        return image(product(x,_pdom),_f); }
    template<class S> decltype(auto) friend image(ParametrisedMultifunction<F,PDOM> const& mf, S const& s) {
        return mf._image(s); }

    friend OutputStream& operator<<(OutputStream& os, ParametrisedMultifunction<F,PDOM> const& mf) {
        return os << "ParametrisedMultifunction<" << class_name<F>() << "," << class_name<PDOM>() << ">"
                  << "(function="<<mf._f<<", parameter_domain="<<mf._pdom<<")"; }
  private:
    template<class S> decltype(auto) _image(S const& s) const {
        return image(product(s,this->_pdom),this->_f); }
};

static_assert(AMultifunction<ParametrisedMultifunction<ValidatedVectorMultivariateFunction,BoxDomainType>,ValidatedTag,RealVector(RealVector),LocatedSet>);
//! \ingroup Multifunction
//! A multifunction \f$F:\R^m \mvto \R^n\f$ defined as \f$F(x) = \{ f(x,p) \mid p \in P \}\f$ for some set \f$P\f$.
template<class DOM, class F> class FunctionImageSet {
    DOM _dom; F _f;
    typedef BoxDomainType BasicSetType ;
    typedef BoxRangeType BoundingSetType;
  public:
    FunctionImageSet(DOM dom, F f) : _dom(dom), _f(f) { }
    SizeType dimension() const { return this->_f.argument_size(); }
    ValidatedLowerKleenean separated(BasicSetType const& bs) const;
    ValidatedLowerKleenean inside(BasicSetType const& bs) const;
    BoundingSetType bounding_box() const;
};

template<class P, class SIG> using ParametrisedPatchMultifunction = ParametrisedMultifunction<Function<P,SIG>,typename SignatureTraits<SIG>::BoundedDomainType>;
using ValidatedVectorMultivariateParametrisedPatchMultifunction = ParametrisedPatchMultifunction<ValidatedTag,RealVector(RealVector)>;



//! \ingroup Multifunction
//! \brief A set of the form \f$\{ F(x) \mid x\in C \}\f$ where \f$F:\mathbb{X} \to \mathcal{S}(\mathbb{Y})\f$ and \f$C:\mathcal{S}(\mathbb{X})\f$ for some /set monad \f$\mathcal{S}\f$ given by \a SET, where \f$\mathbb{X}\f$ is \a ARG and \f$\mathbb{Y}\f$ is \a RES.
template<class P, class RES, class ARG, template<class,class>class SET> class MultiImageSet;

template<class P, class RES, class ARG> using CompactMultiImageSet = MultiImageSet<P,RES,ARG,CompactSet>;
template<class RES, class ARG> using EffectiveCompactMultiImageSet = MultiImageSet<EffectiveTag,RES,ARG,CompactSet>;
template<class RES, class ARG> using ValidatedCompactMultiImageSet = MultiImageSet<ValidatedTag,RES,ARG,CompactSet>;


/*
template<class P> class MultifunctionPatch<P,RealVector(RealVector)> {
    using ARG=RealVector; using RES=RealVector; using SIG=RES(ARG);
    Function<P,SIG> _f;
    BoxDomainType _pdom;
  public:
    MultifunctionPatch(VectorMultivariateFunction<P> f, BoxDomainType pdom) : _f(f), _pdom(pdom) { }

    SizeType argument_size() const { return this->_f.argument_size()-this->_pdom.dimension(); }
    SizeType result_size() const { return this->_f.result_size(); }
    Function<P,SIG> argument_parameter_function() const { return this->_f; }
    BoxDomainType parameter_domain() const { return this->_pdom; }

    ValidatedImageSet operator() (Vector<ValidatedNumber> const& x) const;

    friend OutputStream& operator<<(OutputStream& os, MultifunctionPatch<P,SIG> const& fm) {
        return os << "MultifunctionPatch(function="<<fm._f<<", parameter_domain="<<fm._pdom<<")"; }
};
using ValidatedVectorMultivariateMultifunctionPatch = MultifunctionPatch<ValidatedTag,VectorMultivariate>;


template<class P> auto
MultifunctionPatch<P,RealVector(RealVector)>::operator() (Vector<ValidatedNumber> const& x) const -> ValidatedImageSet
{
    auto fx=ValidatedFunction<SIG>::constant(parameter_domain().dimension(),x);
    auto fid=ValidatedFunction<SIG>::identity(parameter_domain().dimension());
    Function<P,SIG> fim=compose(_f,join(fx,fid));
    return ValidatedImageSet(this->parameter_domain(),fim);
}

static_assert(AMultifunction<MultifunctionPatch<ValidatedTag,RealVector(RealVector)>,ValidatedTag,RealVector(RealVector),LocatedSet>);

#warning Unparametrised image function
ValidatedConstrainedImageSet image(ValidatedConstrainedImageSet const&, ValidatedVectorMultivariateMultifunctionPatch const&);
*/

template<class PR> class Sweeper;

template<class P, class SIG, class PR> class ParametrisedMultifunctionModel;
template<class P, class PR> using VectorMultivariateParametrisedMultifunctionModel = ParametrisedMultifunctionModel<P,RealVector(RealVector),PR>;
template<class PR> using ValidatedVectorMultivariateParametrisedMultifunctionModel = VectorMultivariateParametrisedMultifunctionModel<ValidatedTag,PR>;

using ValidatedImageSet = ValidatedConstrainedImageSet;

//! \brief A multifunction \f$F:\R^m \mvto \R^n\f$ defined as \f$F(x) = \{ f(x,p) \mid p \in P \}\f$ for some set \f$P\f$, where \f$f\f$ is a function model defined on a bounded domain and using numeric type \f$Float<PR>\f$.
//! \deprecated Probably unnecessary; should be able to use ParametrisedMultifunction
template<class P, class PR> class ParametrisedMultifunctionModel<P,RealVector(RealVector),PR> {
    using ARG=RealVector; using RES=RealVector; using SIG=RES(ARG);
    using FLT=RawFloatType<PR>;
    FunctionModel<P,SIG,PR> _f;
    SizeType _as;
  public:
    ParametrisedMultifunctionModel(BoxDomainType dom, VectorMultivariateFunction<P> f, BoxDomainType params, Sweeper<FLT> swp);
    ParametrisedMultifunctionModel(SizeType as, FunctionModel<P,SIG,PR> const& f);

    SizeType argument_size() const { return this->_as; }
    BoxDomainType domain() const { return project(this->_f.domain(),Range(0u,this->argument_size())); }
    BoxDomainType error_domain() const { return project(this->_f.domain(),Range(this->argument_size(),this->_f.argument_size())); }

    ValidatedImageSet operator() (Vector<ValidatedNumber> const& x) const;

    friend OutputStream& operator<<(OutputStream& os, ParametrisedMultifunctionModel<P,SIG,PR> const& fm) {
        return os << "MultifunctionModel(model="<<fm._f<<", argument_size="<<fm._as<<")"; }
};


template<class P, class PR> auto
ParametrisedMultifunctionModel<P,RealVector(RealVector),PR>::operator() (Vector<ValidatedNumber> const& x) const -> ValidatedImageSet
{
    auto fx=factory(this->_f).create_constants(error_domain(),x);
    auto fid=factory(this->_f).create_identity(error_domain());
    FunctionModel<P,SIG,PR> fim=compose(_f,join(fx,fid));
    return ValidatedImageSet(this->error_domain(),fim);
}

static_assert(AMultifunction<ParametrisedMultifunctionModel<ValidatedTag,RealVector(RealVector),DoublePrecision>,ValidatedTag,RealVector(RealVector),LocatedSet>);


template<class P, class RES, class ARG> class Function<P,CompactSet<P,RES>(ARG)>
    : public Multifunction<P,RES(ARG),CompactSet>
{
    using Multifunction<P,RES(ARG),CompactSet>::Multifunction;
};
template<class P, class RES, class ARG> class Function<P,LocatedSet<P,RES>(ARG)>
    : public Multifunction<P,RES(ARG),LocatedSet>
{
    using Multifunction<P,RES(ARG),LocatedSet>::Multifunction;
};



} // namespace Ariadne



#include "../geometry/set_wrapper.hpp"

namespace Ariadne {

//! \brief A set of class \a SET of functions of type \a SIG, providing information of kind \a P.
//! Hence models the concept \a SET<P,Function<P,SIG>>.
template<class P, class SIG, template<class,class>class SET=LocatedSet> class FunctionSet;

template<class F, class P, class SIG, template<class,class>class SET> struct IsFunctionSet {
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;

    template<class FF, class=decltype(declval<SET<P,RES>>()=declval<FF>()(declval<Argument<Number<P>>>()))>
        static std::true_type test(int);
    template<class FF>
        static std::false_type test(...);
    static const bool value = decltype(test<F>(1))::value;
};

template<class F, class P, class SIG, template<class,class>class SET> concept AFunctionSet = IsFunctionSet<F,P,SIG,SET>::value;


template<class P, class SIG, template<class,class>class SET> class FunctionSetInterface {
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    virtual ~FunctionSetInterface() = default;
    virtual SET<P,RES> _call(Argument<Number<P>> const& x) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

//! \brief A set of functions of type \a SET, signature \a SIG=RES(ARG), and information tag \a P.
//! i.e. \f$F\subset\mathcal{C}(\tpX;\tpY)\f$ for \f$\mathbb{X}\f$ being \a ARG and \f$\mathbb{Y}\f$ being \a RES.
template<class P, class SIG, template<class,class>class SET> class FunctionSet
    : public Handle<FunctionSetInterface<P,SIG,SET>>
{
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    typedef FunctionSetInterface<P,SIG,SET> Interface;
    typedef D DomainType;
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;

    using Handle<Interface>::Handle;

    //! \brief <p/>
    template<AFunctionSet<P,SIG,SET> FS> explicit FunctionSet(FS const& fs);
    //! \brief <p/>
    SET<P,RES> operator() (Argument<Number<P>> const& x) const { return this->reference()._call(x); }
    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os,FunctionSet<P,SIG,SET> const& fs) { return fs.reference()._write(os); }
};

template<class FS, class P, class SIG, template<class,class>class SET> class FunctionSetWrapper
    : public virtual FunctionSetInterface<P,SIG,SET>
{
    using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    using RES=typename SignatureTraits<SIG>::ResultKind;
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;

    FunctionSetWrapper(FS mf) : _fs(mf) { }
    FS const& base() const { return this->_fs; }
    virtual SET<P,RES> _call(Argument<Number<P>> const& x) const final override {
        return static_cast<SET<P,RES>>(this->_fs(x)); }
    virtual OutputStream& _write(OutputStream& os) const  final override {
        return os << this->base(); }
  private:
    FS _fs;
};

template<class P, class SIG,template<class,class>class SET> template<AFunctionSet<P,SIG,SET> FS>
FunctionSet<P,SIG,SET>::FunctionSet(FS const& fs)
    : Handle<Interface>(std::make_shared<FunctionSetWrapper<FS,P,SIG,SET>>(fs)) { }

template<class RES, class ARG> class CompactSet<ValidatedTag,RES(ARG)>
    : public FunctionSet<ValidatedTag,RES(ARG),CompactSet>
{
    using FunctionSet<ValidatedTag,RES(ARG),CompactSet>::FunctionSet;
};
template<class RES, class ARG> class LocatedSet<ValidatedTag,RES(ARG)>
    : public FunctionSet<ValidatedTag,RES(ARG),LocatedSet>
{
    using FunctionSet<ValidatedTag,RES(ARG),LocatedSet>::FunctionSet;
};




template<class A1, class A2> using JoinType = decltype(join(declval<A1>(),declval<A2>()));

//! brief A multivalued function of type \f$F:\mathbb{X}_1 \to (\mathbb{X}_2\to\mathcal{S}(\mathbb{Y}))\f$. where \f$\mathcal{S}\f$ is the \a SET type,
//! \f$\mathbb{X}_1\f$ is \a ARG1, \f$\mathbb{X}_2\f$ is \a ARG2, and \f$\mathbb{Y}\f$ is \a RES.
//! Equivalent to functions of the form \f$\mathbb{X}_1\times\mathbb{X}_2 \to \mathbb{Y}\f$.
template<class P, class RES, class ARG1, class ARG2, template<class,class>class SET> class FunctionSetFunction {
    typedef JoinType<ARG1,ARG2> ARG1S;
  public:
    using D=typename SignatureTraits<RES(ARG1)>::DomainType;
  public:
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    typedef SET<P,RES(ARG1)> FunctionSetType;
    class Interface {
      private: public:
        virtual SET<P,RES(ARG1)> _call(Argument<Number<P>> const&) const = 0;
        virtual Function<P,SET<P,RES>(ARG1S)> _curry() const = 0;
        virtual OutputStream& _write(OutputStream&) const = 0;
    };
    SET<P,RES(ARG1)> operator() (Argument<Number<P>> const& x) const { return this->_handle._ptr->_call(x); }
    Function<P,SET<P,RES>(ARG1S)> curry() const { return this->_handle._ptr->_curry(); }
  private:
    Handle<Interface> _handle;
};


template<class P, class RES, class ARG1, class ARG2> class Function<P,CompactSet<P,RES(ARG2)>(ARG1)>
    : public FunctionSetFunction<P,RES,ARG1,ARG2,CompactSet>
{ };

template<class P, template<class,class>class SET=LocatedSet> using Multiflow = Function<P,SET<P,RealVector(Real)>(RealVector)>;
template<class P> using CompactMultiflow = Function<P,CompactSet<P,RealVector(Real)>(RealVector)>;
template<class P> using LocatedMultiflow = Multiflow<P,LocatedSet>;

using ValidatedCompactMultiflow = Multiflow<ValidatedTag,CompactSet>;



} // namespace Ariadne

#endif

