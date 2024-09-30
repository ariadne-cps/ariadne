/***************************************************************************
 *            function/restricted_function.hpp
 *
 *  Copyright  2024  Pieter Collins
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

/*! \file function/restricted_function.hpp
 *  \brief Restrictions of functions to subdomains
 */

#ifndef ARIADNE_RESTRICTED_FUNCTION_HPP
#define ARIADNE_RESTRICTED_FUNCTION_HPP

#include "../numeric/number.decl.hpp"
#include "../function/function.decl.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"

#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/range.hpp"
#include "../algebra/operations.hpp"
#include "../algebra/algebra_interface.hpp"
#include "../algebra/algebra_mixin.hpp"

#include "../function/domain.hpp"
#include "../function/function_interface.hpp"
#include "../function/function_patch_interface.hpp"
#include "../function/function_patch_mixin.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function.hpp"

#include "../function/function_mixin.tpl.hpp"


namespace Ariadne {

template<class P, class SIG> class RestrictedFunction;

using ApproximateScalarUnivariateRestrictedFunction = RestrictedFunction<ApproximateTag,RealScalar(RealScalar)>;
using ApproximateScalarMultivariateRestrictedFunction = RestrictedFunction<ApproximateTag,RealScalar(RealVector)>;
using ApproximateVectorUnivariateRestrictedFunction = RestrictedFunction<ApproximateTag,RealVector(RealScalar)>;
using ApproximateVectorMultivariateRestrictedFunction = RestrictedFunction<ApproximateTag,RealVector(RealVector)>;
using ValidatedScalarUnivariateRestrictedFunction = RestrictedFunction<ValidatedTag,RealScalar(RealScalar)>;
using ValidatedScalarMultivariateRestrictedFunction = RestrictedFunction<ValidatedTag,RealScalar(RealVector)>;
using ValidatedVectorUnivariateRestrictedFunction = RestrictedFunction<ValidatedTag,RealVector(RealScalar)>;
using ValidatedVectorMultivariateRestrictedFunction = RestrictedFunction<ValidatedTag,RealVector(RealVector)>;

template<class P, class RES, class ARG> FunctionPatchCreator<FunctionPatchFactory<P>,ARG> factory(FunctionPatchMixin<RestrictedFunction<P,RES(ARG)>,P,RES(ARG)> const& rf) {
    return factory(static_cast<RestrictedFunction<P,RES(ARG)>const&>(rf)); }

template<class P, class... ARGS> struct Element<RestrictedFunction<P,RealVector(ARGS...)>> {
    typedef RestrictedFunction<P,RealScalar(ARGS...)> Type; };

template<class P> class RestrictedFunctionFactory;

template<class P, class... ARGS> class RestrictedFunction<P,RealScalar(ARGS...)>
    : public FunctionPatchMixin<RestrictedFunction<P,RealScalar(ARGS...)>,P,RealScalar(ARGS...)>
{
    using RES=RealScalar;
    using SIG=RES(ARGS...);
  public:
    typedef FunctionPatchInterface<P,SIG> Interface;
    using typename FunctionPatchInterface<P,SIG>::ArgumentSizeType;
    using typename FunctionPatchInterface<P,SIG>::ArgumentIndexType;
    using typename FunctionPatchInterface<P,SIG>::ResultSizeType;
    using typename FunctionPatchInterface<P,SIG>::ResultIndexType;
    typedef P Paradigm;
    typedef typename SignatureTraits<SIG>::BoundedDomainType DomainType;
    typedef typename SignatureTraits<SIG>::BoundedCodomainType CodomainType;
    typedef typename SignatureTraits<SIG>::BoundedRangeType RangeType;
    template<class X> using Result = typename SignatureTraits<SIG>::template Result<X>;
    template<class X> using Argument = typename SignatureTraits<SIG>::template Argument<X>;
  private:
    Function<P,SIG> _f;
    DomainType _dom;
  public:
    explicit RestrictedFunction(Function<P,SIG> f, DomainType dom) : _f(f), _dom(dom) {
        ARIADNE_PRECONDITION(f.argument_size()==dom.dimension()); }

    RestrictedFunction<P,SIG> create_zero() const {
        return RestrictedFunction(Function<P,SIG>::zero(this->argument_size()),this->domain()); }
    RestrictedFunction<P,SIG> create_constant(Number<P> c) const {
        return RestrictedFunction(Function<P,SIG>::constant(this->argument_size(),c),this->domain()); }

    explicit operator FunctionPatch<P,SIG> () const {
        return FunctionPatch<P,SIG>(new RestrictedFunction<P,SIG>(*this)); }
    friend Function<P,SIG> const& cast_unrestricted(RestrictedFunction<P,SIG> const& rf) {
        return rf._f; }

    ArgumentSizeType argument_size() const override {
        return this->_f.argument_size(); }
    ResultSizeType result_size() const override {
        return this->_f.result_size(); }

    const DomainType domain() const override {
        return this->_dom; }
    const CodomainType codomain() const override {
        return cast_exact_interval(this->range()); }
    const RangeType range() const {
        return apply(this->_f,BoxRangeType(this->_dom)); }

    Void clobber() const {
        ARIADNE_NOT_IMPLEMENTED; }
    ValidatedErrorNumber error() const {
        ARIADNE_NOT_IMPLEMENTED; }
    ValidatedErrorNumber norm() const {
        auto rng=this->range(); return cast_positive(ValidatedUpperNumber(max(-rng.lower_bound(),rng.upper_bound()))); }
    friend ValidatedErrorNumber norm(RestrictedFunction<P,SIG> const& rf) {
        return rf.norm(); }

    template<class XX> Result<XX> operator() (Argument<XX> const& x) const {
        check_domain(this->_dom,x); return this->_f(x); }

    friend RestrictedFunction<P,SIG> partial_evaluate(RestrictedFunction<P,SIG> const& rf, SizeType j, Number<P> x) {
        ARIADNE_NOT_IMPLEMENTED; }
    template<class XX> friend Result<XX> unchecked_evaluate(RestrictedFunction<P,SIG> const& rf, Argument<XX> const& x) {
        return rf._f(x); }

    friend RestrictedFunction<P,RealScalar(ARGS...)> compose(Function<P,RealScalar(RealScalar)> g, RestrictedFunction<P,RealScalar(ARGS...)> rf) {
        return RestrictedFunction<P,SIG>::_compose(g,rf); }
    friend RestrictedFunction<P,RealVector(ARGS...)> compose(Function<P,RealVector(RealScalar)> g, RestrictedFunction<P,RealScalar(ARGS...)> rf) {
        return RestrictedFunction<P,SIG>::_compose(g,rf); }

    friend RestrictedFunction<P,SIG> restriction(RestrictedFunction<P,SIG> rf, BoxDomainType dom) {
        ARIADNE_PRECONDITION(subset(dom,rf._dom)); return RestrictedFunction<P,SIG>(rf._f,dom); }
    friend RestrictedFunction<P,SIG> embed(BoxDomainType, RestrictedFunction<P,SIG>, BoxDomainType) {
        ARIADNE_NOT_IMPLEMENTED; }

    friend RestrictedFunction<P,SIG> nul(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Nul(),rf); }
    friend RestrictedFunction<P,SIG> pos(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Pos(),rf); }
    friend RestrictedFunction<P,SIG> neg(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Neg(),rf); }
    friend RestrictedFunction<P,SIG> sqr(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Sqr(),rf); }
    friend RestrictedFunction<P,SIG> hlf(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Hlf(),rf); }
    friend RestrictedFunction<P,SIG> rec(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Rec(),rf); }
    friend RestrictedFunction<P,SIG> add(RestrictedFunction<P,SIG> rf1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Add(),rf1,rf2); }
    friend RestrictedFunction<P,SIG> sub(RestrictedFunction<P,SIG> rf1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Sub(),rf1,rf2); }
    friend RestrictedFunction<P,SIG> mul(RestrictedFunction<P,SIG> rf1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Mul(),rf1,rf2); }
    friend RestrictedFunction<P,SIG> div(RestrictedFunction<P,SIG> rf1, RestrictedFunction<P,SIG> rf2) {
         return RestrictedFunction<P,SIG>::_apply(Div(),rf1,rf2); }
    friend RestrictedFunction<P,SIG> add(RestrictedFunction<P,SIG> rf1, Number<P> c2) {
        return RestrictedFunction<P,SIG>::_apply(Add(),rf1,c2); }
    friend RestrictedFunction<P,SIG> sub(RestrictedFunction<P,SIG> rf1, Number<P> c2) {
        return RestrictedFunction<P,SIG>::_apply(Sub(),rf1,c2); }
    friend RestrictedFunction<P,SIG> mul(RestrictedFunction<P,SIG> rf1, Number<P> c2) {
        return RestrictedFunction<P,SIG>::_apply(Mul(),rf1,c2); }
    friend RestrictedFunction<P,SIG> div(RestrictedFunction<P,SIG> rf1, Number<P> c2) {
        return RestrictedFunction<P,SIG>::_apply(Div(),rf1,c2); }
    friend RestrictedFunction<P,SIG> add(Number<P> c1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Add(),c1,rf2); }
    friend RestrictedFunction<P,SIG> sub(Number<P> c1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Sub(),c1,rf2); }
    friend RestrictedFunction<P,SIG> mul(Number<P> c1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Mul(),c1,rf2); }
    friend RestrictedFunction<P,SIG> div(Number<P> c1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Div(),c1,rf2); }
    friend RestrictedFunction<P,SIG> pow(RestrictedFunction<P,SIG> rf1, Int n2) {
        return RestrictedFunction<P,SIG>::_apply(Pow(),rf1,n2); }
    friend RestrictedFunction<P,SIG> sqrt(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Sqrt(),rf); }
    friend RestrictedFunction<P,SIG> exp(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Exp(),rf); }
    friend RestrictedFunction<P,SIG> log(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Log(),rf); }
    friend RestrictedFunction<P,SIG> sin(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Sin(),rf); }
    friend RestrictedFunction<P,SIG> cos(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Cos(),rf); }
    friend RestrictedFunction<P,SIG> tan(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Tan(),rf); }
    friend RestrictedFunction<P,SIG> asin(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Asin(),rf); }
    friend RestrictedFunction<P,SIG> acos(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Acos(),rf); }
    friend RestrictedFunction<P,SIG> atan(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Atan(),rf); }

    friend RestrictedFunction<P,SIG> operator+(RestrictedFunction<P,SIG>, RestrictedFunction<P,SIG>);
    friend RestrictedFunction<P,SIG> operator-(RestrictedFunction<P,SIG>, RestrictedFunction<P,SIG>);
    friend RestrictedFunction<P,SIG> operator*(RestrictedFunction<P,SIG>, RestrictedFunction<P,SIG>);
    friend RestrictedFunction<P,SIG> operator/(RestrictedFunction<P,SIG>, RestrictedFunction<P,SIG>);
    friend RestrictedFunction<P,SIG> operator+(RestrictedFunction<P,SIG>, Number<P>);
    friend RestrictedFunction<P,SIG> operator-(RestrictedFunction<P,SIG>, Number<P>);
    friend RestrictedFunction<P,SIG> operator*(RestrictedFunction<P,SIG>, Number<P>);
    friend RestrictedFunction<P,SIG> operator/(RestrictedFunction<P,SIG>, Number<P>);
    friend RestrictedFunction<P,SIG> operator+(Number<P>, RestrictedFunction<P,SIG>);
    friend RestrictedFunction<P,SIG> operator-(Number<P>, RestrictedFunction<P,SIG>);
    friend RestrictedFunction<P,SIG> operator*(Number<P>, RestrictedFunction<P,SIG>);
    friend RestrictedFunction<P,SIG> operator/(Number<P>, RestrictedFunction<P,SIG>);

    friend RestrictedFunction<P,SIG> abs(RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>::_apply(Abs(),rf); }
    friend RestrictedFunction<P,SIG> max(RestrictedFunction<P,SIG> rf1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Max(),rf1,rf2); }
    friend RestrictedFunction<P,SIG> min(RestrictedFunction<P,SIG> rf1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Min(),rf1,rf2); }
    friend RestrictedFunction<P,SIG> max(RestrictedFunction<P,SIG> rf1, Number<P> c2) {
        return RestrictedFunction<P,SIG>::_apply(Max(),rf1,c2); }
    friend RestrictedFunction<P,SIG> min(RestrictedFunction<P,SIG> rf1, Number<P> c2) {
        return RestrictedFunction<P,SIG>::_apply(Min(),rf1,c2); }
    friend RestrictedFunction<P,SIG> max(Number<P> c1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Max(),c1,rf2); }
    friend RestrictedFunction<P,SIG> min(Number<P> c1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>::_apply(Min(),c1,rf2); }

    friend RestrictedFunction<P,SIG> derivative(RestrictedFunction<P,SIG> const& rf, ArgumentIndexType k) {
        return RestrictedFunction<P,SIG>(derivative(rf._f,k),rf._dom); }
    friend RestrictedFunction<P,SIG> antiderivative(RestrictedFunction<P,SIG> const& rf, ArgumentIndexType k) {
        return antiderivative(rf,k,0); }
    friend RestrictedFunction<P,SIG> antiderivative(RestrictedFunction<P,SIG> const& rf, ArgumentIndexType k, Number<P> c) {
        ARIADNE_NOT_IMPLEMENTED; }

    friend Bool refines(RestrictedFunction<ValidatedTag,SIG>, RestrictedFunction<ValidatedTag,SIG>) {
        ARIADNE_NOT_IMPLEMENTED; }

    friend FunctionPatchCreator<RestrictedFunctionFactory<P>,ARGS...> factory(RestrictedFunction<P,SIG> const& rf) {
        return FunctionPatchCreator<RestrictedFunctionFactory<P>,ARGS...>(rf.domain(),RestrictedFunctionFactory<P>()); }
    friend OutputStream& operator<<(OutputStream& os, RestrictedFunction<P,SIG> const& rf) {
        return os << "RestrictedFunction(" << rf._f << ", " << rf._dom << ")"; }
  private:
    template<class OP> static RestrictedFunction<P,SIG> _apply(OP op, RestrictedFunction<P,SIG> rf) {
        return RestrictedFunction<P,SIG>(op(rf._f),rf._dom); }
    template<class OP> static RestrictedFunction<P,SIG> _apply(OP op, RestrictedFunction<P,SIG> rf1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>(op(rf1._f,rf2._f),intersection(rf1._dom,rf2._dom)); }
    template<class OP> static RestrictedFunction<P,SIG> _apply(OP op, RestrictedFunction<P,SIG> rf1, Number<P> c2) {
        return RestrictedFunction<P,SIG>(op(rf1._f,c2),rf1._dom); }
    template<class OP> static RestrictedFunction<P,SIG> _apply(OP op, RestrictedFunction<P,SIG> rf, Int n) {
        return RestrictedFunction<P,SIG>(op(rf._f,n),rf._dom); }
    template<class OP> static RestrictedFunction<P,SIG> _apply(OP op, Number<P> c1, RestrictedFunction<P,SIG> rf2) {
        return RestrictedFunction<P,SIG>(op(c1,rf2._f),rf2._dom); }

    template<class R> static RestrictedFunction<P,R(ARGS...)> _compose(Function<P,R(RealScalar)> const& g, RestrictedFunction<P,RealScalar(ARGS...)> rf) {
        return RestrictedFunction<P,R(ARGS...)>(compose(g,rf._f),rf._dom); }

    template<class XX> Void _check_domain(Argument<XX> const& x) const {
        return check_domain(this->domain(),x); }
  private:
    virtual FunctionPatchInterface<P,SIG>* _create() const override {
        return new RestrictedFunction<P,SIG>(*this); }
};

template class RestrictedFunction<ValidatedTag,RealScalar(RealVector)>;

template<class VF> struct Element {
    using I=SizeType;
    using SF=std::remove_cv_t<decltype(declval<VF>().get(declval<I>()))>;
    VF& _vf; I _i;
    Element(VF& vf, I i) : _vf(vf), _i(i) { }
    Element<VF>& operator=(SF const& sf) { _vf.set(_i,sf); return *this; }
    operator SF const& () const { return _vf.get(_i); }
};

template<class P, class... ARGS> class RestrictedFunction<P,RealVector(ARGS...)>
    : public FunctionPatchMixin<RestrictedFunction<P,RealVector(ARGS...)>,P,RealVector(ARGS...)>
{
    using RES=RealVector;
    using SIG=RES(ARGS...);
  public:
    typedef FunctionPatchInterface<P,SIG> Interface;
    typedef typename FunctionPatchInterface<P,SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename FunctionPatchInterface<P,SIG>::ArgumentIndexType ArgumentIndexType;
    typedef typename FunctionPatchInterface<P,SIG>::ResultSizeType ResultSizeType;
    typedef typename FunctionPatchInterface<P,SIG>::ResultIndexType ResultIndexType;
/*
    using typename FunctionPatchInterface<P,SIG>::ArgumentSizeType;
    using typename FunctionPatchInterface<P,SIG>::ArgumentIndexType;
    using typename FunctionPatchInterface<P,SIG>::ResultSizeType;
    using typename FunctionPatchInterface<P,SIG>::ResultIndexType;
*/
    typedef P Paradigm;
    typedef typename SignatureTraits<SIG>::BoundedDomainType DomainType;
    typedef typename SignatureTraits<SIG>::BoundedCodomainType CodomainType;
    typedef typename SignatureTraits<SIG>::BoundedRangeType RangeType;
    typedef PositiveUpperNumber<P> NormType;
    template<class X> using Result = typename SignatureTraits<SIG>::template Result<X>;
    template<class X> using Argument = typename SignatureTraits<SIG>::template Argument<X>;
  private:
    Function<P,SIG> _f;
    DomainType _dom;
  public:
    explicit RestrictedFunction(Function<P,SIG> f, DomainType dom) : _f(f), _dom(dom) {
        ARIADNE_PRECONDITION(f.argument_size()==dom.dimension()); }

//    RestrictedFunction<P,SIG> create_zero() const;
//    RestrictedFunction<P,SIG> create_constant(Number<P>) const;

    explicit operator FunctionPatch<P,SIG> () const {
        return FunctionPatch<P,SIG>(new RestrictedFunction<P,SIG>(*this)); }
    friend Function<P,SIG> const& cast_unrestricted(RestrictedFunction<P,SIG> const& rf) {
        return rf._f; }

    ArgumentSizeType argument_size() const {
        return this->_f.argument_size(); };
    ResultSizeType result_size() const {
        return this->_f.result_size(); };
    const DomainType domain() const {
        return this->_dom; }
    const CodomainType codomain() const {
        return cast_exact_box(this->range()); }
    const RangeType range() const {
        return apply(this->_f,BoxRangeType(this->_dom)); }

    NormType norm() const {
        UpperNumber<P> nrm(0u); RangeType rng=this->range();
        for (SizeType i=0; i!=rng.dimension(); ++i) {
            nrm=max(nrm,max(-rng[i].lower_bound(),+rng[i].upper_bound())); }
        return cast_positive(nrm); }
    friend NormType norm(RestrictedFunction<P,SIG> const& rf) {
        return rf.norm(); }

    Vector<ValidatedErrorNumber> errors() const {
        ARIADNE_NOT_IMPLEMENTED; }
    ValidatedErrorNumber error() const {
        Vector<ValidatedErrorNumber> es=this->errors(); ValidatedErrorNumber err=0u;
        for (SizeType i=0; i!=es.size(); ++i) {err=max(err,es[i]); } return err; }

    Void clobber() {
        ARIADNE_NOT_IMPLEMENTED; }

    friend Bool refines(RestrictedFunction<P,SIG> const& rf1, RestrictedFunction<P,SIG> const& rf2) {
        ARIADNE_NOT_IMPLEMENTED; };

    friend RestrictedFunction<P,SIG> restriction(RestrictedFunction<P,SIG> const& rf, DomainType const& dom) {
        ARIADNE_PRECONDITION(subset(dom,rf.domain())); return RestrictedFunction(rf._f,dom); }

    template<class XX> Result<XX> operator() (Argument<XX> const& x) const {
        this->_check_domain(x); return this->_f(x); }

    template<class XX> Result<XX> _unchecked_call(Argument<XX> const& x) const {
        return this->_f(x); }

    template<class XX> friend Result<XX> evaluate(RestrictedFunction<P,SIG> const& rf, Argument<XX> const& x) {
        return rf.operator()(x); }
    template<class XX> friend Result<XX> unchecked_evaluate(RestrictedFunction<P,SIG> const& rf, Argument<XX> const& x) {
        return rf._unchecked_call(x); }

    friend RestrictedFunction<P,SIG> partial_evaluate(RestrictedFunction<P,SIG> const& rf, SizeType j, Number<P> const& x) {
        ARIADNE_NOT_IMPLEMENTED; }

//    Element<const RestrictedFunction<P,SIG>> operator[](SizeType i) const;
    RestrictedFunction<P,RealScalar(ARGS...)> const operator[](ResultIndexType i) const {
        return this->get(i); }
    Element<RestrictedFunction<P,SIG>> operator[](ResultIndexType i) {
        return Element<RestrictedFunction<P,SIG>>(*this,i); }
    RestrictedFunction<P,RealScalar(ARGS...)> get(ResultIndexType i) const {
        return RestrictedFunction<P,RealScalar(ARGS...)>(this->_f[i],this->_dom); }
    Void set(ResultIndexType i, Function<P,RealScalar(ARGS...)> sf) {
        this -> _f[i]=sf; }

    friend RestrictedFunction<P,RealScalar(ARGS...)> compose(Function<P,RealScalar(RealVector)> g, RestrictedFunction<P,RealVector(ARGS...)> rf) {
        return RestrictedFunction<P,SIG>::_compose(g,rf); }
    friend RestrictedFunction<P,RealVector(ARGS...)> compose(Function<P,RealVector(RealVector)> g, RestrictedFunction<P,RealVector(ARGS...)> rf) {
        return RestrictedFunction<P,SIG>::_compose(g,rf); }

    RestrictedFunction<P,RealVector(ARGS...)>& adjoin(RestrictedFunction<P,RealScalar(ARGS...)> const& rf) {
        ARIADNE_PRECONDITION(this->domain()==rf.domain());
        this->_f=join(this->_f,cast_unrestricted(rf)); return *this; }
    friend RestrictedFunction<P,RealVector(ARGS...)> join(RestrictedFunction<P,RealScalar(ARGS...)> rf1,  RestrictedFunction<P,RealVector(ARGS...)> rf2) {
        ARIADNE_PRECONDITION(rf1.domain()==rf2.domain());
        return RestrictedFunction<P,RealVector(ARGS...)>(join(rf1._f,rf2._f),rf1._dom); }
    friend RestrictedFunction<P,RealVector(ARGS...)> join(RestrictedFunction<P,RealVector(ARGS...)> rf1,  RestrictedFunction<P,RealScalar(ARGS...)> rf2) {
        ARIADNE_PRECONDITION(rf1.domain()==rf2.domain());
        return RestrictedFunction<P,RealVector(ARGS...)>(join(rf1._f,rf2._f),rf1._dom); }
    friend RestrictedFunction<P,RealVector(ARGS...)> join(RestrictedFunction<P,RealVector(ARGS...)> rf1,  RestrictedFunction<P,RealVector(ARGS...)> rf2) {
        ARIADNE_PRECONDITION(rf1.domain()==rf2.domain());
        return RestrictedFunction<P,RealVector(ARGS...)>(join(rf1._f,rf2._f),rf1._dom); }
    friend RestrictedFunction<P,RealVector(ARGS...)> combine(RestrictedFunction<P,RealVector(ARGS...)> rf1,  RestrictedFunction<P,RealVector(ARGS...)> rf2) {
        return RestrictedFunction<P,RealVector(ARGS...)>(combine(rf1._f,rf2._f),product(rf1._dom,rf2._dom)); }

    friend RestrictedFunction<P,SIG> derivative(RestrictedFunction<P,SIG> const& f, ArgumentIndexType j) {
        return RestrictedFunction<P,SIG>(derivative(f._f,j),f._dom); }

    friend RestrictedFunction<P,SIG> embed(BoxDomainType dom1, RestrictedFunction<P,SIG> rf2, BoxDomainType dom3) {
        ARIADNE_NOT_IMPLEMENTED; }
        //return RestrictedFunction(product(dom1,rf2.domain(),dom3),embed(dom1.dimension(),rf2._f,dom3.dimension())); }

    friend FunctionPatchCreator<RestrictedFunctionFactory<P>,ARGS...> factory(RestrictedFunction<P,SIG> const& rf);
    friend OutputStream& operator<<(OutputStream& os, RestrictedFunction<P,SIG> const& rf) {
        return os << "RestrictedFunction(" << rf._f << ", " << rf._dom << ")"; }
  private:
    template<class R> static RestrictedFunction<P,R(ARGS...)> _compose(Function<P,R(RealVector)> const& g, RestrictedFunction<P,RealVector(ARGS...)> rf) {
        return RestrictedFunction<P,R(ARGS...)>(compose(g,rf._f),rf._dom); }

    template<class XX> Void _check_domain(Argument<XX> const& x) const {
        return check_domain(this->domain(),x); }
  private:
    virtual ScalarFunctionPatchInterface<P,ARGS...>* _get(ResultIndexType i) const {
        return new RestrictedFunction<P,RealScalar(ARGS...)>(this->get(i)); }
    virtual FunctionPatchInterface<P,SIG>* _create() const {
        return new RestrictedFunction<P,SIG>(*this); }
};



template<class P> class RestrictedFunctionFactory
    : public FunctionPatchFactoryMixin<RestrictedFunctionFactory<P>,P>
{
  public:
    typedef P Paradigm;
    operator FunctionPatchFactory<P> () const { return FunctionPatchFactory<P>(new RestrictedFunctionFactory<P>(*this)); }
    RestrictedFunctionFactory() { }
    template<class SIG> RestrictedFunctionFactory(RestrictedFunction<P,SIG> const&) { }
    RestrictedFunction<P,RealScalar(RealVector)> create(BoxDomainType const& dom, Function<P,RealScalar(RealVector)> const& sf) const {
        return RestrictedFunction<P,RealScalar(RealVector)>(sf,dom); }
    RestrictedFunction<P,RealVector(RealVector)> create(BoxDomainType const& dom, Function<P,RealVector(RealVector)> const& vf) const {
        return RestrictedFunction<P,RealVector(RealVector)>(vf,dom); }

    RestrictedFunction<P,RealScalar(RealVector)> create_zero(BoxDomainType const& dom) const {
        return this->create(dom,Function<P,RealScalar(RealVector)>::zero(dom.dimension())); }
    RestrictedFunction<P,RealScalar(RealVector)> create_constant(BoxDomainType const& dom, Number<P> const& c) const {
        return this->create(dom,Function<P,RealScalar(RealVector)>::constant(dom.dimension(),c)); }
    RestrictedFunction<P,RealScalar(RealVector)> create_coordinate(BoxDomainType const& dom, SizeType i) const {
        return this->create(dom,Function<P,RealScalar(RealVector)>::coordinate(dom.dimension(),i)); }

    RestrictedFunction<P,RealVector(RealVector)> create_zeros(SizeType rs, BoxDomainType const& dom) const {
        return this->create(dom,Function<P,RealVector(RealVector)>::zeros(rs,dom.dimension())); }
    RestrictedFunction<P,RealVector(RealVector)> create_constants(BoxDomainType const& dom, Vector<Number<P>> c) const {
        return this->create(dom,Function<P,RealVector(RealVector)>::constant(dom.dimension(),c)); }
    RestrictedFunction<P,RealVector(RealVector)> create_projection(BoxDomainType const& dom, Range rng) const {
        ARIADNE_NOT_IMPLEMENTED; }
    RestrictedFunction<P,RealVector(RealVector)> create_identity(BoxDomainType const& dom) const {
        return this->create(dom,Function<P,RealVector(RealVector)>::identity(dom.dimension())); }
};

FunctionPatchCreator<RestrictedFunctionFactory<ValidatedTag>,RealVector> factory(RestrictedFunction<ValidatedTag,RealVector(RealVector)> const& rf) {
    return FunctionPatchCreator<RestrictedFunctionFactory<ValidatedTag>,RealVector>(rf.domain(),RestrictedFunctionFactory<ValidatedTag>()); }

template<> inline String class_name<RestrictedFunction<ValidatedTag,RealVector(RealVector)>>() {
    return "ValidatedVectorMultivariateRestrictedFunction"; }

} // namespace Ariadne

#endif
