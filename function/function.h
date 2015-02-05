/***************************************************************************
 *            function.h
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

/*! \file function.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/declarations.h"
#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/container.h"
#include "utility/metaprogramming.h"

#include "function/function_interface.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/differential.h"

#include "geometry/interval.h"
#include "geometry/box.h"

namespace Ariadne {

template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;

template<class P, class D=BoxDomain> struct VectorFunctionElementReference;

class RealDomain : public IntervalDomain {
  public:
    RealDomain() : IntervalDomain(-inf,+inf) { }
};

class EuclideanDomain : public BoxDomain {
  public:
    EuclideanDomain(SizeType n) : BoxDomain(n,RealDomain()) { }
};


//! \ingroup FunctionModule
//! \brief A generic scalar function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X\f$.
template<class P, class D, class C>
class Function {
    static_assert(IsStronger<P,ApproximateTag>::value,"P must be an information level/paradigm.");
    typedef CanonicalNumberType<P> Y;
  protected:
    std::shared_ptr< const FunctionInterface<P,D,C> > _ptr;
  public:
    typedef P InformationTag;
    typedef P Paradigm;
    typedef Y NumericType;
    typedef D DomainType;
    typedef C CodomainType;

    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = typename ElementTraits<C>::template Type<Y>;

    explicit Function(DomainType dom);
    explicit Function(DomainType dom, Result<Formula<Y>>const& e);
    explicit Function(SizeType as, Result<Formula<Y>>const& e);
    explicit Function(SizeType as, List<Formula<Y>>const& e);
    explicit Function(SizeType rs, DomainType dom);
    explicit Function(SizeType rs, ScalarFunction<P,D> sf);

    Function(InitializerList<ScalarFunction<P,D>> const& lsf);
    Function(List<ScalarFunction<P,D>> const& lsf);
    Function(Vector<ScalarFunction<P,D>> const& lsf);

    static ScalarFunction<P> zero(SizeType as);
    static ScalarFunction<P> constant(SizeType as, NumericType c);
    static ScalarFunction<P> coordinate(SizeType as, SizeType j);

    static VectorFunction<P> zeros(SizeType rs, SizeType as);
    static VectorFunction<P> identity(SizeType n);

    ScalarFunction<P> create_zero() const { return ScalarFunction<P>::zero(this->argument_size()); }
    ScalarFunction<P> create_constant(NumericType c) const { return ScalarFunction<P>::constant(this->argument_size(),c); }
    ScalarFunction<P> create_coordinate(SizeType j) const { return ScalarFunction<P>::coordinate(this->argument_size(),j); }

    Function();
    explicit Function(FunctionInterface<P,D,C>* p) : _ptr(p) { }
    explicit Function(SharedPointer<FunctionInterface<P,D,C>> p) : _ptr(p) { }
    Function(const FunctionInterface<P,D,C>& t) : _ptr(t._clone()) { }
    Function<P,D,C>& operator=(const FunctionInterface<P,D,C>& f) {
        _ptr=std::shared_ptr< FunctionInterface<P,D,C> >(f._clone()); return *this; }

    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
    Function(const Function<PP,D,C>& f)
        : _ptr(std::dynamic_pointer_cast< const FunctionInterface<P,D,C> >(f.managed_pointer())) { }
    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
        Function<P,D,C>& operator=(Result<NumericType> const& c); // { return *this=this->create_constant(c); }

    shared_ptr< const FunctionInterface<P,D,C> > managed_pointer() const  { return _ptr; }
    const FunctionInterface<P,D,C>* raw_pointer() const  { return _ptr.operator->(); }
    const FunctionInterface<P,D,C>& reference() const  { return _ptr.operator*(); }
    operator const FunctionInterface<P,D,C>& () const { return _ptr.operator*(); }

    DomainType domain() const {
        return this->reference().domain(); }
    CodomainType codomain() const {
        return this->reference().codomain(); }
    SizeType argument_size() const {
        return this->reference().argument_size(); }
    SizeType result_size() const {
        return this->reference().result_size(); }
    template<class X> auto operator() (const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }
    template<class X> auto evaluate(const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }
    template<class X> friend auto evaluate(const Function<P,D,C>& f, const Argument<X>& x) -> decltype(f(x)) {
        return f(x); }

    Function<P,D,C> derivative(SizeType k) const {
        return Function<P,D,C>(this->reference()._derivative(k)); }
    friend Function<P,D,C> derivative(Function<P,D,C> const& f, SizeType k) {
        return Function<P,D,C>(f.reference()._derivative(k)); }

    template<class X> Result<Differential<ArithmeticType<X,Y>>> differential(const Argument<X>& x, DegreeType d) const {
        return this->_ptr->_evaluate(Differential<X>::variables(x.size(),x.size(),d,x)); }

    Void set(SizeType i, ScalarFunction<P,D>);
    Function<P,D,IntervalDomain> get(SizeType i) const;
    Function<P,D,IntervalDomain> operator[](SizeType i) const;
    VectorFunctionElementReference<P,D> operator[](SizeType i);
    template<class X> Matrix<ArithmeticType<X,Y>> jacobian(const Argument<X>& x) const;

    OutputStream& write(OutputStream& os) const { return this->_ptr->write(os); }

    friend OutputStream& operator<<(OutputStream& os, Function<P,D,C> const& f) { return f._ptr->write(os); }
};

template<class P, class D, class C> inline OutputStream&
operator<<(OutputStream& os, const Function<P,D,C>& f) {
    return f.write(os); }

template<class P, class D, class C, class X> inline auto
evaluate(const Function<P,D,C>& f, const ElementType<D,X>& x) -> decltype(f(x)) {
    return f(x); }


template<class P, class C, class X> inline auto
differential(const Function<P,IntervalDomain,C>& f, const X& x, DegreeType d)
    -> ElementType<C,Differential<ArithmeticType<Number<P>,X>>> {
    auto dx=Differential<X>::variable(1u,d,x,0u); return f(dx);
}

template<class P, class C, class X> inline auto
differential(const Function<P,BoxDomain,C>& f, const Vector<X>& x, DegreeType d)
    -> ElementType<C,Differential<ArithmeticType<Number<P>,X>>> {
    auto dx=Differential<X>::variables(d,x); return f(dx);
}

/*
template<class P, class D, class C, class X> inline auto
differential(const Function<P,D,C>& f, const ElementType<D,X>& x, DegreeType d)
    -> ElementType<C,Differential<ArithmeticType<Number<P>,X>>> {
    auto dx=Differential<X>::create_identity(x,d); return f(dx);
}
*/


template<class P, class X> ArithmeticType<Number<P>,X>
derivative(const ScalarUnivariateFunction<P>& f, const X& x) {
    return differential(f,x,1u).gradient(); }

template<class P, class X> Vector<ArithmeticType<Number<P>,X>>
tangent(const VectorUnivariateFunction<P>& f, const X& x) {
    return column(differential(f,x,1u).jacobian(),0u); }

template<class P, class X> Covector<ArithmeticType<Number<P>,X>>
gradient(const ScalarMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).gradient(); }

template<class P, class X> Matrix<ArithmeticType<Number<P>,X>>
jacobian(const VectorMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).jacobian(); }

/*
template<class P, class D> class ScalarFunction : public Function<P,D,IntervalDomain> {
    typedef IntervalDomain C;
    typedef CanonicalNumberType<P> X;
  public:
    typedef D DomainType;
    using Function<P,D,C>::Function;
    ScalarFunction<P,D>() : ScalarFunction(D()) { }
    ScalarFunction<P,D>(Function<P,D,C>const& f) : Function<P,D,C>(f) { }
    explicit ScalarFunction<P,D>(SizeType as) : ScalarFunction(D(as)) { }
    ScalarFunction<P,D>(SizeType as, Formula<X> const& e) : ScalarFunction(D(as),e) { }
    explicit ScalarFunction<P,D>(D dom);
    ScalarFunction<P,D>(D dom, Formula<X> const& e);
};

template<class P, class D> class VectorFunction : public Function<P,D,BoxDomain> {
    typedef BoxDomain C;
    typedef CanonicalNumberType<P> X;
  public:
    using Function<P,D,BoxDomain>::Function;
    VectorFunction();
    VectorFunction<P,D>(Function<P,D,C>const& f) : Function<P,D,C>(f) { }
    VectorFunction(SizeType rs, SizeType as);
    VectorFunction(InitializerList<ScalarFunction<P,D>> const& lsf);
    VectorFunction(List<ScalarFunction<P,D>> const& lsf);
    VectorFunction(Vector<ScalarFunction<P,D>> const& vsf);
    VectorFunction(SizeType rs, ScalarFunction<P,D> const& sf);
    VectorFunction(SizeType as, List<Formula<X>> const& le);

    Void set(SizeType i, ScalarFunction<P,D> const& f);
    ScalarFunction<P,D> get(SizeType i) const;

    VectorFunctionElementReference<P,D> operator[](SizeType i);
    ScalarFunction<P,D> operator[](SizeType i) const;

    template<class XX> Matrix<XX> jacobian(const Vector<XX>& x) const;
};
*/

/*
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator*(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator/(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator*(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator/(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator+(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator*(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator/(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator*(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator/(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator+(double, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(double, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator*(double, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator/(double, const ScalarFunction<X>&);

template<class X> ScalarFunction<X> pow(const ScalarFunction<X>&, Int);
template<class X> ScalarFunction<X> neg(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> rec(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> sqr(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> sqrt(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> exp(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> log(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> sin(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> cos(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> tan(const ScalarFunction<X>&);
*/

EffectiveScalarFunction operator+(const EffectiveScalarFunction&);
EffectiveScalarFunction operator-(const EffectiveScalarFunction&);
EffectiveScalarFunction operator+(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator-(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator*(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator/(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator+(const EffectiveScalarFunction&, const EffectiveNumber&);
EffectiveScalarFunction operator-(const EffectiveScalarFunction&, const EffectiveNumber&);
EffectiveScalarFunction operator*(const EffectiveScalarFunction&, const EffectiveNumber&);
EffectiveScalarFunction operator/(const EffectiveScalarFunction&, const EffectiveNumber&);
EffectiveScalarFunction operator+(const EffectiveNumber&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator-(const EffectiveNumber&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator*(const EffectiveNumber&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator/(const EffectiveNumber&, const EffectiveScalarFunction&);

EffectiveScalarFunction pow(const EffectiveScalarFunction&, Int);
EffectiveScalarFunction neg(const EffectiveScalarFunction&);
EffectiveScalarFunction rec(const EffectiveScalarFunction&);
EffectiveScalarFunction sqr(const EffectiveScalarFunction&);
EffectiveScalarFunction sqrt(const EffectiveScalarFunction&);
EffectiveScalarFunction exp(const EffectiveScalarFunction&);
EffectiveScalarFunction log(const EffectiveScalarFunction&);
EffectiveScalarFunction sin(const EffectiveScalarFunction&);
EffectiveScalarFunction cos(const EffectiveScalarFunction&);
EffectiveScalarFunction tan(const EffectiveScalarFunction&);
EffectiveScalarFunction atan(const EffectiveScalarFunction&);

ValidatedScalarFunction operator+(const ValidatedScalarFunction&, const ValidatedScalarFunction&);
ValidatedScalarFunction operator-(const ValidatedScalarFunction&, const ValidatedScalarFunction&);
ValidatedScalarFunction operator*(const ValidatedScalarFunction&, const ValidatedScalarFunction&);
ValidatedScalarFunction operator/(const ValidatedScalarFunction&, const ValidatedScalarFunction&);


/*
//! \ingroup FunctionModule
//! \brief A generic vector function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X^m\f$.
template<class P>
class VectorFunction
{
    static_assert(IsStronger<P,ApproximateTag>::value,"P must be an information level/paradigm.");
    typedef CanonicalNumberType<P> X;
  private:
    std::shared_ptr< const VectorFunctionInterface<P> > _ptr;
  public:
    typedef X NumericType;

    static VectorFunction<P> zeros(SizeType rs, SizeType as);
    static VectorFunction<P> identity(SizeType n);

    VectorFunction();
    VectorFunction(SizeType rs, SizeType as);
    VectorFunction(SizeType as, const List< Formula<X> >& flst);
    VectorFunction(SizeType rs, ScalarFunction<P> sf);

    explicit VectorFunction(VectorFunctionInterface<P>* fptr) : _ptr(fptr) { }
    VectorFunction(std::shared_ptr< VectorFunctionInterface<P> > fptr) : _ptr(fptr) { }
    VectorFunction(const VectorFunctionInterface<P>& fref) : _ptr(fref._clone()) { }
    std::shared_ptr< const VectorFunctionInterface<P> > managed_pointer() const { return this->_ptr; }
    const VectorFunctionInterface<P>* raw_pointer() const { return this->_ptr.operator->(); }
    const VectorFunctionInterface<P>& reference() const { return this->_ptr.operator*(); }
    operator const VectorFunctionInterface<P>& () const { return *this->_ptr; }

    VectorFunction(const List< ScalarFunction<P> >& lsf);
    VectorFunction(const Vector< ScalarFunction<P> >& lsf);
    VectorFunction(InitializerList< ScalarFunction<P> > lsf);
    template<class PP> VectorFunction(const List< ScalarFunction<PP> >& lsf, EnableIf< IsStronger<PP,P>, Void >* = 0) {
        *this=VectorFunction<P>(List< ScalarFunction<P> >(lsf)); }
    template<class PP> VectorFunction(const VectorFunction<PP>& vf, EnableIf< IsStronger<PP,P>, Void >* = 0)
        : _ptr(std::dynamic_pointer_cast< const VectorFunctionInterface<P> >(vf.managed_pointer())) { }

    ScalarFunction<P> get(SizeType i) const { return ScalarFunction<P>(this->_ptr->_get(i)); }
    //Void set(SizeType i, ScalarFunction<P> f) { this->_ptr->_set(i,f); };
    Void set(SizeType i, ScalarFunction<P> f);

    ScalarFunction<P> operator[](SizeType i) const { return this->get(i); }
    VectorFunctionElementReference<P> operator[](SizeType i);

    SizeType result_size() const { return this->_ptr->result_size(); }
    SizeType argument_size() const { return this->_ptr->argument_size(); }

    template<class XX> auto evaluate(const Vector<XX>& x) const -> decltype(this->_ptr->evaluate(x)) { return this->_ptr->evaluate(x); }
    template<class XX> auto operator() (const Vector<XX>& x) const -> decltype(this->_ptr->evaluate(x)) { return this->_ptr->evaluate(x); }

    template<class XX> auto jacobian(const Vector<XX>& x) const -> decltype(this->_ptr->jacobian(x)) {
        return this->_ptr->jacobian(x); }
    template<class XX> auto differentials(const Vector<XX>& x, DegreeType d) const -> decltype(this->_ptr->differentials(x,d)) {
        return this->_ptr->differentials(x,d); }

    OutputStream& write(OutputStream& os) const { return this->_ptr->write(os); }
};

template<class P> inline OutputStream& operator<<(OutputStream& os, const VectorFunction<P>& f) { return f.write(os); }

template<class P, class XX> inline Vector<XX> evaluate(const VectorFunction<P>& f, const Vector<XX>& x) { return f(x); }
template<class P, class XX> inline Matrix<XX> jacobian(const VectorFunction<P>& f, const Vector<XX>& x) { return f.jacobian(x); }

inline Matrix<ValidatedNumber> VectorFunctionInterface<ValidatedTag>::jacobian(const Vector<ExactNumber>& x) const {
    return this->jacobian(Vector<ValidatedNumber>(x)); }
*/

/*
template<class X> VectorFunction<P> operator*(const ScalarFunction<X>& sf, const Vector<X>& e);
template<class X> VectorFunction<X> operator+(const VectorFunction<X>& f1, const VectorFunction<X>& f2);
template<class X> VectorFunction<X> operator-(const VectorFunction<X>& f1, const VectorFunction<X>& f2);
template<class X> VectorFunction<X> operator*(const VectorFunction<X>& vf, const ScalarFunction<X>& sf);
template<class X> VectorFunction<X> operator*(const ScalarFunction<X>& sf, const VectorFunction<X>& vf);

template<class X> ScalarFunction<X> embed(SizeType as1, const ScalarFunction<X>& f2, SizeType as3);
template<class X> VectorFunction<X> embed(SizeType as1, const VectorFunction<X>& f2, SizeType as3);

template<class X> VectorFunction<X> join(const ScalarFunction<X>& f1, const ScalarFunction<X>& f2);
template<class X> VectorFunction<X> join(const ScalarFunction<X>& f1, const VectorFunction<X>& f2);
template<class X> VectorFunction<X> join(const VectorFunction<X>& f1, const ScalarFunction<X>& f2);
template<class X> VectorFunction<X> join(const VectorFunction<X>& f1, const VectorFunction<X>& f2);

template<class X> ScalarFunction<X> compose(const ScalarFunction<X>& f, const VectorFunction<X>& g);
template<class X> VectorFunction<X> compose(const VectorFunction<X>& f, const VectorFunction<X>& g);

template<class X> ScalarFunction<X> lie_derivative(const ScalarFunction<X>& g, const VectorFunction<X>& f);
*/

EffectiveVectorFunction operator*(const EffectiveScalarFunction& sf, const Vector<EffectiveNumber>& e);
EffectiveVectorFunction operator+(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2);
EffectiveVectorFunction operator-(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2);
EffectiveVectorFunction operator*(const EffectiveVectorFunction& vf, const EffectiveScalarFunction& sf);
EffectiveVectorFunction operator*(const EffectiveScalarFunction& sf, const EffectiveVectorFunction& vf);
EffectiveVectorFunction operator*(const EffectiveNumber& c, const EffectiveVectorFunction& vf);

EffectiveScalarFunction embed(SizeType as1, const EffectiveScalarFunction& f2, SizeType as3);
EffectiveVectorFunction embed(SizeType as1, const EffectiveVectorFunction& f2, SizeType as3);

EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2);
EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveVectorFunction& f2);
EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveScalarFunction& f2);
EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2);

EffectiveScalarFunction compose(const EffectiveScalarFunction& f, const EffectiveVectorFunction& g);
EffectiveVectorFunction compose(const EffectiveVectorFunction& f, const EffectiveVectorFunction& g);

EffectiveScalarFunction lie_derivative(const EffectiveScalarFunction& g, const EffectiveVectorFunction& f);

Formula<EffectiveNumber> formula(const EffectiveScalarFunction& f);
Vector< Formula<EffectiveNumber> > formula(const EffectiveVectorFunction& f);


ValidatedScalarFunction operator-(const ValidatedScalarFunction&, const ValidatedScalarFunction&);
ValidatedScalarFunction operator-(const ValidatedScalarFunction&, const ValidatedNumber&);
ValidatedScalarFunction operator-(const ValidatedNumber&, const ValidatedScalarFunction&);
ValidatedVectorFunction operator-(const ValidatedVectorFunction&, const ValidatedVectorFunction&);
ValidatedVectorFunction join(const ValidatedVectorFunction& f1, const ValidatedVectorFunction& f2);
ValidatedScalarFunction compose(const ValidatedScalarFunction& f, const ValidatedVectorFunction& g);
ValidatedVectorFunction compose(const ValidatedVectorFunction& f, const ValidatedVectorFunction& g);



template<class P, class D>
struct VectorFunctionElementReference {
    typedef IntervalDomain SC; typedef BoxDomain VC;
    Function<P,D,VC>& _vf; SizeType _i;
    VectorFunctionElementReference<P,D>(Function<P,D,VC>& vf, SizeType i) : _vf(vf), _i(i) { }
    template<class WP> operator Function<WP,D,SC> () const;
    Void operator=(const Function<P,D,SC>& sf);
    VectorFunctionElementReference<P,D>& operator=(const VectorFunctionElementReference<P,D>& sfr);
    template<class XX> XX evaluate(const Vector<XX> & x) const;
    template<class XX> XX operator()(const Vector<XX> & x) const;
};

template<class P, class D> template<class WP> inline
VectorFunctionElementReference<P,D>::operator Function<WP,D,SC> () const {
    return this->_vf.get(this->_i);
}

template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(ScalarFunction<P,D>& sf, SizeType i) {
    ARIADNE_ASSERT(false); }
template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(VectorFunction<P,D>& vf, SizeType i) {
    return VectorFunctionElementReference<P,D>(vf,i); }

template<class P, class D, class C> inline VectorFunctionElementReference<P,D> Function<P,D,C>::operator[](SizeType i) {
    return make_element_reference(*this,i); }

template<class P, class D, class C> inline ScalarFunction<P,D> Function<P,D,C>::operator[](SizeType i) const {
    return this->get(i); }

template<class P, class D> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionElementReference<P,D>& vfe) {
    return  os << static_cast< Function<P,D,IntervalDomain> >(vfe); }

template<class P, class D> inline Void VectorFunctionElementReference<P,D>::operator=(const Function<P,D,SC>& sf) {
    _vf.set(_i,sf); }
template<class P, class D> inline VectorFunctionElementReference<P,D>& VectorFunctionElementReference<P,D>::operator=(const VectorFunctionElementReference<P,D>& sfr) {
    _vf.set(_i,static_cast<Function<P,D,SC>>(sfr)); return *this; }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::evaluate(const Vector<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::operator()(const Vector<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }



inline UpperInterval apply(ScalarFunction<ValidatedTag>const& f, const Vector<UpperInterval>& x) {
    return static_cast<UpperInterval>(f(reinterpret_cast<Vector<ValidatedNumber>const&>(x))); }
inline Box<UpperInterval> apply(VectorFunction<ValidatedTag>const& f, const Vector<UpperInterval>& x) {
    return static_cast<Vector<UpperInterval>>(f(reinterpret_cast<Vector<ValidatedNumber>const&>(x))); }
inline Vector<Differential<UpperInterval>> apply_derivative(VectorFunction<ValidatedTag>const& f, const Vector<Differential<UpperInterval>>& x) {
    return static_cast<Vector<Differential<UpperInterval>>>(f(reinterpret_cast<Vector<Differential<ValidatedNumber>>const&>(x))); }


inline Covector<UpperInterval> gradient_range(ValidatedScalarFunction const& f, const Vector<UpperInterval>& x) {
    return static_cast<Covector<UpperInterval>>(static_cast<Covector<ValidatedFloat>>(gradient(f,reinterpret_cast<Vector<ValidatedFloat>const&>(x)))); }
inline Matrix<UpperInterval> jacobian_range(ValidatedVectorFunction const& f, const Vector<UpperInterval>& x) {
    return static_cast<Matrix<UpperInterval>>(static_cast<Matrix<ValidatedFloat>>(jacobian(f,reinterpret_cast<Vector<ValidatedFloat>const&>(x)))); }

/*
inline Matrix<UpperInterval> jacobian(VectorFunction<ValidatedTag>const& f, const Vector<UpperInterval>& x) {
    return static_cast<Matrix<UpperInterval>>(f.jacobian(reinterpret_cast<Vector<ValidatedFloat>const&>(x))); }
inline Matrix<UpperInterval> jacobian(VectorFunction<ValidatedTag>const& f, const Vector<ExactInterval>& x) {
    return static_cast<Matrix<UpperInterval>>(f.jacobian(reinterpret_cast<Vector<ValidatedFloat>const&>(x))); }
inline Matrix<UpperInterval> jacobian_range(VectorFunction<ValidatedTag>const& f, const Vector<UpperInterval>& x) {
    return static_cast<Matrix<UpperInterval>>(f.jacobian(reinterpret_cast<Vector<ValidatedFloat>const&>(x))); }

// FIXME: Needed to override templated gradient and jacobian
inline Covector<UpperInterval> gradient(ScalarFunction<EffectiveTag>const& f, const Vector<UpperInterval>& x) {
    return static_cast<Covector<UpperInterval>>(gradient(f,reinterpret_cast<Vector<ValidatedFloat>const&>(x))); }
inline Covector<UpperInterval> gradient(ScalarFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return gradient(f,static_cast<Vector<UpperInterval>>(x)); }
inline Matrix<UpperInterval> jacobian(VectorFunction<EffectiveTag>const& f, const Vector<UpperInterval>& x) {
    return static_cast<Matrix<UpperInterval>>(f.jacobian(reinterpret_cast<Vector<ValidatedFloat>const&>(x))); }
inline Matrix<UpperInterval> jacobian(VectorFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return jacobian(f,static_cast<Vector<UpperInterval>>(x)); }
*/

template<class P> class FunctionFactory;
typedef FunctionFactory<ValidatedTag> ValidatedFunctionFactory;

template<>
class FunctionFactory<ValidatedTag>
{
    std::shared_ptr< const FunctionFactoryInterface<ValidatedTag> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>* ptr) : _ptr(ptr) { }
    FunctionFactory(std::shared_ptr< const FunctionFactoryInterface<ValidatedTag> > ptr) : _ptr(ptr) { }
    inline ValidatedScalarFunction create(const ExactBox& d, const ValidatedScalarFunctionInterface& f) const;
    inline ValidatedVectorFunction create(const ExactBox& d, const ValidatedVectorFunctionInterface& f) const;
    inline ValidatedScalarFunction create_zero(const ExactBox& d) const;
    inline ValidatedVectorFunction create_identity(const ExactBox& d) const;
    friend OutputStream& operator<<(OutputStream& os, const FunctionFactory<ValidatedTag>& factory);
};

inline ValidatedScalarFunction FunctionFactoryInterface<ValidatedTag>::create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const {
    return ValidatedScalarFunction(SharedPointer<ValidatedScalarFunctionInterface>(this->_create(domain,function))); }
inline ValidatedVectorFunction FunctionFactoryInterface<ValidatedTag>::create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const {
    return ValidatedVectorFunction(SharedPointer<ValidatedVectorFunctionInterface>(this->_create(domain,function))); }
inline ValidatedScalarFunction FunctionFactoryInterface<ValidatedTag>::create_zero(const ExactBox& domain) const {
    return this->create(domain,EffectiveScalarFunction::zero(domain.dimension())); }
inline ValidatedVectorFunction FunctionFactoryInterface<ValidatedTag>::create_identity(const ExactBox& domain) const {
    return this->create(domain,EffectiveVectorFunction::identity(domain.dimension())); }

inline ValidatedScalarFunction FunctionFactory<ValidatedTag>::create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedVectorFunction FunctionFactory<ValidatedTag>::create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedScalarFunction FunctionFactory<ValidatedTag>::create_zero(const ExactBox& domain) const {
    return this->_ptr->create_zero(domain); }
inline ValidatedVectorFunction FunctionFactory<ValidatedTag>::create_identity(const ExactBox& domain) const {
    return this->_ptr->create_identity(domain); }

inline OutputStream& operator<<(OutputStream& os, const ValidatedFunctionFactory& factory) {
    factory._ptr->write(os); return os; }

} // namespace Ariadne

#endif
