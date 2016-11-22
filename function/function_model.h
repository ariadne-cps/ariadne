/***************************************************************************
 *            function_model.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file function_model.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_MODEL_H
#define ARIADNE_FUNCTION_MODEL_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function_model_interface.h"

#include "numeric/operators.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/operations.h"
#include "geometry/box.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function.h"

namespace Ariadne {

template<class X> class ScalarFunctionModel;
template<class X> class VectorFunctionModel;

typedef ScalarFunctionModelInterface<ValidatedTag> ValidatedScalarFunctionModelInterface;
typedef VectorFunctionModelInterface<ValidatedTag> ValidatedVectorFunctionModelInterface;

typedef ScalarFunctionModel<ValidatedTag> ValidatedScalarFunctionModel;
typedef VectorFunctionModel<ValidatedTag> ValidatedVectorFunctionModel;

template<class P> class FuntionModelFactoryInterface;
typedef FuntionModelFactoryInterface<ValidatedTag> ValidatedFuntionModelFactoryInterface;
template<class P> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedTag> ValidatedFunctionModelFactory;


template<class P, class F> class TaylorModel;
template<class M> class FunctionPatch;
typedef FunctionPatch<TaylorModel<ValidatedTag,Float64>> ScalarTaylorFunction;


template<class P> class FunctionModelBuilder {
    ScalarFunctionModel<P> const& _prototype;
  public:
    FunctionModelBuilder(ScalarFunctionModel<P> const& f);

    ScalarFunctionModel<P> create(ScalarFunction<P> const& f);
    VectorFunctionModel<P> create(VectorFunction<P> const& f);
    CanonicalNumericType<P> create(Number<P> const& c);

    ScalarFunctionModel<P> const& create(ScalarFunctionModel<P> const& f) { return f; }
    CanonicalNumericType<P> const& create(CanonicalNumericType<P> const& c) { return c; }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on singleton domains.
template<class P> class ScalarFunctionModel
    : public DispatchAlgebraOperators<ScalarFunctionModel<P>, CanonicalNumericType<P>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P>>
{
  public:
    typedef ScalarFunction<P> GenericType;
    typedef ExactBoxType DomainType;
    typedef CanonicalCoefficientType<P> CoefficientType;
    typedef CanonicalErrorType<P> ErrorType;
    typedef CanonicalNumericType<P> NumericType;
  public:
    clone_on_copy_ptr< ScalarFunctionModelInterface<P> > _ptr;
  public:
    ScalarFunctionModel() : _ptr() { }
    ScalarFunctionModel(ScalarFunctionModelInterface<P>* p) : _ptr(p) { }
    ScalarFunctionModel(const shared_ptr<const ScalarFunctionModelInterface<P>> p) : _ptr(p->_clone()) { }
    ScalarFunctionModel(const ScalarFunctionModel<P>& f) : _ptr(f._ptr) { }
    ScalarFunctionModel(const ScalarFunctionModelInterface<P>& f) : _ptr(f._clone()) { }
    ScalarFunctionModel(const ScalarFunction<P>& f) : _ptr(dynamic_cast<ScalarFunctionModelInterface<P>*>(f.raw_pointer()->_clone())) { }
    operator ScalarFunction<P>() const { return ScalarFunction<P>(this->_ptr->_clone()); }
    operator ScalarFunctionModelInterface<P>& () { return *_ptr; }
    operator const ScalarFunctionModelInterface<P>& () const { return *_ptr; }
    const ScalarFunctionModelInterface<P>* raw_pointer() const { return _ptr.operator->(); }
    ScalarFunctionModelInterface<P>& reference() { return *_ptr; }
    const ScalarFunctionModelInterface<P>& reference() const { return *_ptr; }
    ScalarFunctionModel<P> create_zero() const;
    ScalarFunctionModel<P> create_constant(const CanonicalNumericType<P>& c) const;
    ScalarFunctionModel<P> create_coordinate(SizeType j) const;
    VectorFunctionModel<P> create_identity() const;
    ScalarFunctionModel<P> create(const ScalarFunction<P>& f) const;
    CanonicalNumericType<P> create(const Number<P>& f) const;
    Vector<ScalarFunctionModel<P>> create_coordinates(DomainType const&) const;
    ScalarFunctionModel<P>& operator=(const Number<P>& c);
    ScalarFunctionModel<P>& operator=(const CanonicalNumericType<P>& c);
    ScalarFunctionModel<P>& operator=(const ScalarFunction<P>& f);
    ScalarFunctionModel<P>& operator=(const ScalarTaylorFunction& f);
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class XX> XX operator() (const Vector<XX>& x) const {
        return this->_ptr->_evaluate(x); }
    template<class XX> XX evaluate(const Vector<XX>& x) const {
        return this->_ptr->_evaluate(x); }
    inline ExactBoxType const domain() const { return this->_ptr->domain(); }
    inline ExactIntervalType const codomain() const { return this->_ptr->codomain(); }
    inline UpperIntervalType const range() const { return this->_ptr->range(); }

    inline CoefficientType value() const { return this->_ptr->value(); }
    inline CoefficientType gradient_value(SizeType j) const { return this->_ptr->gradient_value(j); }
    inline ErrorType error() const { return this->_ptr->error(); }

    inline Void set_error(const ErrorType& e) { return this->_ptr->set_error(e); }
    inline Void set_error(Nat e) { return this->_ptr->set_error(ErrorType(e,this->error().precision())); }
    inline Void clobber() { return this->_ptr->clobber(); }

    inline ScalarFunctionModel<P> apply(Operator op) const { return this->_ptr->_apply(op); }
    inline Void restrict(const ExactBoxType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelBuilder<P> factory(ScalarFunctionModel<P> const& f) {
        return FunctionModelBuilder<P>(f); }
  public:
    friend inline ScalarFunctionModel<P> operator/(const ScalarFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2) {
        return mul(f1,rec(f2)); }
    friend inline ScalarFunctionModel<P> neg(ScalarFunctionModel<P> f) {
        f._ptr->_imul(CanonicalNumericType<P>(-1)); return std::move(f); }
    friend inline ScalarFunctionModel<P> add(ScalarFunctionModel<P> f1, const ScalarFunctionModel<P>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P>(+1),f2); return std::move(f1); }
    friend inline ScalarFunctionModel<P> sub(ScalarFunctionModel<P> f1, const ScalarFunctionModel<P>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P>(-1),f2); return std::move(f1); }
    friend inline ScalarFunctionModel<P> mul(const ScalarFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2) {
        ScalarFunctionModel<P> r=f1.create_zero(); r._ptr->_ifma(f1,f2); return r; }
    friend inline ScalarFunctionModel<P> div(const ScalarFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2) {
        return mul(f1,rec(f2)); }
    friend inline ScalarFunctionModel<P> add(ScalarFunctionModel<P> f1, const CanonicalNumericType<P>& c2) {
        f1._ptr->_iadd(c2); return std::move(f1); }
    friend inline ScalarFunctionModel<P> mul(ScalarFunctionModel<P> f1, const CanonicalNumericType<P>& c2) {
        f1._ptr->_imul(c2); return std::move(f1); }
    friend inline ScalarFunctionModel<P> add(const CanonicalNumericType<P>& c1, ScalarFunctionModel<P> f2) {
        f2._ptr->_iadd(c1); return std::move(f2); }
    friend inline ScalarFunctionModel<P> mul(const CanonicalNumericType<P>& c1, ScalarFunctionModel<P> f2) {
        f2._ptr->_imul(c1); return std::move(f2); }

    friend inline ScalarFunctionModel<P> sub(ScalarFunctionModel<P> f1, const CanonicalNumericType<P>& c2) {
        return add(std::move(f1),neg(c2)); }
    friend inline ScalarFunctionModel<P> div(ScalarFunctionModel<P> f1, const CanonicalNumericType<P>& c2) {
        return mul(std::move(f1),rec(c2)); }
    friend inline ScalarFunctionModel<P> sub(const CanonicalNumericType<P>& c1, ScalarFunctionModel<P> f2) {
        return add(neg(std::move(f2)),c1); }
    friend inline ScalarFunctionModel<P> div(const CanonicalNumericType<P>& c1, ScalarFunctionModel<P> f2) {
        return mul(rec(std::move(f2)),c1); }

    friend inline ScalarFunctionModel<P> pow(const ScalarFunctionModel<P>& f1, Int n2) {
        return generic_pow(f1,n2); }

    friend inline ScalarFunctionModel<P> rec(const ScalarFunctionModel<P>& f) {
        return f.apply(Rec()); }

  public:
    friend NormType norm(const ScalarFunctionModel<P>& f) { return f._ptr->_norm(); }
    friend ScalarFunctionModel<P> derivative(const ScalarFunctionModel<P>& f, SizeType j) { return f._ptr->_derivative(j); }
    friend ScalarFunctionModel<P> antiderivative(const ScalarFunctionModel<P>& f, SizeType j) { return f._ptr->_antiderivative(j); }
    friend ScalarFunctionModel<P> antiderivative(const ScalarFunctionModel<P>& f, SizeType j, CanonicalNumericType<P> c) { return f._ptr->_antiderivative(j,c); }

    friend ScalarFunctionModel<P> embed(const ExactBoxType& d1, const ScalarFunctionModel<P>& f, const ExactBoxType& d2) {
        return f._ptr->_embed(d1,d2); }
    friend ScalarFunctionModel<P> embed(const ExactBoxType& d, const ScalarFunctionModel<P>& f) {
        return embed(d,f,ExactBoxType()); }
    friend ScalarFunctionModel<P> embed(const ScalarFunctionModel<P>& f, const ExactBoxType& d) {
        return embed(ExactBoxType(),f,d); }
    friend ScalarFunctionModel<P> embed(const ScalarFunctionModel<P>& f, const ExactIntervalType& d) {
        return embed(f,ExactBoxType(1,d)); }
    friend ScalarFunctionModel<P> restrict(const ScalarFunctionModel<P>& f, const ExactBoxType& d) {
        return f._ptr->_restriction(d); }
    friend ScalarFunctionModel<P> restriction(const ScalarFunctionModel<P>& f, const ExactBoxType& d) {
        return f._ptr->_restriction(d); }

    friend VectorFunctionModel<P> join(const ScalarFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2);
    friend VectorFunctionModel<P> combine(const ScalarFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2);
  public:
    friend OutputStream& operator<<(OutputStream& os, const ScalarFunctionModel<P>& f) {
        return os <<  f.operator ScalarFunction<P>(); }
};

inline ScalarFunctionModel<ValidatedTag> refinement(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_refinement(f2); }
inline Boolean inconsistent(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_inconsistent(f2); }
inline Boolean refines(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_refines(f2); }

// FIXME: Should not be needed since ScalarFunctionModel has a representation
template<class P> inline ScalarFunctionModel<P> embed(const ScalarFunction<P>& f, const ExactIntervalType& d) {
    return embed(ScalarFunctionModel<P>(f),d); }

template<class P> inline ScalarFunctionModel<P>& ScalarFunctionModel<P>::operator=(const CanonicalNumericType<P>& c) {
        (*this)*=CanonicalNumericType<P>(0); (*this)+=c; return *this; }
template<class P> inline ScalarFunctionModel<P>& ScalarFunctionModel<P>::operator=(const ScalarFunction<P>& f) {
        return (*this)=factory(*this).create(f); }



template<class V> struct Element;

template<class M> class FunctionPatch;
template<class M> class VectorFunctionPatch;
template<class M> struct Element<VectorFunctionPatch<M>> { typedef FunctionPatch<M> Type; };

typedef FunctionPatch<ValidatedTaylorModel> ScalarTaylorFunction;

template<class P> class VectorFunctionModelElement
    : public DispatchAlgebraOperators<ScalarFunctionModel<P>, CanonicalNumericType<P>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P>>
{
    VectorFunctionModel<P>* _p; SizeType _i;
  public:
    typedef typename ScalarFunctionModel<P>::GenericType GenericType;
    operator const ScalarFunctionModel<P> () const;
    VectorFunctionModelElement(VectorFunctionModel<P>* p, SizeType i) : _p(p), _i(i) { }
    VectorFunctionModelElement<P>& operator=(const ScalarFunctionModel<P>& sf) { _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<P>& operator=(const VectorFunctionModelElement<P>& sf) { return this->operator=(static_cast<ScalarFunctionModel<P>const>(sf)); }
    Void clobber() { ScalarFunctionModel<P> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    Void set_error(CanonicalErrorType<P> e) const { ScalarFunctionModel<P> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    Void set_error(Nat e) const { ScalarFunctionModel<P> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    friend ScalarFunctionModel<P> antiderivative(VectorFunctionModelElement<P> const& f, SizeType k) {
        return antiderivative(ScalarFunctionModel<P>(f),k); }
};
template<class P> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<P>& function) {
    return os << static_cast< const ScalarFunctionModel<P> >(function);
}

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on singleton domains.
template<class P> class VectorFunctionModel
{
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<P> > _ptr;
  public:
    typedef CanonicalCoefficientType<P> CoefficientType;
    typedef CanonicalErrorType<P> ErrorType;
    typedef CanonicalNumericType<P> NumericType;
  public:
    inline VectorFunctionModel() : _ptr() { }
    inline VectorFunctionModel(SizeType n, const ScalarFunctionModelInterface<P>& sf)
        : _ptr(sf._create_vector(n)) { for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline VectorFunctionModel(Array<ScalarFunctionModel<P>> const& asf)
        : VectorFunctionModel(asf.size(),asf[0]) { for(SizeType i=0; i!=asf.size(); ++i) { (*this)[i]=asf[i]; } }
    inline VectorFunctionModel(List<ScalarFunctionModel<P>> const& lsf)
        : VectorFunctionModel(lsf.size(),lsf[0]) { for(SizeType i=0; i!=lsf.size(); ++i) { (*this)[i]=lsf[i]; } }
    inline VectorFunctionModel(VectorFunctionModelInterface<P>* p) : _ptr(p) { }
    inline VectorFunctionModel(const VectorFunctionModelInterface<P>& f) : _ptr(f._clone()) { }
    inline VectorFunctionModel(const VectorFunctionModel<P>& f) : _ptr(f._ptr) { }
    inline operator const VectorFunctionModelInterface<P>& () const { return *_ptr; }
    inline operator ValidatedVectorFunction () const { return ValidatedVectorFunction(*_ptr); }
    inline const VectorFunctionModelInterface<P>* raw_pointer() const { return _ptr.operator->(); }
    inline const VectorFunctionModelInterface<P>& reference() const { return *_ptr; }
    inline VectorFunctionModelInterface<P>& reference() { return *_ptr; }
    inline ScalarFunctionModel<P> create_zero() const { return this->_ptr->_create_zero(); }
    inline VectorFunctionModel<P> create_identity() const { return this->_ptr->_create_identity(); }
    inline VectorFunctionModel<P> create(VectorFunction<P> const& vf) const { return compose(vf,this->create_identity()); }
    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    inline ScalarFunctionModel<P> const get(SizeType i) const { return this->_ptr->_get(i); }
    inline Void set(SizeType i, ScalarFunctionModel<P> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<P> const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionModelElement<P> operator[](SizeType i) { return VectorFunctionModelElement<P>(this,i); }
    inline ExactBoxType const domain() const { return this->_ptr->domain(); }
    inline ExactBoxType const codomain() const { return this->_ptr->codomain(); }
    inline UpperBoxType const range() const { return this->_ptr->range(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->errors(); }
    inline ErrorType const error() const { return this->_ptr->error(); }
    inline Void clobber() { this->_ptr->clobber(); }
    inline Matrix<NumericType> const jacobian(const Vector<NumericType>& x) const {
        Vector<Differential<NumericType>> dfx=this->_ptr->_evaluate(Differential<NumericType>::variables(1u,x));
        return dfx.jacobian(); }

    inline Void restrict(const ExactBoxType& d) { this->_ptr->restrict(d); }
  public:
    friend inline ScalarFunctionModel<P> compose(const ScalarFunction<P>& f, const VectorFunctionModel<P>& g) {
        return g._ptr->_compose(f); }
    friend inline ScalarFunctionModel<P> compose(const ScalarFunctionModel<P>& f, const VectorFunctionModel<P>& g) {
        return g._ptr->_compose(f); }
    friend inline VectorFunctionModel<P> compose(const VectorFunction<P>& f, const VectorFunctionModel<P>& g) {
        return g._ptr->_compose(f); }
    friend inline VectorFunctionModel<P> compose(const VectorFunctionModel<P>& f, const VectorFunctionModel<P>& g) {
        return g._ptr->_compose(f); }

    friend inline ScalarFunctionModel<P> unchecked_compose(const ScalarFunctionModel<P>& f, const VectorFunctionModel<P>& g) {
        return g._ptr->_unchecked_compose(f); }
    friend inline VectorFunctionModel<P> unchecked_compose(const VectorFunctionModel<P>& f, const VectorFunctionModel<P>& g) {
        return g._ptr->_unchecked_compose(f); }

    friend inline VectorFunctionModel<P> operator+(const VectorFunctionModel<P>& f) {
        return f._ptr->_clone(); }
    friend inline VectorFunctionModel<P> operator-(const VectorFunctionModel<P>& f) {
        VectorFunctionModel<P> r=f; for(SizeType i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
    friend inline VectorFunctionModel<P> operator+(const VectorFunctionModel<P>& f1, const VectorFunctionModel<P>& f2) {
        VectorFunctionModel<P> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+f2[i]; } return r; }
    friend inline VectorFunctionModel<P> operator-(const VectorFunctionModel<P>& f1, const VectorFunctionModel<P>& f2) {
        VectorFunctionModel<P> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
    friend inline VectorFunctionModel<P> operator+(const VectorFunctionModel<P>& f1, const Vector<CanonicalNumericType<P>>& c2) {
        VectorFunctionModel<P> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
    friend inline VectorFunctionModel<P> operator-(const VectorFunctionModel<P>& f1, const Vector<CanonicalNumericType<P>>& c2) {
        VectorFunctionModel<P> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
    friend inline VectorFunctionModel<P> operator+(const Vector<CanonicalNumericType<P>>& c1, const VectorFunctionModel<P>& f2);
    friend inline VectorFunctionModel<P> operator-(const Vector<CanonicalNumericType<P>>& c1, const VectorFunctionModel<P>& f2);
    friend inline VectorFunctionModel<P> operator*(const VectorFunctionModel<P>& f1, const CanonicalNumericType<P>& c2) {
        VectorFunctionModel<P> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
    friend inline VectorFunctionModel<P> operator*(const CanonicalNumericType<P>& c1, const VectorFunctionModel<P>& f2) {
        VectorFunctionModel<P> r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

    friend inline VectorFunctionModel<P> operator+(const VectorFunctionModel<P>& f1, const VectorFunction<P>& f2) {
        return f1+f1.create(f2); }
    friend inline VectorFunctionModel<P> operator-(const VectorFunctionModel<P>& f1, const VectorFunction<P>& f2) {
        return f1-f1.create(f2); }
    friend inline VectorFunctionModel<P> operator+(const VectorFunction<P>& f1, const VectorFunctionModel<P>& f2) {
        return f2.create(f1)+f2; }
    friend inline VectorFunctionModel<P> operator-(const VectorFunction<P>& f1, const VectorFunctionModel<P>& f2) {
        return f2.create(f1)-f2; }


  public:
    friend NormType norm(const VectorFunctionModel<P>& f) {
        return f._ptr->_norm(); }
    friend VectorFunctionModel<P> embed(const ExactBoxType& d1, const VectorFunctionModel<P>& f, const ExactBoxType& d2) {
        return f._ptr->_embed(d1,d2); }
    friend VectorFunctionModel<P> embed(const VectorFunctionModel<P>& f, const ExactBoxType& d) {
        return embed(ExactBoxType(),f,d); }
    friend VectorFunctionModel<P> embed(const VectorFunctionModel<P>& f, const ExactIntervalType& d) {
        return embed(f,ExactBoxType(1,d)); }
    friend VectorFunctionModel<P> restrict(const VectorFunctionModel<P>& f, const ExactBoxType& d) {
        VectorFunctionModelInterface<P>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }
    friend VectorFunctionModel<P> restriction(const VectorFunctionModel<P>& f, const ExactBoxType& d) {
        VectorFunctionModelInterface<P>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }

    friend Vector<CanonicalNumericType<P>> unchecked_evaluate(const VectorFunctionModel<P>& f, const Vector<CanonicalNumericType<P>>& x) {
        return f._ptr->_unchecked_evaluate(x); }

    friend ScalarFunctionModel<P> unchecked_compose(const ScalarFunction<P>& f, const VectorFunctionModel<P>& g) {
        ScalarFunctionModelInterface<P> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(ScalarFunctionModel<P>(*fptr),g); } else { return compose(f,g); } }
    friend VectorFunctionModel<P> unchecked_compose(const VectorFunction<P>& f, const VectorFunctionModel<P>& g) {
        VectorFunctionModelInterface<P> const* fptr = dynamic_cast<VectorFunctionModelInterface<P> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(VectorFunctionModel<P>(*fptr),g); } else { return compose(f,g); } }

    friend VectorFunctionModel<P> join(const ScalarFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2);
    friend VectorFunctionModel<P> join(const ScalarFunctionModel<P>& f1, const VectorFunctionModel<P>& f2);
    friend VectorFunctionModel<P> join(const VectorFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2) {
        VectorFunctionModel<P> r=f1._ptr->_clone(); r._ptr->_adjoin(f2); return r; }
    friend VectorFunctionModel<P> join(const VectorFunctionModel<P>& f1, const VectorFunctionModel<P>& f2) {
        return f1._ptr->_join(f2); }

    friend VectorFunctionModel<P> combine(const ScalarFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2) {
        return VectorFunctionModel<P>(1,f1)._ptr->_combine(VectorFunctionModel<P>(1,f2)); };;
    friend VectorFunctionModel<P> combine(const ScalarFunctionModel<P>& f1, const VectorFunctionModel<P>& f2) {
        return VectorFunctionModel<P>(1,f1)._ptr->_combine(f2); };;
    friend VectorFunctionModel<P> combine(const VectorFunctionModel<P>& f1, const ScalarFunctionModel<P>& f2) {
        return f1._ptr->_combine(VectorFunctionModel<P>(1,f2)); };
    friend VectorFunctionModel<P> combine(const VectorFunctionModel<P>& f1, const VectorFunctionModel<P>& f2) {
        return f1._ptr->_combine(f2); }

    friend inline VectorFunctionModel<ValidatedTag> refinement(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
        ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"refinement(f1,f2): f1="<<f1<<", f2="<<f2<<")");
        VectorFunctionModel<ValidatedTag> r=+f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=refinement(f1[i],f2[i]); } return r; }

    friend VectorFunctionModel<P> antiderivative(const VectorFunctionModel<P>& f, SizeType j) {
        VectorFunctionModel<P> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
        return r;
    }

    friend VectorFunctionModel<P> antiderivative(const VectorFunctionModel<P>& f, SizeType j, const CanonicalNumericType<P>& c) {
        VectorFunctionModel<P> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j,c); }
        return r;
    }

    friend VectorFunctionModel<P> partial_evaluate(const VectorFunctionModel<P>& f, SizeType j, const CanonicalNumericType<P>& c) {
        return f._ptr->_partial_evaluate(j,c); }

    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionModel<P>& f) {
        return os <<  f.operator VectorFunction<P>(); }

};

template<class P> inline CanonicalNumericType<P> evaluate(const ScalarFunctionModel<P>& f, const Vector<CanonicalNumericType<P>>& x) {
    return f(x); }
template<class P> inline CanonicalNumericType<P> unchecked_evaluate(const ScalarFunctionModel<P>& f, const Vector<CanonicalNumericType<P>>& x) {
    return f._ptr->_unchecked_evaluate(x); }


template<class P> inline CanonicalNumericType<P> unchecked_evaluate(const ScalarFunction<P>& f, const Vector<CanonicalNumericType<P>>& x) {
    ScalarFunctionModelInterface<P> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarFunctionModel<P>(*fptr),x); } else { return evaluate(f,x); } }
template<class P> inline Vector<CanonicalNumericType<P>> unchecked_evaluate(const VectorFunction<P>& f, const Vector<CanonicalNumericType<P>>& x) {
    VectorFunctionModelInterface<P> const* fptr = dynamic_cast<VectorFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorFunctionModel<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P> inline ScalarFunctionModel<P> unchecked_compose(const ScalarFunction<P>& f, const VectorFunctionModel<P>& g) {
    ScalarFunctionModelInterface<P> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<P>(*fptr),g); } else { return compose(f,g); } }
template<class P> inline VectorFunctionModel<P> unchecked_compose(const VectorFunction<P>& f, const VectorFunctionModel<P>& g) {
    VectorFunctionModelInterface<P> const* fptr = dynamic_cast<VectorFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<P>(*fptr),g); } else { return compose(f,g); } }

// FIXME: Declaration should be unneeded
VectorFunctionModel<ValidatedTag> unchecked_compose(const VectorFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g);
VectorFunctionModel<ValidatedTag> restrict(const VectorFunctionModel<ValidatedTag>& f, const ExactBoxType& g);
ScalarFunctionModel<ValidatedTag> restrict(const ScalarFunctionModel<ValidatedTag>& f, const ExactBoxType& g);



template<class P> VectorFunctionModelElement<P>::operator const ScalarFunctionModel<P> () const {
    return _p->get(_i); }



template<class P> inline FunctionModelBuilder<P>::FunctionModelBuilder(ScalarFunctionModel<P> const& f)
    : _prototype(f) { }
template<class P> inline ScalarFunctionModel<P> FunctionModelBuilder<P>::create(ScalarFunction<P> const& f) {
    return _prototype.create(f); }
template<class P> inline VectorFunctionModel<P> FunctionModelBuilder<P>::create(VectorFunction<P> const& f) {
    return _prototype.create(f); }
template<class P> inline CanonicalNumericType<P> FunctionModelBuilder<P>::create(Number<P> const& c) {
    return CanonicalNumericType<P>(c,_prototype.value().precision()); }


// Full output
template<class T> struct Representation { const T* pointer; Representation(const T& t) : pointer(&t) { } const T& reference() const { return *pointer; } };
template<class T> inline Representation<T> representation(const T& t) { return Representation<T>(t); }
template<class T> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { obj.reference().repr(os); return os; }

} // namespace Ariadne

#endif
