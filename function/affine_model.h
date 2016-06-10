/***************************************************************************
 *            affine_model.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file affine_model.h
 *  \brief Affine models defined on the unit box
 */

#ifndef ARIADNE_AFFINE_MODEL_H
#define ARIADNE_AFFINE_MODEL_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/declarations.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/matrix.h"
#include "algebra/algebra_operations.h"

namespace Ariadne {

template<class X> class Affine;
typedef Affine<Float64> FloatAffine;
typedef Affine<ExactIntervalType> IntervalAffine;
typedef Affine<ApproximateNumericType> ApproximateAffine;
typedef Affine<ValidatedNumericType> ValidatedAffine;

template<class X> class AffineModel;
typedef AffineModel<Float64> FloatAffineModel;
typedef AffineModel<ApproximateNumericType> ApproximateAffineModel;
typedef AffineModel<ValidatedNumericType> ValidatedAffineModel;

template<class P, class F> class TaylorModel;
typedef TaylorModel<ApproximateTag,Float64> ApproximateTaylorModel;
typedef TaylorModel<ValidatedTag,Float64> ValidatedTaylorModel;


AffineModel<ValidatedNumericType> affine_model(const Affine<ValidatedNumericType>& affine);
AffineModel<ValidatedNumericType> affine_model(const Affine<EffectiveNumericType>& affine);
AffineModel<ValidatedNumericType> affine_model(const TaylorModel<ValidatedTag,Float64>& taylor_model);
AffineModel<ValidatedNumericType> affine_model(const ExactBoxType& domain, const ValidatedScalarFunction& function);
Vector< AffineModel<ValidatedNumericType> > affine_models(const Vector< TaylorModel<ValidatedTag,Float64> >& taylor_models);
Vector< AffineModel<ValidatedNumericType> > affine_models(const ExactBoxType& domain, const ValidatedVectorFunction& function);

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x) \approx \sum_{i=0}^{n-1} a_i x_i + b\f$.
template<>
class AffineModel<ApproximateNumericType>
    : DispatchAlgebraOperations<AffineModel<ApproximateNumericType>,ApproximateNumericType>
{
  public:
    typedef Float64Approximation CoefficientType;

    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(Nat n) : _c(0.0), _g(n,0.0) { }
    explicit AffineModel(const ApproximateNumericType& c, const Covector<ApproximateNumericType>& g) : _c(c), _g(g) { }
    explicit AffineModel(ApproximateNumericType c, InitializerList<ApproximateNumericType> g) : _c(c), _g(g) { }

    AffineModel<ApproximateNumericType>& operator=(const ApproximateNumericType& c) {
        this->_c=c; for(Nat i=0; i!=this->_g.size(); ++i) { this->_g[i]=0.0; } return *this; }
    ApproximateNumericType& operator[](Nat i) { return this->_g[i]; }
    const ApproximateNumericType& operator[](Nat i) const { return this->_g[i]; }
    static AffineModel<ApproximateNumericType> constant(Nat n, const ApproximateNumericType& c) {
        return AffineModel<ApproximateNumericType>(c,Covector<ApproximateNumericType>(n,0.0)); }
    static AffineModel<ApproximateNumericType> variable(Nat n, Nat j) {
        return AffineModel<ApproximateNumericType>(0.0,Covector<ApproximateNumericType>::unit(n,j)); }

    const Covector<ApproximateNumericType>& a() const { return this->_g; }
    const ApproximateNumericType& b() const { return this->_c; }

    const Covector<ApproximateNumericType>& gradient() const { return this->_g; }
    const ApproximateNumericType& gradient(Nat i) const { return this->_g[i]; }
    const ApproximateNumericType& value() const { return this->_c; }

    Void resize(Nat n) { this->_g.resize(n); }

    Nat argument_size() const { return this->_g.size(); }
  private:
    ApproximateNumericType _c;
    Covector<ApproximateNumericType> _g;
};


//! An affine expression \f$f:[-1,+1]^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b \pm e\f$.
template<>
class AffineModel<ValidatedNumericType>
    : DispatchAlgebraOperations<AffineModel<ValidatedNumericType>,ValidatedNumericType>
{
  public:
    typedef Float64Value CoefficientType;
    typedef Float64Error ErrorType;

    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(Nat n) : _c(0.0), _g(n,Float64Value(0.0)), _e(0u) { }
    explicit AffineModel(const Float64Value& c, const Covector<Float64Value>& g, const Float64Error& e) : _c(c), _g(g), _e(e) { }
    explicit AffineModel(Float64Value c, InitializerList<Float64Value> g) : _c(c), _g(g), _e(0u) { }

    AffineModel<ValidatedNumericType>& operator=(const Float64Value& c) {
        this->_c=c; for(Nat i=0; i!=this->_g.size(); ++i) { this->_g[i]=0; } this->_e=0u; return *this; }
    static AffineModel<ValidatedNumericType> constant(Nat n, const Float64Value& c) {
        return AffineModel<ValidatedNumericType>(c,Covector<Float64Value>(n,0),Float64Error(0u)); }
    static AffineModel<ValidatedNumericType> variable(Nat n, Nat j) {
        return AffineModel<ValidatedNumericType>(0,Covector<Float64Value>::unit(n,j),0u); }


    const Covector<Float64Value>& a() const { return this->_g; }
    const Float64Value& b() const { return this->_c; }
    const Float64Error& e() const { return this->_e; }

    Float64Value& operator[](Nat i) { return this->_g[i]; }
    const Float64Value& operator[](Nat i) const { return this->_g[i]; }
    const Covector<Float64Value>& gradient() const { return this->_g; }
    const Float64Value& gradient(Nat i) const { return this->_g[i]; }
    const Float64Value& value() const { return this->_c; }
    const Float64Error& error() const { return this->_e; }

    Void set_value(const Float64Value& c) { _c=c; }
    Void set_gradient(Nat j, const Float64Value& g) { _g[j]=g; }
    Void set_error(const Float64Error& e) { _e=e; }

    Void resize(Nat n) { this->_g.resize(n); }

    Nat argument_size() const { return this->_g.size(); }
    template<class X> X evaluate(const Vector<X>& v) const {
        X r=v.zero_element()+static_cast<X>(this->_c);
        for(Nat i=0; i!=this->_g.size(); ++i) {
            r+=X(this->_g[i])*v[i]; }
        r+=ValidatedNumericType(-_e,+_e);
        return r;
    }

  private:
    Float64Value _c;
    Covector<Float64Value> _g;
    Float64Error _e;
};


//! \relates AffineModel
//! \brief Write to an output stream.
OutputStream& operator<<(OutputStream& os, const AffineModel<ValidatedNumericType>& f);

//! \relates AffineModel
//! \brief Create from a Taylor model.
AffineModel<ValidatedNumericType> affine_model(const TaylorModel<ValidatedTag,Float64>& tm);




} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
