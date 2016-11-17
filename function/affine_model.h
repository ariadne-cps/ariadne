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
#include "algebra/operations.h"

namespace Ariadne {

template<class X> class Affine;
typedef Affine<Float64> FloatAffine;
typedef Affine<ExactIntervalType> IntervalAffine;
typedef Affine<ApproximateNumericType> ApproximateAffine;
typedef Affine<ValidatedNumericType> ValidatedAffine;

template<class P, class F> class AffineModel;
typedef AffineModel<ApproximateTag,Float64> ApproximateAffineModel;
typedef AffineModel<ValidatedTag,Float64> ValidatedAffineModel;

template<class P, class F> class TaylorModel;
typedef TaylorModel<ApproximateTag,Float64> ApproximateTaylorModel;
typedef TaylorModel<ValidatedTag,Float64> ValidatedTaylorModel;


AffineModel<ValidatedTag,Float64> affine_model(const Affine<Float64Bounds>& affine);
AffineModel<ValidatedTag,Float64> affine_model(const Affine<ValidatedNumber>& affine);
AffineModel<ValidatedTag,Float64> affine_model(const TaylorModel<ValidatedTag,Float64>& taylor_model);
AffineModel<ValidatedTag,Float64> affine_model(const ExactBoxType& domain, const ValidatedScalarFunction& function);
Vector< AffineModel<ValidatedTag,Float64> > affine_models(const Vector< TaylorModel<ValidatedTag,Float64> >& taylor_models);
Vector< AffineModel<ValidatedTag,Float64> > affine_models(const ExactBoxType& domain, const ValidatedVectorFunction& function);

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x) \approx \sum_{i=0}^{n-1} a_i x_i + b\f$.
template<class F>
class AffineModel<ApproximateTag,F>
    : DispatchAlgebraOperations<AffineModel<ApproximateTag,F>,ApproximateNumericType>
{
    typedef typename F::PrecisionType PR;
  public:
    typedef PR PrecisionType;
    typedef FloatApproximation<PR> CoefficientType;

    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(SizeType n, PrecisionType prec) : _c(0,prec), _g(n,CoefficientType(0,prec)) { }
    explicit AffineModel(SizeType n, const CoefficientType& c) : _c(c), _g(n,nul(c)) { }
    explicit AffineModel(const CoefficientType& c, const Covector<CoefficientType>& g) : _c(c), _g(g) { }
    explicit AffineModel(CoefficientType c, InitializerList<CoefficientType> g) : _c(c), _g(g) { }

    AffineModel<ApproximateTag,F>& operator=(const CoefficientType& c) {
        this->_c=c; for(SizeType i=0; i!=this->_g.size(); ++i) { this->_g[i]=0.0; } return *this; }
    CoefficientType& operator[](SizeType i) { return this->_g[i]; }
    const CoefficientType& operator[](SizeType i) const { return this->_g[i]; }
    static AffineModel<ApproximateTag,F> constant(SizeType n, const CoefficientType& c) {
        return AffineModel<ApproximateTag,F>(c,Covector<CoefficientType>(n,CoefficientType(0.0))); }
    static AffineModel<ApproximateTag,F> variable(SizeType n, SizeType j) {
        return AffineModel<ApproximateTag,F>(CoefficientType(0.0),Covector<CoefficientType>::unit(n,j)); }

    const Covector<CoefficientType>& a() const { return this->_g; }
    const CoefficientType& b() const { return this->_c; }

    const Covector<CoefficientType>& gradient() const { return this->_g; }
    const CoefficientType& gradient(SizeType i) const { return this->_g[i]; }
    const CoefficientType& value() const { return this->_c; }

    Void resize(SizeType n) { this->_g.resize(n); }

    SizeType argument_size() const { return this->_g.size(); }
  private:
    CoefficientType _c;
    Covector<CoefficientType> _g;
};


//! An affine expression \f$f:[-1,+1]^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b \pm e\f$.
template<class F>
class AffineModel<ValidatedTag,F>
    : DispatchAlgebraOperations<AffineModel<ValidatedTag,F>,ValidatedNumericType>
{
    typedef typename ValidatedNumericType::PrecisionType PR;
  public:
    typedef PR PrecisionType;
    typedef Float64Value CoefficientType;
    typedef Float64Error ErrorType;
    typedef AffineModel<ValidatedTag,F> AffineModelType;

    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(SizeType n, PrecisionType prec) : _c(0,prec), _g(n,CoefficientType(0,prec)), _e(0u,prec) { }
    explicit AffineModel(SizeType n, const CoefficientType& c) : _c(c), _g(n,nul(c)), _e(nul(c)) { }
    explicit AffineModel(const CoefficientType& c, const Covector<CoefficientType>& g, const ErrorType& e) : _c(c), _g(g), _e(e) { }
    explicit AffineModel(CoefficientType c, InitializerList<CoefficientType> g) : _c(c), _g(g), _e(0u) { }

    AffineModelType& operator=(const CoefficientType& c) {
        this->_c=c; for(SizeType i=0; i!=this->_g.size(); ++i) { this->_g[i]=0; } this->_e=0u; return *this; }
    static AffineModelType constant(SizeType n, const CoefficientType& c) {
        PrecisionType prec=c.precision(); return AffineModelType(c,Covector<CoefficientType>(n,CoefficientType(0,prec)),ErrorType(0u,prec)); }
    static AffineModelType variable(SizeType n, SizeType j, PrecisionType prec) {
        return AffineModelType(CoefficientType(0,prec),Covector<CoefficientType>::unit(n,j),ErrorType(0u,prec)); }


    PrecisionType precision() const { return this->_c.precision(); }
    const Covector<CoefficientType>& a() const { return this->_g; }
    const CoefficientType& b() const { return this->_c; }
    const ErrorType& e() const { return this->_e; }

    CoefficientType& operator[](SizeType i) { return this->_g[i]; }
    const CoefficientType& operator[](SizeType i) const { return this->_g[i]; }
    const Covector<CoefficientType>& gradient() const { return this->_g; }
    const CoefficientType& gradient(SizeType i) const { return this->_g[i]; }
    const CoefficientType& value() const { return this->_c; }
    const ErrorType& error() const { return this->_e; }

    Void set_value(const CoefficientType& c) { _c=c; }
    Void set_gradient(SizeType j, const CoefficientType& g) { _g[j]=g; }
    Void set_error(const ErrorType& e) { _e=e; }
    Void set_error(Nat m) { _e=m; }

    Void resize(SizeType n) { this->_g.resize(n); }

    SizeType argument_size() const { return this->_g.size(); }
    template<class X> X evaluate(const Vector<X>& v) const {
        X r=v.zero_element()+static_cast<X>(this->_c);
        for(SizeType i=0; i!=this->_g.size(); ++i) {
            r+=X(this->_g[i])*v[i]; }
        r+=ValidatedNumericType(-_e,+_e);
        return r;
    }

  private:
    CoefficientType _c;
    Covector<CoefficientType> _g;
    ErrorType _e;
};


//! \relates AffineModel
//! \brief Write to an output stream.
OutputStream& operator<<(OutputStream& os, const AffineModel<ValidatedTag,Float64>& f);

//! \relates AffineModel
//! \brief Create from a Taylor model.
AffineModel<ValidatedTag,Float64> affine_model(const TaylorModel<ValidatedTag,Float64>& tm);




} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
