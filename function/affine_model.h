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

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/matrix.h"

namespace Ariadne {

template<class X> class Affine;
typedef Affine<Float> FloatAffine;
typedef Affine<ExactInterval> IntervalAffine;
typedef Affine<ApproximateNumber> ApproximateAffine;
typedef Affine<ValidatedNumber> ValidatedAffine;

template<class X> class AffineModel;
typedef AffineModel<Float> FloatAffineModel;
typedef AffineModel<ApproximateNumber> ApproximateAffineModel;
typedef AffineModel<ValidatedNumber> ValidatedAffineModel;

template<class X> class ScalarFunction;
template<class X> class TaylorModel;
typedef TaylorModel<ApproximateNumber> ApproximateTaylorModel;
typedef TaylorModel<ValidatedNumber> ValidatedTaylorModel;


AffineModel<ValidatedNumber> affine_model(const Affine<ValidatedNumber>& affine);
AffineModel<ValidatedNumber> affine_model(const Affine<EffectiveNumber>& affine);
AffineModel<ValidatedNumber> affine_model(const TaylorModel<ValidatedNumber>& taylor_model);
AffineModel<ValidatedNumber> affine_model(const ExactBox& domain, const ValidatedScalarFunction& function);
Vector< AffineModel<ValidatedNumber> > affine_models(const Vector< TaylorModel<ValidatedNumber> >& taylor_models);
Vector< AffineModel<ValidatedNumber> > affine_models(const ExactBox& domain, const ValidatedVectorFunction& function);

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x) \approx \sum_{i=0}^{n-1} a_i x_i + b\f$.
template<>
class AffineModel<ApproximateNumber>
{
  public:
    typedef ApproximateFloat CoefficientType;

    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(Nat n) : _c(0.0), _g(n,0.0) { }
    explicit AffineModel(const ApproximateNumber& c, const Covector<ApproximateNumber>& g) : _c(c), _g(g) { }
    explicit AffineModel(ApproximateNumber c, InitializerList<ApproximateNumber> g) : _c(c), _g(g) { }

    AffineModel<ApproximateNumber>& operator=(const ApproximateNumber& c) {
        this->_c=c; for(Nat i=0; i!=this->_g.size(); ++i) { this->_g[i]=0.0; } return *this; }
    ApproximateNumber& operator[](Nat i) { return this->_g[i]; }
    const ApproximateNumber& operator[](Nat i) const { return this->_g[i]; }
    static AffineModel<ApproximateNumber> constant(Nat n, const ApproximateNumber& c) {
        return AffineModel<ApproximateNumber>(c,Covector<ApproximateNumber>(n,0.0)); }
    static AffineModel<ApproximateNumber> variable(Nat n, Nat j) {
        return AffineModel<ApproximateNumber>(0.0,Covector<ApproximateNumber>::unit(n,j)); }

    const Covector<ApproximateNumber>& a() const { return this->_g; }
    const ApproximateNumber& b() const { return this->_c; }

    const Covector<ApproximateNumber>& gradient() const { return this->_g; }
    const ApproximateNumber& gradient(Nat i) const { return this->_g[i]; }
    const ApproximateNumber& value() const { return this->_c; }

    Void resize(Nat n) { this->_g.resize(n); }

    Nat argument_size() const { return this->_g.size(); }
  private:
    ApproximateNumber _c;
    Covector<ApproximateNumber> _g;
};

//! \relates AffineModel
//! \brief Negation of an affine model.
AffineModel<ApproximateNumber> operator-(const AffineModel<ApproximateNumber>& f);
//! \relates AffineModel
//! \brief Addition of two affine models.
AffineModel<ApproximateNumber> operator+(const AffineModel<ApproximateNumber>& f1, const AffineModel<ApproximateNumber>& f2);
//! \relates AffineModel
//! \brief Subtraction of two affine models.
AffineModel<ApproximateNumber> operator-(const AffineModel<ApproximateNumber>& f1, const AffineModel<ApproximateNumber>& f2);
//! \relates AffineModel
//! \brief Multiplication of two affine models.
AffineModel<ApproximateNumber> operator*(const AffineModel<ApproximateNumber>& f1, const AffineModel<ApproximateNumber>& f2);
//! \relates AffineModel
//! \brief Addition of a constant to an affine model.
AffineModel<ApproximateNumber>& operator+=(AffineModel<ApproximateNumber>& f1, const ApproximateNumber& c2);
//! \relates AffineModel
//! \brief Scalar multiplication of an affine model.
AffineModel<ApproximateNumber>& operator*=(AffineModel<ApproximateNumber>& f1, const ApproximateNumber& c2);

//! \relates AffineModel
//! \brief Write to an output stream.
OutputStream& operator<<(OutputStream& os, const AffineModel<ApproximateNumber>& f);


//! An affine expression \f$f:[-1,+1]^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b \pm e\f$.
template<>
class AffineModel<ValidatedNumber>
{
  public:
    typedef ExactFloat CoefficientType;
    typedef ErrorFloat ErrorType;

    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(Nat n) : _c(0.0), _g(n,ExactFloat(0.0)), _e(0.0) { }
    explicit AffineModel(const ExactFloat& c, const Covector<ExactFloat>& g, const ErrorFloat& e) : _c(c), _g(g), _e(e) { }
    explicit AffineModel(ExactFloat c, InitializerList<ExactFloat> g) : _c(c), _g(g), _e(0.0) { }

    AffineModel<ValidatedNumber>& operator=(const ExactFloat& c) {
        this->_c=c; for(Nat i=0; i!=this->_g.size(); ++i) { this->_g[i]=0; } this->_e=0; return *this; }
    static AffineModel<ValidatedNumber> constant(Nat n, const ExactFloat& c) {
        return AffineModel<ValidatedNumber>(c,Covector<ExactFloat>(n,0),ErrorFloat(0)); }
    static AffineModel<ValidatedNumber> variable(Nat n, Nat j) {
        return AffineModel<ValidatedNumber>(0,Covector<ExactFloat>::unit(n,j),0); }


    const Covector<ExactFloat>& a() const { return this->_g; }
    const ExactFloat& b() const { return this->_c; }
    const ErrorFloat& e() const { return this->_e; }

    ExactFloat& operator[](Nat i) { return this->_g[i]; }
    const ExactFloat& operator[](Nat i) const { return this->_g[i]; }
    const Covector<ExactFloat>& gradient() const { return this->_g; }
    const ExactFloat& gradient(Nat i) const { return this->_g[i]; }
    const ExactFloat& value() const { return this->_c; }
    const ErrorFloat& error() const { return this->_e; }

    Void set_value(const ExactFloat& c) { _c=c; }
    Void set_gradient(Nat j, const ExactFloat& g) { _g[j]=g; }
    Void set_error(const ErrorFloat& e) { _e=e; }

    Void resize(Nat n) { this->_g.resize(n); }

    Nat argument_size() const { return this->_g.size(); }
    template<class X> X evaluate(const Vector<X>& v) const {
        X r=v.zero_element()+static_cast<X>(this->_c);
        for(Nat i=0; i!=this->_g.size(); ++i) {
            r+=X(this->_g[i])*v[i]; }
        r+=ValidatedNumber(-_e,+_e);
        return r;
    }

  private:
    ExactFloat _c;
    Covector<ExactFloat> _g;
    ErrorFloat _e;
};

//! \relates AffineModel
//! \brief Negation of an affine model.
AffineModel<ValidatedNumber> operator-(const AffineModel<ValidatedNumber>& f);
//! \relates AffineModel
//! \brief Addition of two affine models.
AffineModel<ValidatedNumber> operator+(const AffineModel<ValidatedNumber>& f1, const AffineModel<ValidatedNumber>& f2);
//! \relates AffineModel
//! \brief Subtraction of two affine models.
AffineModel<ValidatedNumber> operator-(const AffineModel<ValidatedNumber>& f1, const AffineModel<ValidatedNumber>& f2);
//! \relates AffineModel
//! \brief Multiplication of two affine models.
AffineModel<ValidatedNumber> operator*(const AffineModel<ValidatedNumber>& f1, const AffineModel<ValidatedNumber>& f2);
//! \relates AffineModel
//! \brief Addition of a constant to an affine model.
AffineModel<ValidatedNumber>& operator+=(AffineModel<ValidatedNumber>& f1, const ValidatedNumber& c2);
//! \relates AffineModel
//! \brief Scalar multiplication of an affine model.
AffineModel<ValidatedNumber>& operator*=(AffineModel<ValidatedNumber>& f1, const ValidatedNumber& c2);

//! \relates AffineModel \brief Scalar addition to an affine model.
AffineModel<ValidatedNumber> operator+(const ValidatedNumber& c1, const AffineModel<ValidatedNumber>& f2);
AffineModel<ValidatedNumber> operator+(const AffineModel<ValidatedNumber>& f1, const ValidatedNumber& c2);
//! \relates AffineModel \brief Subtraction of an affine model from a scalar.
AffineModel<ValidatedNumber> operator-(const ValidatedNumber& c1, const AffineModel<ValidatedNumber>& f2);
//! \relates AffineModel \brief Subtraction of a scalar from an affine model.
AffineModel<ValidatedNumber> operator-(const AffineModel<ValidatedNumber>& f1, const ValidatedNumber& c2);
//! \relates AffineModel \brief Scalar multiplication of an affine model.
AffineModel<ValidatedNumber> operator*(const ValidatedNumber& c1, const AffineModel<ValidatedNumber>& f2);

//! \relates AffineModel
//! \brief Write to an output stream.
OutputStream& operator<<(OutputStream& os, const AffineModel<ValidatedNumber>& f);

//! \relates AffineModel
//! \brief Create from a Taylor model.
AffineModel<ValidatedNumber> affine_model(const TaylorModel<ValidatedNumber>& tm);




} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
