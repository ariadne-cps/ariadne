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

#include "macros.h"
#include "pointer.h"

#include "numeric.h"
#include "vector.h"
#include "matrix.h"

namespace Ariadne {

template<class X> class Affine;
typedef Affine<Float> FloatAffine;
typedef Affine<ExactInterval> IntervalAffine;
typedef Affine<ApproximateNumberType> ApproximateAffine;
typedef Affine<ValidatedNumberType> ValidatedAffine;

template<class X> class AffineModel;
typedef AffineModel<Float> FloatAffineModel;
typedef AffineModel<ExactInterval> IntervalAffineModel;
typedef AffineModel<ApproximateTag> ApproximateAffineModel;
typedef AffineModel<ValidatedTag> ValidatedAffineModel;

template<class X> class ScalarFunction;
template<class X> class TaylorModel;
typedef TaylorModel<ApproximateTag> ApproximateTaylorModel;
typedef TaylorModel<ValidatedTag> ValidatedTaylorModel;


AffineModel<ValidatedTag> affine_model(const Affine<ValidatedNumberType>& affine);
AffineModel<ValidatedTag> affine_model(const Affine<EffectiveNumberType>& affine);
AffineModel<ValidatedTag> affine_model(const TaylorModel<ValidatedTag>& taylor_model);
AffineModel<ValidatedTag> affine_model(const ExactBox& domain, const ValidatedScalarFunction& function);
Vector< AffineModel<ValidatedTag> > affine_models(const Vector< TaylorModel<ValidatedTag> >& taylor_models);
Vector< AffineModel<ValidatedTag> > affine_models(const ExactBox& domain, const ValidatedVectorFunction& function);

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x) \approx \sum_{i=0}^{n-1} a_i x_i + b\f$.
template<>
class AffineModel<ApproximateTag>
{
  public:
    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(uint n) : _c(0.0), _g(n,0.0) { }
    explicit AffineModel(const ApproximateNumberType& c, const Vector<ApproximateNumberType>& g) : _c(c), _g(g) { }
    explicit AffineModel(ApproximateNumberType c, std::initializer_list<ApproximateNumberType> g) : _c(c), _g(g) { }

    AffineModel<ApproximateTag>& operator=(const ApproximateNumberType& c) {
        this->_c=c; for(uint i=0; i!=this->_g.size(); ++i) { this->_g[i]=0.0; } return *this; }
    ApproximateNumberType& operator[](uint i) { return this->_g[i]; }
    const ApproximateNumberType& operator[](uint i) const { return this->_g[i]; }
    static AffineModel<ApproximateTag> constant(uint n, const ApproximateNumberType& c) {
        return AffineModel<ApproximateTag>(c,Vector<ApproximateNumberType>(n,0.0)); }
    static AffineModel<ApproximateTag> variable(uint n, uint j) {
        return AffineModel<ApproximateTag>(0.0,Vector<ApproximateNumberType>::unit(n,j)); }

    const Vector<ApproximateNumberType>& a() const { return this->_g; }
    const ApproximateNumberType& b() const { return this->_c; }

    const Vector<ApproximateNumberType>& gradient() const { return this->_g; }
    const ApproximateNumberType& gradient(uint i) const { return this->_g[i]; }
    const ApproximateNumberType& value() const { return this->_c; }

    void resize(uint n) { this->_g.resize(n); }

    uint argument_size() const { return this->_g.size(); }
  private:
    ApproximateNumberType _c;
    Vector<ApproximateNumberType> _g;
};

//! \relates AffineModel
//! \brief Negation of an affine model.
AffineModel<ApproximateTag> operator-(const AffineModel<ApproximateTag>& f);
//! \relates AffineModel
//! \brief Addition of two affine models.
AffineModel<ApproximateTag> operator+(const AffineModel<ApproximateTag>& f1, const AffineModel<ApproximateTag>& f2);
//! \relates AffineModel
//! \brief Subtraction of two affine models.
AffineModel<ApproximateTag> operator-(const AffineModel<ApproximateTag>& f1, const AffineModel<ApproximateTag>& f2);
//! \relates AffineModel
//! \brief Multiplication of two affine models.
AffineModel<ApproximateTag> operator*(const AffineModel<ApproximateTag>& f1, const AffineModel<ApproximateTag>& f2);
//! \relates AffineModel
//! \brief Addition of a constant to an affine model.
AffineModel<ApproximateTag>& operator+=(AffineModel<ApproximateTag>& f1, const ApproximateNumberType& c2);
//! \relates AffineModel
//! \brief Scalar multiplication of an affine model.
AffineModel<ApproximateTag>& operator*=(AffineModel<ApproximateTag>& f1, const ApproximateNumberType& c2);

//! \relates AffineModel
//! \brief Write to an output stream.
std::ostream& operator<<(std::ostream& os, const AffineModel<ApproximateTag>& f);


//! An affine expression \f$f:[-1,+1]^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b \pm e\f$.
template<>
class AffineModel<ValidatedTag>
{
  public:
    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(uint n) : _c(0.0), _g(n,ExactFloatType(0.0)), _e(0.0) { }
    explicit AffineModel(const ExactFloatType& c, const Vector<ExactFloatType>& g, const ErrorFloatType& e) : _c(c), _g(g), _e(e) { }
    explicit AffineModel(ExactFloatType c, std::initializer_list<ExactFloatType> g) : _c(c), _g(g), _e(0.0) { }

    AffineModel<ValidatedTag>& operator=(const ExactFloatType& c) {
        this->_c=c; for(uint i=0; i!=this->_g.size(); ++i) { this->_g[i]=0; } this->_e=0; return *this; }
    static AffineModel<ValidatedTag> constant(uint n, const ExactFloatType& c) {
        return AffineModel<ValidatedTag>(c,Vector<ExactFloatType>(n,0),ErrorFloatType(0)); }
    static AffineModel<ValidatedTag> variable(uint n, uint j) {
        return AffineModel<ValidatedTag>(0,Vector<ExactFloatType>::unit(n,j),0); }


    const Vector<ExactFloatType>& a() const { return this->_g; }
    const ExactFloatType& b() const { return this->_c; }
    const ErrorFloatType& e() const { return this->_e; }

    ExactFloatType& operator[](uint i) { return this->_g[i]; }
    const ExactFloatType& operator[](uint i) const { return this->_g[i]; }
    const Vector<ExactFloatType>& gradient() const { return this->_g; }
    const ExactFloatType& gradient(uint i) const { return this->_g[i]; }
    const ExactFloatType& value() const { return this->_c; }
    const ErrorFloatType& error() const { return this->_e; }

    void set_value(const ExactFloatType& c) { _c=c; }
    void set_gradient(uint j, const ExactFloatType& g) { _g[j]=g; }
    void set_error(const ErrorFloatType& e) { _e=e; }

    void resize(uint n) { this->_g.resize(n); }

    uint argument_size() const { return this->_g.size(); }
    template<class X> X evaluate(const Vector<X>& v) const {
        X r=v.zero_element()+static_cast<X>(this->_c);
        for(uint i=0; i!=this->_g.size(); ++i) {
            r+=X(this->_g[i])*v[i]; }
        r+=ValidatedNumberType(-_e,+_e);
        return r;
    }

  private:
    ExactFloatType _c;
    Vector<ExactFloatType> _g;
    ErrorFloatType _e;
};

//! \relates AffineModel
//! \brief Negation of an affine model.
AffineModel<ValidatedTag> operator-(const AffineModel<ValidatedTag>& f);
//! \relates AffineModel
//! \brief Addition of two affine models.
AffineModel<ValidatedTag> operator+(const AffineModel<ValidatedTag>& f1, const AffineModel<ValidatedTag>& f2);
//! \relates AffineModel
//! \brief Subtraction of two affine models.
AffineModel<ValidatedTag> operator-(const AffineModel<ValidatedTag>& f1, const AffineModel<ValidatedTag>& f2);
//! \relates AffineModel
//! \brief Multiplication of two affine models.
AffineModel<ValidatedTag> operator*(const AffineModel<ValidatedTag>& f1, const AffineModel<ValidatedTag>& f2);
//! \relates AffineModel
//! \brief Addition of a constant to an affine model.
AffineModel<ValidatedTag>& operator+=(AffineModel<ValidatedTag>& f1, const ValidatedNumberType& c2);
//! \relates AffineModel
//! \brief Scalar multiplication of an affine model.
AffineModel<ValidatedTag>& operator*=(AffineModel<ValidatedTag>& f1, const ValidatedNumberType& c2);

//! \relates AffineModel \brief Scalar addition to an affine model.
AffineModel<ValidatedTag> operator+(const ValidatedNumberType& c1, const AffineModel<ValidatedTag>& f2);
AffineModel<ValidatedTag> operator+(const AffineModel<ValidatedTag>& f1, const ValidatedNumberType& c2);
//! \relates AffineModel \brief Subtraction of an affine model from a scalar.
AffineModel<ValidatedTag> operator-(const ValidatedNumberType& c1, const AffineModel<ValidatedTag>& f2);
//! \relates AffineModel \brief Subtraction of a scalar from an affine model.
AffineModel<ValidatedTag> operator-(const AffineModel<ValidatedTag>& f1, const ValidatedNumberType& c2);
//! \relates AffineModel \brief Scalar multiplication of an affine model.
AffineModel<ValidatedTag> operator*(const ValidatedNumberType& c1, const AffineModel<ValidatedTag>& f2);

//! \relates AffineModel
//! \brief Write to an output stream.
std::ostream& operator<<(std::ostream& os, const AffineModel<ValidatedTag>& f);

//! \relates AffineModel
//! \brief Create from a Taylor model.
AffineModel<ValidatedTag> affine_model(const TaylorModel<ValidatedTag>& tm);




} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
