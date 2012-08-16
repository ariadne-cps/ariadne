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
typedef Affine<Interval> IntervalAffine;

template<class X> class AffineModel;
typedef AffineModel<Float> FloatAffineModel;
typedef AffineModel<Interval> IntervalAffineModel;

template<class F, class B> class Constraint;
typedef Constraint<IntervalAffineModel,Float> IntervalAffineModelConstraint;

template<class X> class ScalarFunction;
template<class X> class TaylorModel;

AffineModel<Interval> affine_model(const Affine<Interval>& affine);
AffineModel<Interval> affine_model(const Affine<Real>& affine);
AffineModel<Interval> affine_model(const TaylorModel<Interval>& taylor_model);
AffineModel<Interval> affine_model(const IntervalVector& domain, const ScalarFunction<Interval>& function);
Vector< AffineModel<Interval> > affine_models(const Vector< TaylorModel<Interval> >& taylor_models);
Vector< AffineModel<Interval> > affine_models(const IntervalVector& domain, const VectorFunction<Interval>& function);

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x) \approx \sum_{i=0}^{n-1} a_i x_i + b\f$.
template<>
class AffineModel<Float>
{
  public:
    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(uint n) : _c(0.0), _g(n,0.0) { }
    explicit AffineModel(const Float& c, const Vector<Float>& g) : _c(c), _g(g) { }
    explicit AffineModel(Float c, std::initializer_list<Float> g) : _c(c), _g(g) { }

    AffineModel<Float>& operator=(const Float& c) {
        this->_c=c; for(uint i=0; i!=this->_g.size(); ++i) { this->_g[i]=0.0; } return *this; }
    Float& operator[](uint i) { return this->_g[i]; }
    const Float& operator[](uint i) const { return this->_g[i]; }
    static AffineModel<Float> constant(uint n, const Float& c) {
        return AffineModel<Float>(c,Vector<Float>(n,0.0)); }
    static AffineModel<Float> variable(uint n, uint j) {
        return AffineModel<Float>(0.0,Vector<Float>::unit(n,j)); }

    const Vector<Float>& a() const { return this->_g; }
    const Float& b() const { return this->_c; }

    const Vector<Float>& gradient() const { return this->_g; }
    const Float& gradient(uint i) const { return this->_g[i]; }
    const Float& value() const { return this->_c; }

    void resize(uint n) { this->_g.resize(n); }

    uint argument_size() const { return this->_g.size(); }
  private:
    Float _c;
    Vector<Float> _g;
};

//! \relates AffineModel
//! \brief Negation of an affine model.
AffineModel<Float> operator-(const AffineModel<Float>& f);
//! \relates AffineModel
//! \brief Addition of two affine models.
AffineModel<Float> operator+(const AffineModel<Float>& f1, const AffineModel<Float>& f2);
//! \relates AffineModel
//! \brief Subtraction of two affine models.
AffineModel<Float> operator-(const AffineModel<Float>& f1, const AffineModel<Float>& f2);
//! \relates AffineModel
//! \brief Multiplication of two affine models.
AffineModel<Float> operator*(const AffineModel<Float>& f1, const AffineModel<Float>& f2);
//! \relates AffineModel
//! \brief Addition of a constant to an affine model.
AffineModel<Float>& operator+=(AffineModel<Float>& f1, const Float& c2);
//! \relates AffineModel
//! \brief Scalar multiplication of an affine model.
AffineModel<Float>& operator*=(AffineModel<Float>& f1, const Float& c2);

//! \relates AffineModel
//! \brief Write to an output stream.
std::ostream& operator<<(std::ostream& os, const AffineModel<Float>& f);


//! An affine expression \f$f:[-1,+1]^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b \pm e\f$.
template<>
class AffineModel<Interval>
{
  public:
    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(uint n) : _c(0.0), _g(n,0.0), _e(0.0) { }
    explicit AffineModel(const Float& c, const Vector<Float>& g, const Float& e) : _c(c), _g(g), _e(e) { }
    explicit AffineModel(Float c, std::initializer_list<Float> g) : _c(c), _g(g), _e(0.0) { }

    AffineModel<Interval>& operator=(const Float& c) {
        this->_c=c; for(uint i=0; i!=this->_g.size(); ++i) { this->_g[i]=0.0; } this->_e=0.0; return *this; }
    static AffineModel<Interval> constant(uint n, const Float& c) {
        return AffineModel<Interval>(c,Vector<Float>(n,0.0),0.0); }
    static AffineModel<Interval> variable(uint n, uint j) {
        return AffineModel<Interval>(0.0,Vector<Float>::unit(n,j),0.0); }


    const Vector<Float>& a() const { return this->_g; }
    const Float& b() const { return this->_c; }
    const Float& e() const { return this->_e; }

    Float& operator[](uint i) { return this->_g[i]; }
    const Float& operator[](uint i) const { return this->_g[i]; }
    const Vector<Float>& gradient() const { return this->_g; }
    const Float& gradient(uint i) const { return this->_g[i]; }
    const Float& value() const { return this->_c; }
    const Float& error() const { return this->_e; }

    void set_value(const Float& c) { _c=c; }
    void set_gradient(uint j, const Float& g) { _g[j]=g; }
    void set_error(const Float& e) { _e=e; }

    void resize(uint n) { this->_g.resize(n); }

    uint argument_size() const { return this->_g.size(); }
    template<class X> X evaluate(const Vector<X>& v) const {
        X r=v.zero_element()+Interval(this->_c);
        for(uint i=0; i!=this->_g.size(); ++i) {
            r+=Interval(this->_g[i])*v[i]; }
        r+=Interval(-_e,+_e);
        return r;
    }

  private:
    Float _c;
    Vector<Float> _g;
    Float _e;
};

//! \relates AffineModel
//! \brief Negation of an affine model.
AffineModel<Interval> operator-(const AffineModel<Interval>& f);
//! \relates AffineModel
//! \brief Addition of two affine models.
AffineModel<Interval> operator+(const AffineModel<Interval>& f1, const AffineModel<Interval>& f2);
//! \relates AffineModel
//! \brief Subtraction of two affine models.
AffineModel<Interval> operator-(const AffineModel<Interval>& f1, const AffineModel<Interval>& f2);
//! \relates AffineModel
//! \brief Multiplication of two affine models.
AffineModel<Interval> operator*(const AffineModel<Interval>& f1, const AffineModel<Interval>& f2);
//! \relates AffineModel
//! \brief Addition of a constant to an affine model.
AffineModel<Interval>& operator+=(AffineModel<Interval>& f1, const Interval& c2);
//! \relates AffineModel
//! \brief Scalar multiplication of an affine model.
AffineModel<Interval>& operator*=(AffineModel<Interval>& f1, const Interval& c2);

//! \relates AffineModel \brief Scalar addition to an affine model.
AffineModel<Interval> operator+(const Interval& c1, const AffineModel<Interval>& f2);
AffineModel<Interval> operator+(const AffineModel<Interval>& f1, const Interval& c2);
//! \relates AffineModel \brief Subtraction of an affine model from a scalar.
AffineModel<Interval> operator-(const Interval& c1, const AffineModel<Interval>& f2);
//! \relates AffineModel \brief Subtraction of a scalar from an affine model.
AffineModel<Interval> operator-(const AffineModel<Interval>& f1, const Interval& c2);
//! \relates AffineModel \brief Scalar multiplication of an affine model.
AffineModel<Interval> operator*(const Interval& c1, const AffineModel<Interval>& f2);

//! \relates AffineModel
//! \brief Write to an output stream.
std::ostream& operator<<(std::ostream& os, const AffineModel<Interval>& f);

//! \relates AffineModel
//! \brief Create from a Taylor model.
AffineModel<Interval> affine_model(const TaylorModel<Interval>& tm);




} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
