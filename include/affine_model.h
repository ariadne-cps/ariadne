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

#include "vector.h"
#include "matrix.h"

namespace Ariadne {

template<class X> class AffineModel;


//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
template<>
class AffineModel<Interval>
{
  public:
    explicit AffineModel() : _c(), _g() { }
    explicit AffineModel(uint n) : _c(0.0), _g(n,0.0), _e(0.0) { }
    explicit AffineModel(const Float& c, const Vector<Float>& g, const Float& e) : _c(c), _g(g), _e(e) { }
    explicit AffineModel(uint as, double c, double g0, ...) : _c(static_cast<X>(c)), _g(as), _e(0.0) {
        _g[0]=static_cast<X>(g0); va_list args; va_start(args,g0);
        for(uint i=1; i!=as; ++i) { _g[i]=static_cast<X>(va_arg(args,double)); } }

    AffineModel<Interval>& operator=(const Float& c) {
        this->_c=c; for(uint i=0; i!=this->_g.size(); ++i) { this->_g[i]=0.0; } this->_e=0.0; return *this; }
    static AffineModel<Interval> constant(uint n, const Float& c) {
        return AffineModel<Interval>(c,Vector<Float>(n,0.0)); }
    static AffineModel<Interval> variable(uint n, uint j) {
        return AffineModel<Interval>(0.0,Vector<Float>::unit(n,j)); }


    const Vector<Float>& a() const { return this->_g; }
    const Float& b() const { return this->_c; }

    const Vector<Float>& gradient() const { return this->_g; }
    const Float& gradient(uint i) const { return this->_g[i]; }
    const Float& value() const { return this->_c; }
    const Float& error() const { return this->_e; }

    uint argument_size() const { return this->_g.size(); }
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

//! \relates AffineModel
//! \brief Write to an output stream.
std::ostream& operator<<(std::ostream& os, const AffineModel<Interval>& f);





} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
