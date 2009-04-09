/***************************************************************************
 *            transformation.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file transformation.h
 *  \brief Concrete transformations for low-level manipulation.
 */

#ifndef ARIADNE_TRANSFORMATION_H
#define ARIADNE_TRANSFORMATION_H

//#include <vararg>
#include <iosfwd>
#include <iostream>
#include <map>
#include "array.h"

#include "macros.h"

#include "vector.h"
#include "matrix.h"
#include "polynomial.h"

namespace Ariadne {

static const int SMOOTH=255;


struct ConstantTransformation 
{
    ConstantTransformation(const Vector<Float>& c, uint as)
        : _as(as), _c(c) { }
    const Vector<Float>& c() const { return _c; }
    const uint result_size() const { return _c.size(); }
    const uint argument_size() const { return _as; }
    const int smoothness() const { return SMOOTH; }
    template<class R, class A>
    void compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_c[i]; } }
  private:
    uint _as;
    Vector<Float> _c;
};

struct IdentityTransformation 
{
    IdentityTransformation(uint n)
        : _n(n) { }
    const uint result_size() const { return _n; }
    const uint argument_size() const { return _n; }
    const int smoothness() const { return SMOOTH; }
    template<class R, class A>
    void compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=x[i]; } }
  private:
    uint _n;
};

struct ProjectionTransformation 
{
    ProjectionTransformation(const array<uint>& p, uint as) 
        : _as(as), _p(p) 
    { for(uint i=0; i!=_p.size(); ++i) { ARIADNE_ASSERT(p[i]<as); } }
    ProjectionTransformation(uint rs, uint as, uint i) 
        : _as(as), _p(rs)
    { ARIADNE_ASSERT(rs+i<=as);
        for(uint j=0; j!=_p.size(); ++j) { _p[j]=i+j; } }
    ProjectionTransformation(const Range& rng, uint as) 
        : _as(as), _p(rng.size())
    { ARIADNE_ASSERT(rng.start()+rng.size()<=as);
        for(uint i=0; i!=_p.size(); ++i) { _p[i]=rng.start()+i; } }
    const array<uint>& p() const { return _p; }
    const uint operator[](uint i) const { return _p[i]; }
    const uint result_size() const { return _p.size(); }
    const uint argument_size() const { return _as; }
    const int smoothness() const { return SMOOTH; }
    template<class X> Vector<X> operator() (const Vector<X>& v) const {
        ARIADNE_ASSERT(v.size()==this->argument_size());
        Vector<X> r(result_size());    
        for(uint i=0; i!=result_size(); ++i) { r[i]=v[_p[i]]; }
        return r; }
    template<class R, class A>
    void compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=x[_p[i]]; } }
  private:
    uint _as;
    array<uint> _p;
};

struct ScalingTransformation 
{
    ScalingTransformation(const Vector<Float>& origin, 
                          const Vector<Float>& lengths)
        : _o(origin), _l(lengths) { ARIADNE_ASSERT(origin.size()==lengths.size()); }
    explicit ScalingTransformation(const Vector<Interval>& range)
        : _o(midpoint(range)), _l(range.size()) { for(uint i=0; i!=_l.size(); ++i) { _l[i]=range[i].radius(); } }
    const Vector<Float>& origin() const { return _o; }
    const Vector<Float>& lengths() const { return _l; }
    const uint result_size() const { return _l.size(); }
    const uint argument_size() const { return _l.size(); }
    const int smoothness() const { return SMOOTH; }
    template<class R, class A>
    void compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) {
            r[i]=_o[i]+_l[i]*x[i]; 
        }
    }
  private:
    Vector<Float> _o;
    Vector<Float> _l;
};


struct AffineExpression 
{
    AffineExpression(const Vector<Float>& g, const Float& c) : _c(c), _g(g.begin(),g.end()) { }
    AffineExpression(uint as, double c, double g0, ...) : _c(c), _g(as) {
        _g[0]=g0; va_list args; va_start(args,g0);
        for(uint i=1; i!=as; ++i) { _g[i]=va_arg(args,double); } }
    const uint argument_size() const { return _g.size()-1; }
    const int smoothness() const { return SMOOTH; }
    template<class R, class A>
    void compute(R& r, const A& x) const {
        r=_c; 
        for(uint j=0; j!=argument_size(); ++j) {
            r+=_g[j]*x[j];
        }
    }
  private:
    Float _c;
    array<Float> _g;
};

struct AffineTransformation 
{
    AffineTransformation(const Matrix<Float>& A, const Vector<Float>& b)
        : _A(A), _b(b) { ARIADNE_ASSERT(A.row_size()==b.size()); }
    const Matrix<Float>& A() const { return _A; }
    const Vector<Float>& b() const { return _b; }
    const uint result_size() const { return _A.row_size(); }
    const uint argument_size() const { return _A.column_size(); }
    const int smoothness() const { return SMOOTH; }
    template<class R, class A>
    void compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) {
            r[i]=_b[i]; 
            for(uint j=0; j!=argument_size(); ++j) {
                r[i]+=_A[i][j]*x[j];
            }
        }
    }
  private:
    Matrix<Float> _A;
    Vector<Float> _b;
};



struct PolynomialTransformation
{
    PolynomialTransformation();
    const uint result_size() const { return _p.size(); }
    const uint argument_size() const { return _p[0].argument_size(); }
    const int smoothness() const { return SMOOTH; }
    template<class R, class A>
    void compute(R& r, const A& x) const;
  private:
    std::vector< Polynomial<Float> > _p;
};



} // namespace Ariadne

#endif
