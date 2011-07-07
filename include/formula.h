/***************************************************************************
 *            formula.h
 *
 *  Copyright 2008-9  Pieter Collins
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


/*! \file formula.h
 *  \brief Formulae over variables
 */

#ifndef ARIADNE_FORMULA_H
#define ARIADNE_FORMULA_H

#include <cstdarg>
#include <iostream>
#include <iomanip>
#include <string>


#include "macros.h"
#include "pointer.h"
#include "container.h"
#include "stlio.h"

#include "operators.h"
#include "expansion.h"
#include "numeric.h"
#include "vector.h"
#include "expression.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

//! \brief A formula defining a real function.
//!
//! The Formula class is implemented as a directed acyclic graph, with
//! each node being an atomic operation.
template<class X>
class Formula
    : public Expression<X>
{
    typedef Index I;
  public:
    typedef X NumericType;
    typedef I IndexType;
  public:
    Formula() : Expression<X>(static_cast<X>(0)) { }
    Formula(const Expression<X>& fh) : Expression<X>(fh) { }
    Formula(double c) : Expression<X>(static_cast<X>(c)) { }
    Formula(const X& c) : Expression<X>(c) { }
    Formula(const I& i) : Expression<X>(i) { }
    Formula<X>& operator=(double c) { return *this=Expression<X>(static_cast<X>(c)); }
    Formula<X>& operator=(const X& c) { return *this=Expression<X>(c); }
    static Formula<X> constant(const X& c) { return Expression<X>(c); }
    static Formula<X> constant(double c) { return constant(static_cast<X>(c)); }
    static Formula<X> coordinate(const Nat& j) { return Expression<X>(Index(j)); };
    static Vector< Formula<X> > identity(uint n);

    const Expression<X>& handle() const { return *this; }
};

// Class for which an object x produces coordinate \f$x_j\f$ when calling \c x[j].
struct Coordinate { Formula<Real> operator[](uint j) { return Formula<Real>::coordinate(j); } };

template<class X, class R> inline Formula<X> make_formula(const R& c) {
    return Expression<X>(numeric_cast<X>(c)); }
template<class X> inline Formula<X> make_formula(OperatorCode op, const Formula<X>& f) {
    return Expression<X>(op,f.handle()); }
template<class X> inline Formula<X> make_formula(OperatorCode op, const Formula<X>& f1, const Formula<X>& f2) {
    return make_expression<X>(op,f1.handle(),f2.handle()); }
template<class X> inline Formula<X> make_formula(OperatorCode op, const Formula<X>& f, int n) {
    return make_expression<X>(op,f.handle(),n); }

template<class X> inline Vector< Formula<X> > Formula<X>::identity(uint n) {
    Vector< Formula<X> > r(n); for(uint i=0; i!=n; ++i) { r[i]=Formula<X>::coordinate(i); } return r; }

template<class X> inline Formula<X>& operator+=(Formula<X>& f1, const Formula<X>& f2) { Formula<X> r=f1+f2; return f1=r; }
template<class X> inline Formula<X>& operator*=(Formula<X>& f1, const Formula<X>& f2) { Formula<X> r=f1*f2; return f1=r; }

template<class X> inline Formula<X> operator+(const Formula<X>& f) { return make_formula(POS,f); }
template<class X> inline Formula<X> operator-(const Formula<X>& f) { return make_formula(NEG,f); }
template<class X> inline Formula<X> operator+(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(ADD,f1,f2); }
template<class X> inline Formula<X> operator-(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(SUB,f1,f2); }
template<class X> inline Formula<X> operator*(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(MUL,f1,f2); }
template<class X> inline Formula<X> operator/(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(DIV,f1,f2); }

template<class X> inline Formula<X> neg(const Formula<X>& f) { return make_formula(NEG,f); }
template<class X> inline Formula<X> rec(const Formula<X>& f) { return make_formula(REC,f); }
template<class X> inline Formula<X> sqr(const Formula<X>& f) { return make_formula(SQR,f); }
template<class X> inline Formula<X> pow(const Formula<X>& f, int n) { return make_formula(POW,f,n); }
template<class X> inline Formula<X> sqrt(const Formula<X>& f) { return make_formula(SQRT,f); }
template<class X> inline Formula<X> exp(const Formula<X>& f) { return make_formula(EXP,f); }
template<class X> inline Formula<X> log(const Formula<X>& f) { return make_formula(LOG,f); }
template<class X> inline Formula<X> sin(const Formula<X>& f) { return make_formula(SIN,f); }
template<class X> inline Formula<X> cos(const Formula<X>& f) { return make_formula(COS,f); }
template<class X> inline Formula<X> tan(const Formula<X>& f) { return make_formula(TAN,f); }
template<class X> inline Formula<X> atan(const Formula<X>& f) { return make_formula(ATAN,f); }

template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator+(Formula<X> f, R c) { return f + make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator-(Formula<X> f, R c) { return f - make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator*(Formula<X> f, R c) { return f * make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator/(Formula<X> f, R c) { return f / make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator+(R c, Formula<X> f) { return make_formula<X>(c) + f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator-(R c, Formula<X> f) { return make_formula<X>(c) - f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator*(R c, Formula<X> f) { return make_formula<X>(c) * f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator/(R c, Formula<X> f) { return make_formula<X>(c) / f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type& operator+=(Formula<X>& f, const R& c) { return f+=make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type& operator*=(Formula<X>& f, const R& c) { return f*=make_formula<X>(c); }

/*
template<class X> inline Formula<X> operator*(Formula<X> f, int c) {
    if(c==0) { return Formula<X>::constant(0.0); }
    else { return f * Formula<X>::constant(c); }
}

template<class X> inline Formula<X> operator+(Formula<X> f, int c) {
    if(f.raw_ptr()->op==CNST) { return Formula<X>::constant(*f.raw_ptr()->val+c); }
    else { return f + Formula<X>::constant(c); }
}
*/


template<class X, class T> T evaluate(const Formula<X>& f, const Vector<T>& v) {
    return evaluate(f.handle(),v);
}

template<class X, class T> T map_evaluate(const Formula<X>& f, const Vector<T>& v) {
    return cached_evaluate(f.handle(),v,Map<const Void*,T>());
}

template<class X, class T> Vector<T> map_evaluate(const Vector< Formula<X> >& f, const Vector<T>& v) {
    assert(v.size()!=0);
    Vector<T> r(f.size());
    Map<const void*,T> cache;
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=cached_evaluate(f[i],v,cache);
    }
    return r;
}

template<class X> Formula<X> derivative(const Formula<X>& f, uint j)
{
    switch(f.op()) {
        case CNST:
            return Formula<X>(0.0);
        case IND:
            if(f.ind()==j) { return Formula<Real>(1.0); }
            else { return Formula<Real>(0.0); }
        case ADD:
            return derivative(f.arg1(),j)+derivative(f.arg2(),j);
        case SUB:
            return derivative(f.arg1(),j)-derivative(f.arg2(),j);
        case MUL:
            return f.arg1()*derivative(f.arg2(),j)+derivative(f.arg1(),j)*f.arg2();
        case DIV:
            return derivative(f.arg1() * rec(f.arg2()),j);
        case NEG:
            return  - derivative(f.arg(),j);
        case REC:
            return  - derivative(f.arg(),j) * rec(sqr(f.arg()));
        case SQR:
            return static_cast<X>(2) * derivative(f.arg(),j) * f.arg();
        case EXP:
            return derivative(f.arg(),j) * f.arg();
        case LOG:
            return derivative(f.arg(),j) * rec(f.arg());
        case SIN:
            return derivative(f.arg(),j) * cos(f.arg());
        case COS:
            return -derivative(f.arg(),j) * sin(f.arg());
        case TAN:
            return derivative(f.arg(),j) * (static_cast<X>(1)-sqr(f.arg()));
        default:
            ARIADNE_THROW(std::runtime_error,"derivative(Formual<X>)",
                          "Cannot compute derivative of "<<f<<"\n");
    }
}


//! \ingroup FunctionModule
//! \brief Convert a power-series expansion into a formula using a version of Horner's rule.
//! See J. M. Pena and T. Sauer, "On the multivariate Horner scheme", SIAM J. Numer. Anal. 37(4) 1186-1197, 2000.
//!
template<class X> Formula<X> formula(const Expansion<X>& e)
{
    Vector< Formula<X> > identity(e.argument_size());
    for(uint i=0; i!=identity.size(); ++i) { identity[i]=Formula<X>::coordinate(i); }
    return horner_evaluate(e,identity);
}

} // namespace Ariadne


#endif // ARIADNE_FORMULA_H
