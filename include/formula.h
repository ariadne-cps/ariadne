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
#include <boost/concept_check.hpp>

namespace Ariadne {

template<class X> class FormulaNode;

//! \brief A formula defining a real function.
//!
//! The Formula class is implemented as a directed acyclic graph, with
//! each node being an atomic operation.
template<class X>
class Formula {
    typedef uint I;
  public:
    typedef typename X::NumericType NumericType;
    typedef I IndexType;
  public:
    ~Formula() { --_root->count; if(_root->count==0) { delete _root; } }
    explicit Formula() : _root(new FormulaNode<X>(0.0)) { ++_root->count; }
    explicit Formula(const FormulaNode<X>* fptr) : _root(fptr) { assert(fptr); ++_root->count; }
    Formula(const Formula<X>& f) : _root(f._root) { ++_root->count; }
    Formula(const X& c) : _root(new FormulaNode<X>(c)) { ++_root->count; }
    Formula<X>& operator=(const Formula<X>& f) { if(this != &f) { --_root->count; if(_root->count==0) { delete _root; } _root=f._root; ++_root->count; } return *this; }
    Formula<X>& operator=(const X& c) { --_root->count; if(_root->count==0) { delete _root; } _root=new FormulaNode<X>(c); ++_root->count; return *this; }
    static Formula<X> constant(const X& c);
    static Formula<X> constant(double c);
    static Formula<X> coordinate(const I& j);
    static Vector< Formula<X> > identity(uint n);
    const FormulaNode<X>* raw_ptr() const { return _root; }
  public:
    Formula<X> _arg() const { return Formula<X>(_root->arg); }
    Formula<X> _arg1() const { return Formula<X>(_root->arg1); }
    Formula<X> _arg2() const { return Formula<X>(_root->arg2); }
  public:
    const FormulaNode<X>* _root;
    //boost::intrusive_fptr< const FormulaNode<X> > _root;
};

//! \related Procedure \brief Evaluate a function \a f defined by a formula.
template<class X, class T> T evaluate(const Formula<X>& f, const Vector<T>& v);

// Class for which an object x produces coordinate \f$x_j\f$ when calling \c x[j].
struct Coordinate { Formula<Real> operator[](uint j) { return Formula<Real>::coordinate(j); } };

template<class X>
struct FormulaNode {
    ~FormulaNode() {
        switch(op) {
            case CNST: delete val; case IND: break;
            case POW: --arg->count; if(arg->count==0) { delete arg; } break;
            default: if(arg2) { --arg2->count; if(arg2->count==0) { delete arg2; } } --arg1->count; if(arg1->count==0) { delete arg1; } } }
    explicit FormulaNode(Operator o, const FormulaNode<X>* a) : count(0u), op(o), arg(a), arg2(0) { ++arg->count; }
    explicit FormulaNode(Operator o, const FormulaNode<X>* a, int n) : count(0u), op(o), arg(a), np(n) { ARIADNE_ASSERT(o==POW); ++arg->count; }
    explicit FormulaNode(Operator o, const FormulaNode<X>* a1, const FormulaNode<X>* a2) : count(0u), op(o), arg1(a1), arg2(a2) { ARIADNE_ASSERT(a2); ++arg1->count; ++arg2->count; }
    explicit FormulaNode(const X& x) : count(0u), op(CNST), val(new X(x)) { }
    explicit FormulaNode(double x) : count(0u), op(CNST), val(new X(x)) { }
    explicit FormulaNode(uint i) : count(0u), op(IND), ind(i) { }
    mutable uint count;
    Operator op;
    union {
        uint ind; X* val;
        struct { const FormulaNode* arg; int np; };
        struct { const FormulaNode* arg1; const FormulaNode* arg2; };
    };
};

template<class X> inline void intrusive_fptr_add_ref(const FormulaNode<X>* nodefptr) { ++nodefptr->count; }
template<class X> inline void intrusive_fptr_release(const FormulaNode<X>* nodefptr) { --nodefptr->count; if(nodefptr->count==0) { delete nodefptr; } }

template<class X> inline Formula<X> make_formula(Operator op, const Formula<X>& f) {
    return Formula<X>(new FormulaNode<X>(op,f.raw_ptr())); }
template<class X> inline Formula<X> make_formula(Operator op, const Formula<X>& f, int n) {
    return Formula<X>(new FormulaNode<X>(op,f.raw_ptr(),n)); }
template<class X> inline Formula<X> make_formula(Operator op, const Formula<X>& f1, const Formula<X>& f2) {
    return Formula<X>(new FormulaNode<X>(op,f1.raw_ptr(),f2.raw_ptr())); }


template<class X> inline Formula<X> Formula<X>::constant(const X& c) { return Formula<X>(new FormulaNode<X>(c)); }
template<class X> inline Formula<X> Formula<X>::constant(double c) { return Formula<X>(new FormulaNode<X>(c)); }
template<class X> inline Formula<X> Formula<X>::coordinate(const I& j) { return Formula<X>(new FormulaNode<X>(j)); }
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

template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator+(Formula<X> f, R c) { return f + Formula<X>::constant(numeric_cast<X>(c)); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator-(Formula<X> f, R c) { return f - Formula<X>::constant(numeric_cast<X>(c)); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator*(Formula<X> f, R c) { return f * Formula<X>::constant(numeric_cast<X>(c)); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator/(Formula<X> f, R c) { return f / Formula<X>::constant(numeric_cast<X>(c)); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator+(R c, Formula<X> f) { return Formula<X>::constant(numeric_cast<X>(c)) + f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator-(R c, Formula<X> f) { return Formula<X>::constant(numeric_cast<X>(c)) - f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator*(R c, Formula<X> f) { return Formula<X>::constant(numeric_cast<X>(c)) * f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator/(R c, Formula<X> f) { return Formula<X>::constant(numeric_cast<X>(c)) / f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type& operator+=(Formula<X>& f, const R& c) { return f+=Formula<X>::constant(numeric_cast<X>(c)); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type& operator*=(Formula<X>& f, const R& c) { return f*=Formula<X>::constant(numeric_cast<X>(c)); }

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

template<class X> inline std::ostream& operator<<(std::ostream& os, const FormulaNode<X>* fptr) {
    switch(fptr->op) {
        //case CNST: return os << std::fixed << std::setprecision(4) << fptr->val;
        case CNST:
            if(*fptr->val==0.0) { return os << 0.0; } if(abs(*fptr->val)<1e-4) { os << std::fixed << *fptr->val; } else { os << *fptr->val; } return os;
        case IND:
            return os << "x[" << fptr->ind << "]";
        case ADD:
            return os << fptr->arg1 << '+' << fptr->arg2;
        case SUB:
            os << fptr->arg1 << '-';
            switch(fptr->arg2->op) { case ADD: case SUB: os << '(' << fptr->arg2 << ')'; break; default: os << fptr->arg2; }
            return os;
        case MUL:
            switch(fptr->arg1->op) { case ADD: case SUB: case DIV: os << '(' << fptr->arg1 << ')'; break; default: os << fptr->arg1; }
            os << '*';
            switch(fptr->arg2->op) { case ADD: case SUB: os << '(' << fptr->arg2 << ')'; break; default: os << fptr->arg2; }
            return os;
        case DIV:
            switch(fptr->arg1->op) { case ADD: case SUB: case DIV: os << '(' << fptr->arg1 << ')'; break; default: os << fptr->arg1; }
            os << '/';
            switch(fptr->arg2->op) { case ADD: case SUB: case MUL: case DIV: os << '(' << fptr->arg2 << ')'; break; default: os << fptr->arg2; }
            return os;
        case POW:
            return os << "pow" << '(' << fptr->arg << ',' << fptr->np << ')';
        case EXP:
            return os << "exp(" << fptr->arg << ")";
        default: return os << fptr->op << "(" << fptr->arg << ")";
    }
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const Formula<X>& f) {
    if(f._root) { return os << f.raw_ptr(); } else { return  os << "NULL"; }
}

template<class X, class T> T _evaluate(const FormulaNode<X>* fptr, const T* a, const T& z) {
    switch(fptr->op) {
        //case CNST: return os << std::fixed << std::setprecision(4) << fptr->val;
        case CNST: return z+*fptr->val;
        case IND: return a[fptr->ind];
        case ADD: return _evaluate(fptr->arg1,a,z) + _evaluate(fptr->arg2,a,z);
        case SUB: return _evaluate(fptr->arg1,a,z) - _evaluate(fptr->arg2,a,z);
        case MUL: return _evaluate(fptr->arg1,a,z) * _evaluate(fptr->arg2,a,z);
        case DIV: return _evaluate(fptr->arg1,a,z) / _evaluate(fptr->arg2,a,z);
        case NEG: return neg(_evaluate(fptr->arg,a,z));
        case REC: return rec(_evaluate(fptr->arg,a,z));
        case SQR: return sqr(_evaluate(fptr->arg,a,z));
        case SQRT: return sqrt(_evaluate(fptr->arg,a,z));
        case EXP: return exp(_evaluate(fptr->arg,a,z));
        case LOG: return log(_evaluate(fptr->arg,a,z));
        case SIN: return sin(_evaluate(fptr->arg,a,z));
        case COS: return cos(_evaluate(fptr->arg,a,z));
        case TAN: return tan(_evaluate(fptr->arg,a,z));
        default: assert(false);
    }
}

template<class X, class T> T& _evaluate(const FormulaNode<X>* fptr, const T* a, const T& z, Map<const void*,T>& tmp) {
    if(tmp.has_key(fptr)) { return tmp[fptr]; }
    switch(fptr->op) { // Can't use simple evaluate (above) as we need to pass the cache to subformulae
        case CNST: return tmp[fptr]=z+(*fptr->val);
        case IND: return tmp[fptr]=a[fptr->ind];
        case ADD: return tmp[fptr]=_evaluate(fptr->arg1,a,z,tmp) + _evaluate(fptr->arg2,a,z,tmp);
        case SUB: return tmp[fptr]=_evaluate(fptr->arg1,a,z,tmp) - _evaluate(fptr->arg2,a,z,tmp);
        case MUL: return tmp[fptr]=_evaluate(fptr->arg1,a,z,tmp) * _evaluate(fptr->arg2,a,z,tmp);
        case DIV: return tmp[fptr]=_evaluate(fptr->arg1,a,z,tmp) / _evaluate(fptr->arg2,a,z,tmp);
        case NEG: return tmp[fptr]=neg(_evaluate(fptr->arg,a,z,tmp));
        case REC: return tmp[fptr]=rec(_evaluate(fptr->arg,a,z,tmp));
        case SQRT: return tmp[fptr]=sqrt(_evaluate(fptr->arg,a,z,tmp));
        case EXP: return tmp[fptr]=exp(_evaluate(fptr->arg,a,z,tmp));
        case LOG: return tmp[fptr]=log(_evaluate(fptr->arg,a,z,tmp));
        case SIN: return tmp[fptr]=sin(_evaluate(fptr->arg,a,z,tmp));
        case COS: return tmp[fptr]=cos(_evaluate(fptr->arg,a,z,tmp));
        case TAN: return tmp[fptr]=tan(_evaluate(fptr->arg,a,z,tmp));
        default: assert(false);
    }
}

template<class X, class T> T evaluate(const Formula<X>& f, const Vector<T>& v) {
    assert(v.size()!=0);
    const T z=v[0]*0.0;
    return _evaluate(f.raw_ptr(),&v[0],z);
}

template<class X, class T> T map_evaluate(const Formula<X>& f, const Vector<T>& v) {
    assert(v.size()!=0);
    const T z=v[0]*0.0;
    Map<const void*,T> tmp;
    return _evaluate(f.raw_ptr(),&v[0],z,tmp);
}

template<class X, class T> Vector<T> map_evaluate(const Vector< Formula<X> >& f, const Vector<T>& v) {
    assert(v.size()!=0);
    Vector<T> r(f.size());
    const T z=v[0]*0.0;
    Map<const void*,T> tmp;
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=_evaluate(f[i].raw_ptr(),&v[0],z,tmp);
    }
    return r;
}

template<class X> Formula<X> derivative(const Formula<X>& f, uint j)
{
    switch(f.raw_ptr()->op) {
        case CNST:
            return Formula<X>(0.0);
        case IND:
            if(f.raw_ptr()->ind==j) { return Formula<Real>(1.0); }
            else { return Formula<Real>(0.0); }
        case ADD:
            return derivative(f._arg1(),j)+derivative(f._arg2(),j);
        case SUB:
            return derivative(f._arg1(),j)-derivative(f._arg2(),j);
        case MUL:
            return f._arg1()*derivative(f._arg2(),j)+derivative(f._arg1(),j)*f._arg2();
        case DIV:
            return derivative(f._arg1() * rec(f._arg2()),j);
        case NEG:
            return  - derivative(f._arg(),j);
        case REC:
            return  - derivative(f._arg(),j) * rec(sqr(f._arg()));
        case SQR:
            return 2 * derivative(f._arg(),j) * f._arg();
        case EXP:
            return derivative(f._arg(),j) * f._arg();
        case LOG:
            return derivative(f._arg(),j) * rec(f._arg());
        case SIN:
            return derivative(f._arg(),j) * cos(f._arg());
        case COS:
            return -derivative(f._arg(),j) * sin(f._arg());
        case TAN:
            return derivative(f._arg(),j) * (1-sqr(f._arg()));
        default:
            ARIADNE_THROW(std::runtime_error,"derivative(Formual<X>,Nat)",
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
