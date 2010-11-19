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

template<class X> class Expression;
template<class X> class FormulaNode;

//! \brief A formula defining a real function.
//!
//! The Formula class is implemented as a directed acyclic graph, with
//! each node being an atomic operation.
template<class X>
class Formula {
  public:
    ~Formula() { --_root->count; if(_root->count==0) { delete _root; } }
    explicit Formula(const Expression<Real>& e, const Map<std::string,uint> s);
    explicit Formula() : _root(new FormulaNode<X>(0.0)) { ++_root->count; }
    explicit Formula(const FormulaNode<X>* ptr) : _root(ptr) { assert(ptr); ++_root->count; }
    Formula(const Formula<X>& f) : _root(f._root) { ++_root->count; }
    Formula(const X& c) : _root(new FormulaNode<X>(c)) { ++_root->count; }
    Formula<X>& operator=(const Formula<X>& f) { if(this != &f) { --_root->count; if(_root->count==0) { delete _root; } _root=f._root; ++_root->count; } return *this; }
    Formula<X>& operator=(const X& c) { --_root->count; if(_root->count==0) { delete _root; } _root=new FormulaNode<X>(c); ++_root->count; return *this; }
    static Formula<X> constant(const X& c);
    static Formula<X> constant(double c);
    static Formula<X> coordinate(uint j);
    static Vector< Formula<X> > identity(uint n);
    const FormulaNode<X>* ptr() const { return _root; }
  public:
    const FormulaNode<X>* _root;
    //boost::intrusive_ptr< const FormulaNode<X> > _root;
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

template<class X> inline void intrusive_ptr_add_ref(const FormulaNode<X>* nodeptr) { ++nodeptr->count; }
template<class X> inline void intrusive_ptr_release(const FormulaNode<X>* nodeptr) { --nodeptr->count; if(nodeptr->count==0) { delete nodeptr; } }

template<class X> inline Formula<X> make_formula(Operator op, const Formula<X>& f) {
    return Formula<X>(new FormulaNode<X>(op,f.ptr())); }
template<class X> inline Formula<X> make_formula(Operator op, const Formula<X>& f, int n) {
    return Formula<X>(new FormulaNode<X>(op,f.ptr(),n)); }
template<class X> inline Formula<X> make_formula(Operator op, const Formula<X>& f1, const Formula<X>& f2) {
    return Formula<X>(new FormulaNode<X>(op,f1.ptr(),f2.ptr())); }


template<class X> inline Formula<X> Formula<X>::constant(const X& c) { return Formula<X>(new FormulaNode<X>(c)); }
template<class X> inline Formula<X> Formula<X>::constant(double c) { return Formula<X>(new FormulaNode<X>(c)); }
template<class X> inline Formula<X> Formula<X>::coordinate(uint j) { return Formula<X>(new FormulaNode<X>(j)); }
template<class X> inline Vector< Formula<X> > Formula<X>::identity(uint n) {
    Vector< Formula<X> > r(n); for(uint i=0; i!=n; ++i) { r[i]=Formula<X>::coordinate(i); } return r; }

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

template<class X> inline Formula<X> operator+(double c, Formula<X> f) { return Formula<X>::constant(c) + f; }
template<class X> inline Formula<X> operator-(double c, Formula<X> f) { return Formula<X>::constant(c) - f; }
template<class X> inline Formula<X> operator*(double c, Formula<X> f) { return Formula<X>::constant(c) * f; }
template<class X> inline Formula<X> operator/(double c, Formula<X> f) { return Formula<X>::constant(c) / f; }
template<class X> inline Formula<X> operator+(Formula<X> f, double c) { return f + Formula<X>::constant(c); }
template<class X> inline Formula<X> operator-(Formula<X> f, double c) { return f - Formula<X>::constant(c); }
template<class X> inline Formula<X> operator*(Formula<X> f, double c) { return f * Formula<X>::constant(c); }
template<class X> inline Formula<X> operator/(Formula<X> f, double c) { return f / Formula<X>::constant(c); }

template<class X> inline Formula<X> operator+(Real c, Formula<X> f) { return Formula<X>::constant(c) + f; }
template<class X> inline Formula<X> operator-(Real c, Formula<X> f) { return Formula<X>::constant(c) - f; }
template<class X> inline Formula<X> operator*(Real c, Formula<X> f) { return Formula<X>::constant(c) * f; }
template<class X> inline Formula<X> operator/(Real c, Formula<X> f) { return Formula<X>::constant(c) / f; }
template<class X> inline Formula<X> operator+(Formula<X> f, Real c) { return f + Formula<X>::constant(c); }
template<class X> inline Formula<X> operator-(Formula<X> f, Real c) { return f - Formula<X>::constant(c); }
template<class X> inline Formula<X> operator*(Formula<X> f, Real c) { return f * Formula<X>::constant(c); }
template<class X> inline Formula<X> operator/(Formula<X> f, Real c) { return f / Formula<X>::constant(c); }

template<class X> inline Formula<X>& operator+=(Formula<X>& f1, const Formula<X>& f2) { Formula<X> r=f1+f2; return f1=r; }
template<class X> inline Formula<X>& operator*=(Formula<X>& f1, const Formula<X>& f2) { Formula<X> r=f1*f2; return f1=r; }
template<class X> inline Formula<X>& operator+=(Formula<X>& f, const Real& c) { Formula<X> r=f+Formula<X>::constant(c); return f=r; }
template<class X> inline Formula<X>& operator*=(Formula<X>& f, const Real& c) { Formula<X> r=f*Formula<X>::constant(c); return f=r; }

template<class X> inline Formula<X> operator*(Formula<X> f, int c) {
    if(c==0) { return Formula<X>::constant(0.0); }
    else { return f * Formula<X>::constant(c); }
}

template<class X> inline Formula<X> operator+(Formula<X> f, int c) {
    if(f.ptr()->op==CNST) { return Formula<X>::constant(*f.ptr()->val+c); }
    else { return f + Formula<X>::constant(c); }
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const FormulaNode<X>* f) {
    switch(f->op) {
        //case CNST: return os << std::fixed << std::setprecision(4) << f->val;
        case CNST:
            if(*f->val==0.0) { return os << 0.0; } if(abs(*f->val)<1e-4) { os << std::fixed << *f->val; } else { os << *f->val; } return os;
        case IND:
            return os << "x[" << f->ind << "]";
        case ADD:
            return os << f->arg1 << '+' << f->arg2;
        case SUB:
            os << f->arg1 << '-';
            switch(f->arg2->op) { case ADD: case SUB: os << '(' << f->arg2 << ')'; break; default: os << f->arg2; }
            return os;
        case MUL:
            switch(f->arg1->op) { case ADD: case SUB: case DIV: os << '(' << f->arg1 << ')'; break; default: os << f->arg1; }
            os << '*';
            switch(f->arg2->op) { case ADD: case SUB: os << '(' << f->arg2 << ')'; break; default: os << f->arg2; }
            return os;
        case DIV:
            switch(f->arg1->op) { case ADD: case SUB: case DIV: os << '(' << f->arg1 << ')'; break; default: os << f->arg1; }
            os << '/';
            switch(f->arg2->op) { case ADD: case SUB: case MUL: case DIV: os << '(' << f->arg2 << ')'; break; default: os << f->arg2; }
            return os;
        case POW:
            return os << "pow" << '(' << f->arg << ',' << f->np << ')';
        case EXP:
            return os << "exp(" << f->arg << ")";
        default: return os << f->op << "(" << f->arg << ")";
    }
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const Formula<X>& f) {
    if(f._root) { return os << f.ptr(); } else { return  os << "NULL"; }
}

template<class X, class T> inline T evaluate(const FormulaNode<X>* f, const T* a) {
    switch(f->op) {
        //case CNST: return os << std::fixed << std::setprecision(4) << f->val;
        case CNST: return static_cast<T>(*f->val);
        case IND: return a[f->ind];
        case ADD: return evaluate<T>(f->arg1,a) + evaluate<T>(f->arg2,a);
        case SUB: return evaluate<T>(f->arg1,a) - evaluate<T>(f->arg2,a);
        case MUL: return evaluate<T>(f->arg1,a) * evaluate<T>(f->arg2,a);
        case DIV: return evaluate<T>(f->arg1,a) / evaluate<T>(f->arg2,a);
        case NEG: return neg(evaluate<T>(f->arg,a));
        case REC: return rec(evaluate<T>(f->arg,a));
        case SQR: return sqr(evaluate<T>(f->arg,a));
        case SQRT: return sqrt(evaluate<T>(f->arg,a));
        case EXP: return exp(evaluate<T>(f->arg,a));
        case LOG: return log(evaluate<T>(f->arg,a));
        case SIN: return sin(evaluate<T>(f->arg,a));
        case COS: return cos(evaluate<T>(f->arg,a));
        case TAN: return tan(evaluate<T>(f->arg,a));
        default: assert(false);
    }
}

template<class X, class T> inline T& evaluate(const FormulaNode<X>* f, const T* a, Map<const void*,T>& tmp) {
    if(tmp.has_key(f)) { return tmp[f]; }
    switch(f->op) { // Can't use simple evaluate (above) as we need to pass the cache to subformulae
        case CNST: return tmp[f]=f->val;
        case IND: return tmp[f]=a[f->ind];
        case ADD: return tmp[f]=evaluate<T>(f->arg1,a,tmp) + evaluate<T>(f->arg2,a,tmp);
        case SUB: return tmp[f]=evaluate<T>(f->arg1,a,tmp) - evaluate<T>(f->arg2,a,tmp);
        case MUL: return tmp[f]=evaluate<T>(f->arg1,a,tmp) * evaluate<T>(f->arg2,a,tmp);
        case DIV: return tmp[f]=evaluate<T>(f->arg1,a,tmp) / evaluate<T>(f->arg2,a,tmp);
        case NEG: return tmp[f]=neg(evaluate<T>(f->arg,a,tmp));
        case REC: return tmp[f]=rec(evaluate<T>(f->arg,a,tmp));
        case SQRT: return tmp[f]=sqrt(evaluate<T>(f->arg,a,tmp));
        case EXP: return tmp[f]=exp(evaluate<T>(f->arg,a,tmp));
        case LOG: return tmp[f]=log(evaluate<T>(f->arg,a,tmp));
        case SIN: return tmp[f]=sin(evaluate<T>(f->arg,a,tmp));
        case COS: return tmp[f]=cos(evaluate<T>(f->arg,a,tmp));
        case TAN: return tmp[f]=tan(evaluate<T>(f->arg,a,tmp));
        default: assert(false);
    }
}

template<class X, class T> T evaluate(const Formula<X>& f, const Vector<T>& v) {
    return evaluate(f.ptr(),&v[0]);
}

template<class X, class T> T map_evaluate(const Formula<X>& f, const Vector<T>& v) {
    Map<const void*,T> tmp;
    return evaluate(f.ptr(),&v[0],tmp);
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
