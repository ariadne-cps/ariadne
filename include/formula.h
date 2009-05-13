/***************************************************************************
 *            formula.h
 *
 *  Copyright 2009  Pieter Collins
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
 *  \brief Formulae over variables.
 */

#ifndef ARIADNE_FORMULA_H
#define ARIADNE_FORMULA_H

#include <cstdarg>
#include <iostream>
#include <string>
#include <utility>
#include "expression_interface.h"
#include "function_interface.h"

#include "macros.h"
#include "pointer.h"

#include "polynomial.h"
#include "function.h"

namespace Ariadne {

template<class F> class Formula;
template<class V, class F> class Assignment;
template<class F1,class Op, class F2> class Comparison;

class GtrZero {}; class LessZero {}; class Gtr { }; class Less { };
//struct Add { template<class R,class A1,class A2> R evaluate(const A1& a1,const A2& a2) const { return a1+a2; } };
//struct Sub { template<class R,class A1,class A2> R evaluate(const A1& a1,const A2& a2) const { return a1-a2; } };
//struct Mul { template<class R,class A1,class A2> R evaluate(const A1& a1,const A2& a2) const { return a1*a2; } };
//struct Div { template<class R,class A1,class A2> R evaluate(const A1& a1,const A2& a2) const { return a1/a2; } };
struct Add { template<class T> T operator()(const T& a1,const T& a2) const { return a1+a2; } };
struct Sub { template<class T> T operator()(const T& a1,const T& a2) const { return a1-a2; } };
struct Mul { template<class T> T operator()(const T& a1,const T& a2) const { return a1*a2; } };
struct Div { template<class T> T operator()(const T& a1,const T& a2) const { return a1/a2; } };

template<class F> struct Formula {
    F& operator()() { return static_cast<F&>(*this); }
    const F& operator()() const { return static_cast<const F&>(*this); }
};

template<class C> struct Constant : Formula< Constant<C> > {
    Constant(const C& v) : val(v) { }
    C val;
};

struct Variable : Formula<Variable> {
  public:
    Variable(std::string name) : _name_ptr(new std::string(name)) { }
    const std::string& name() const { return *this->_name_ptr; }
    template<class F> Assignment<Variable,F> operator=(const Formula<F>&) const;
    //bool operator==(const Variable& v) const { return name()==v.name(); }
    bool operator==(const Variable& v) const { return _name_ptr==v._name_ptr; }
    // Use default constructor
  private:
    shared_ptr<const std::string> _name_ptr;
};

struct NextVariable {
    NextVariable(Variable v) : var(v) { }
    std::string name() const { return "next("+var.name()+")"; }
    template<class F> Assignment<NextVariable,F> operator=(const Formula<F>& f) const;;
    Variable var;
};

struct DottedVariable {
    DottedVariable(Variable v) : var(v) { }
    std::string name() const { return "dot("+var.name()+")"; }
    template<class F> Assignment<DottedVariable,F> operator=(const Formula<F>& f) const;
    Variable var;
};

NextVariable next(Variable v) { return NextVariable(v); }
DottedVariable dot(Variable v) { return DottedVariable(v); }

struct Space {
    Space() : _variables() { }
    unsigned int size() const { return _variables.size(); }
    unsigned int dimension() const { return _variables.size(); }
    const Variable& operator[](unsigned int i) const { return _variables.at(i); }
    const Variable& variable(unsigned int i) const { return _variables.at(i); }
    unsigned int index(const Variable& v) const;
    Space& operator,(const Variable& v) { _variables.push_back(v); }
  private:
    std::vector<Variable> _variables;
};

inline std::ostream& operator<<(std::ostream& os, const Space& s) {
    os<<"Space"; for(unsigned int i=0; i!=s.size(); ++i) { os<<(i==0?"(":",")<<s.variable(i).name(); } return os<<")"; }

inline unsigned int Space::index(const Variable& v) const { for(unsigned int i=0; i!=_variables.size(); ++i) {
        if(_variables[i]==v) { return i; } } std::cerr<<"Cannot find variable "<<v.name()<<" in space "<<*this<<"\n"; assert(false); }


template<class Op, class Arg> struct UnaryExpression : Formula< UnaryExpression<Op,Arg> > {
    UnaryExpression(const Op& o, const Arg& a) : op(o), arg(a) { };
    Op op; Arg arg;
};

template<class Op, class Arg1, class Arg2> struct BinaryExpression : Formula< BinaryExpression<Op,Arg1,Arg2> > {
    BinaryExpression(const Op& o, const Arg1& a1, const Arg2& a2) : op(o), arg1(a1), arg2(a2) { };
    Op op; Arg1 arg1; Arg2 arg2;
};

Constant<int> make_constant(const int& x) { return Constant<int>(x); }
Constant<Float> make_constant(const Float& x) { return Constant<Float>(x); }
Constant<Interval> make_constant(const Interval& x) { return Constant<Interval>(x); }

template<class Op, class Arg> inline
UnaryExpression<Op,Arg> make_unary_expression(const Op& op, const Arg& arg) {
    return UnaryExpression<Op,Arg>(op,arg); }

template<class Op, class Arg1, class Arg2> inline
BinaryExpression<Op,Arg1,Arg2> make_binary_expression(const Op& op, const Arg1& arg1, const Arg2& arg2) {
    return BinaryExpression<Op,Arg1,Arg2>(op,arg1,arg2); }

template<class Arg1, class Arg2> inline BinaryExpression<Add,Arg1,Arg2> operator+(const Formula<Arg1>& arg1, const Formula<Arg2>& arg2) {
    return make_binary_expression(Add(),arg1(),arg2()); }
template<class Arg1, class Arg2> inline BinaryExpression<Sub,Arg1,Arg2> operator-(const Formula<Arg1>& arg1, const Formula<Arg2>& arg2) {
    return make_binary_expression(Sub(),arg1(),arg2()); }
template<class Arg1, class Arg2> inline BinaryExpression<Mul,Arg1,Arg2> operator*(const Formula<Arg1>& arg1, const Formula<Arg2>& arg2) {
    return make_binary_expression(Mul(),arg1(),arg2()); }
//template<class Arg1, class Arg2> inline BinaryExpression<Mul,Arg1,Arg2> operator-(const Arg1& arg1, const Arg2& arg2) {
//    return make_binary_expression(Mul(),arg1,arg2); }


template<class Res, class C> inline Res eval(const Space& spc, const Constant<C>& c) {
    unsigned int as=spc.size(); return Res::constant(as,c.val); }
template<class Res> inline Res eval(const Space& spc, Variable v) {
    unsigned int as=spc.size(); unsigned int j=spc.index(v); assert(j<as); return Res::variable(as,j); }
template<class Res, class Spc, class Op, class Arg1,class Arg2> inline
Res eval(const Spc& spc, const BinaryExpression<Op,Arg1,Arg2>& e) {
    Res a1=eval<Res>(spc,e.arg1);
    Res a2=eval<Res>(spc,e.arg2);
    return e.op(a1,a2);
    return e.op.template operator()<Res>(a1,a2);
}

/*
template<class Res, class Arg1,class Arg2> inline Res eval(int as, BinaryExpression<Add,Arg1,Arg2> e);
template<class Res, class Arg1,class Arg2> inline Res eval(int as, BinaryExpression<Sub,Arg1,Arg2> e);
template<class Res, class Arg1,class Arg2> inline Res eval(int as, BinaryExpression<Mul,Arg1,Arg2> e);
template<class Res, class Arg1,class Arg2> inline Res eval(int as, BinaryExpression<Div,Arg1,Arg2> e);
*/

template<class Arg1> inline BinaryExpression<Add,Arg1,Constant<Float> > operator+(const Formula<Arg1>& arg1, const Float& arg2) {
    return make_binary_expression(Add(),arg1(),make_constant(arg2)); }
template<class Arg1> inline BinaryExpression<Sub,Arg1,Constant<Float> > operator-(const Formula<Arg1>& arg1, const Float& arg2) {
    return make_binary_expression(Sub(),arg1(),make_constant(arg2)); }
template<class Arg1> inline BinaryExpression<Mul,Arg1,Constant<Float> > operator*(const Formula<Arg1>& arg1, const Float& arg2) {
    return make_binary_expression(Mul(),arg1(),make_constant(arg2)); }
template<class Arg1> inline BinaryExpression<Div,Arg1,Constant<Float> > operator/(const Formula<Arg1>& arg1, const Float& arg2) {
    return make_binary_expression(Div(),arg1(),make_constant(arg2)); }

template<class Arg2> inline BinaryExpression<Add,Constant<Float>,Arg2 > operator+(const Float& arg1, const Formula<Arg2>& arg2) {
    return make_binary_expression(Add(),make_constant(arg1),arg2()); }
template<class Arg2> inline BinaryExpression<Sub,Constant<Float>,Arg2 > operator-(const Float& arg1, const Formula<Arg2>& arg2) {
    return make_binary_expression(Sub(),make_constant(arg1),arg2()); }
template<class Arg2> inline BinaryExpression<Mul,Constant<Float>,Arg2 > operator*(const Float& arg1, const Formula<Arg2>& arg2) {
    return make_binary_expression(Mul(),make_constant(arg1),arg2()); }
template<class Arg2> inline BinaryExpression<Div,Constant<Float>,Arg2 > operator/(const Float& arg1, const Formula<Arg2>& arg2) {
    return make_binary_expression(Div(),make_constant(arg1),arg2()); }

//template<class Arg1, class Arg2> inline BinaryExpression<Add,Arg1,Constant<Arg2> > operator+(const Formula<Arg1>& arg1, const Arg2& arg2) {
//    return make_binary_expression(Add(),arg1(),make_constant(arg2)); }
//template<class Arg1, class Arg2> inline BinaryExpression<Add,Arg1,Constant<Arg2> > operator*(const Formula<Arg1>& arg1, const Arg2& arg2) {
//    return make_binary_expression(Add(),arg1(),make_constant(arg2)); }

template<class V, class F>
struct Assignment {
    V lhs; F rhs;
};


template<class F1, class Op, class F2>
struct Comparison {
    F1 lhs; Op op; F2 rhs;
};

template<class F> Assignment<Variable,F> Variable::operator=(const Formula<F>& f) const { Assignment<Variable,F> a={*this,f()}; return a; }
template<class F> Assignment<NextVariable,F> NextVariable::operator=(const Formula<F>& f) const { Assignment<NextVariable,F> a={*this,f()}; return a; }
template<class F> Assignment<DottedVariable,F> DottedVariable::operator=(const Formula<F>& f) const { Assignment<DottedVariable,F> a={*this,f()}; return a; }

template<class F1, class F2> Comparison<F1,Gtr,F2> operator>(const Formula<F1>& f1, const Formula<F2>& f2) {
    Comparison<F1,Gtr,F2> c={f1(),Gtr(),f2()}; return c; }
template<class F1, class F2> Comparison<F1,Less,F2> operator<(const Formula<F1>& f1, const Formula<F2>& f2) {
    Comparison<F1,Less,F2> c={f1(),Less(),f2()}; return c; }
//Comparison operator>(const Formula& f1, const double& x2) {  Comparison c={f1,GREATER,0.0}; return c; }

inline std::ostream& operator<<(std::ostream& os, const Add&) { return os<<"+"; }
inline std::ostream& operator<<(std::ostream& os, const Sub&) { return os<<"-"; }
inline std::ostream& operator<<(std::ostream& os, const Mul&) { return os<<"*"; }
inline std::ostream& operator<<(std::ostream& os, const Div&) { return os<<"/"; }
inline std::ostream& operator<<(std::ostream& os, const Gtr&) { return os<<">"; }
inline std::ostream& operator<<(std::ostream& os, const Less&) { return os<<"<"; }

inline std::ostream& operator<<(std::ostream& os, const Variable& v) { return os<<v.name(); }
inline std::ostream& operator<<(std::ostream& os, const DottedVariable& v) { return os<<v.name(); }
inline std::ostream& operator<<(std::ostream& os, const NextVariable& v) { return os<<v.name(); }

template<class C> std::ostream& operator<<(std::ostream& os, const Constant<C>& f) {
    return os << f.val; }

template<class Op, class Arg> std::ostream& operator<<(std::ostream& os, const UnaryExpression<Op,Arg>& f) {
    return os << f.op << "(" << f.arg << ")"; }

template<class Op, class Arg1, class Arg2> std::ostream& operator<<(std::ostream& os, const BinaryExpression<Op,Arg1,Arg2>& f) {
    return os << "(" << f.arg1 << f.op << f.arg2 << ")"; }

template<class V, class F> std::ostream& operator<<(std::ostream& os, const Assignment<V,F>& a) {
    return os << a.lhs << ":=" << a.rhs; }

template<class F1, class Op, class F2> std::ostream& operator<<(std::ostream& os, const Comparison<F1,Op,F2>& c) {
    return os << c.lhs << c.op << c.rhs; }

class Reset
{
  public:
    Reset(const Space& domain, const Space& codomain) : _domain(domain), _codomain(codomain),
        _function(result_size(),Polynomial<Interval>(argument_size())) { }
    int result_size() const { return _codomain.size(); }
    int argument_size() const { return _codomain.size(); }
    template<class F> Reset& operator,(Assignment<NextVariable,F> a) {
        assert(_codomain.index(a.lhs.var)<result_size());
        _function[_codomain.index(a.lhs.var)]=eval<Polynomial<Interval> >(_domain,a.rhs); return *this; }
    Space _domain; Space _codomain;
    std::vector<Polynomial<Interval> > _function;
};

class Dynamic 
{
  public:
    Dynamic(const Space& spc) : _domain(spc), _function(spc.size(),Polynomial<Interval>(spc.size())) { }
    int size() const { return _function.size(); }
    const Polynomial<Interval>& operator[](unsigned int i) const { return _function[i]; }
    Polynomial<Interval>& operator[](unsigned int i) { return _function[i]; }
    //int result_size() const { return _v.size(); }
    //int argument_size() const { return _v[0].argument_size(); }
    template<class F> Dynamic& operator,(Assignment<DottedVariable,F> a) {
        assert(_domain.index(a.lhs.var)<size());
        _function[_domain.index(a.lhs.var)]=eval<Polynomial<Interval> >(_domain,a.rhs); return *this; }
    Space _domain;
    std::vector<Polynomial<Interval> > _function;
};

class Guard
{
  public:
    Guard(const Space& domain) : _domain(domain), _function(Polynomial<Interval>(argument_size())) { }
    int argument_size() const { return _domain.size(); }
    template<class F1,class F2> Guard& operator,(Comparison<F1,Gtr,F2> c) {
        _function=eval<Polynomial<Interval> >(_domain,c.lhs-c.rhs); return *this; }
    template<class F1,class F2> Guard& operator,(Comparison<F1,Less,F2> c) {
        _function=eval<Polynomial<Interval> >(_domain,c.rhs-c.lhs); return *this; }
    Space _domain;
    Polynomial<Interval> _function;
};


std::ostream& operator<<(std::ostream& os, const Reset& r) { os<<"Reset("<<r.result_size()<<","<<r.argument_size()<<")";
    for(unsigned int i=0; i!=r._function.size(); ++i) { os<<(i==0?"[":",")<<r._function[i]; } return os<<"]"; }
std::ostream& operator<<(std::ostream& os, const Dynamic& d) { os<<"Dynamic("<<d.size()<<")";
    for(unsigned int i=0; i!=d.size(); ++i) { os<<(i==0?"[":",")<<d[i]; } return os<<"]"; }
std::ostream& operator<<(std::ostream& os, const Guard& g) { os<<"Guard("<<g.argument_size()<<")";
    os<<"["<<g._function; return os<<"]"; }


template<class X>
std::ostream& operator<<(std::ostream& os, const std::pair<Polynomial<X>,Space>& arg) {
    const Polynomial<X>& p=arg.first;
    const Space& spc=arg.second;
    assert(p.argument_size()==spc.size());
    for(typename Polynomial<X>::const_iterator piter=p.begin(); piter!=p.end(); ++piter) {
        MultiIndex a=piter->key();
        X v=piter->data();
        if(piter==p.begin()) { if(v<0) { os<<"-"; } }
        else { if(v>0) { os << " + "; } else { os << " - "; } }
        if(a.degree()==0 && !(v==1 || v==-1)) { os<<midpoint(v); }
        for(unsigned int j=0; j!=a.size(); ++j) {
            if(a[j]>=1) { os<<spc[j].name(); }
            if(a[j]>=2) { os<<"^"<<int(a[j]); }
        }
    }
    return os;
}


} // namespace Ariadne

#endif // ARIADNE_FORMULA_H
