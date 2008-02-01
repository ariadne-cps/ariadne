/***************************************************************************
 *            numeric/operators.h
 *
 *  Copyright  2007 Pieter Collins
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
 
/*! \file numeric/operators.h
 *  \brief Numerical operators templates for deferred evaluation.
 */

#ifndef ARIADNE_NUMERIC_OPERATORS_H
#define ARIADNE_NUMERIC_OPERATORS_H

#include <iostream>

#include "numeric/rounding.h"
#include "numeric/expression.h"

#define ARIADNE_NULLARY_FUNCTION(Fn,fn,ifn,str)                         \
  class Fn { public: template<class R> void operator()(R& r) const { ifn(r); } }; \
  inline std::ostream& operator<<(std::ostream& os, const Fn&) { \
    return os << str; }                                          \
  inline Expression< Nullary<Fn> > fn() { \
    return make_expression(Fn()); }       \
  template<class R> inline R fn() {      \
    R r; ifn(r); return r; } \

#define ARIADNE_UNARY_FUNCTION(Fn,fn,ifn,str)                                   \
  class Fn { public: template<class R, class A> void operator()(R& r, const A& a) const { ifn(r,a); } }; \
  inline std::ostream& operator<<(std::ostream& os, const Fn&) { return os << str; } \
  template<class R, class A> inline R fn(const A& a) { \
    R r; ifn(r,a); return r; } \
  template<class A> inline Expression< Unary<Fn,A> > fn(const Scalar<A>& a) { \
    return make_expression(Fn(),static_cast<const A&>(a)); } \

#define ARIADNE_BINARY_FUNCTION(Fn,fn,ifn,str)                          \
  class Fn { public: template<class R, class A1, class A2> void operator()(R& r, const A1& a1, const A2& a2) const { ifn(r,a1,a2); } }; \
  inline std::ostream& operator<<(std::ostream& os, const Fn&) { return os << str; } \
  template<class R, class A1, class A2> inline R fn(const A1& a1, const A2& a2) {    \
    R r; ifn(r,a1,a2); return r; }                                        \
  template<class A1,class A2> inline \
  Expression< Binary<Fn,A1,A2> > fn(const Scalar<A1>& a1,const Scalar<A2>& a2) { \
    return make_expression(Fn(),static_cast<const A1&>(a1),static_cast<const A2&>(a2)); } \
  template<class A> inline \
  Expression< Binary<Fn,A,int> > fn(const Scalar<A>& a1,const int& a2) { \
	return make_expression(Fn(),static_cast<const A&>(a1),a2); } \
  template<class A> inline \
  Expression< Binary<Fn,int,A> > fn(const int& a1,const Scalar<A>& a2) {  \
	return make_expression(Fn(),a1,static_cast<const A&>(a2)); } \
  template<class A> inline \
  Expression< Binary<Fn,A,uint> > fn(const Scalar<A>& a1,const uint& a2) { \
	return make_expression(Fn(),static_cast<const A&>(a1),a2); } \
  template<class A> inline \
  Expression< Binary<Fn,uint,A> > fn(const uint& a1,const Scalar<A>& a2) {  \
	return make_expression(Fn(),a1,static_cast<const A&>(a2)); } \
  template<class A> inline \
  Expression< Binary<Fn,A,uint> > fn(const Scalar<A>& a1,const double& a2) { \
	return make_expression(Fn(),static_cast<const A&>(a1),a2); } \
  template<class A> inline \
  Expression< Binary<Fn,double,A> > fn(const uint& a1,const Scalar<A>& a2) {  \
	return make_expression(Fn(),a1,static_cast<const A&>(a2)); } \

#define ARIADNE_INPLACE_OPERATOR(op,ifn) \
  template<class X, class Y> inline X& op(Value<X>& x, const Expression<Y>& y) { \
	ifn(x.promote(),x.promote(),X(y)); return x.promote(); } \
  template<class X, class Y> inline X& op(Value<X>& x, const Y& y) { \
	ifn(x.promote(),x.promote(),y); return x.promote(); } \
  template<class X, class Y> inline X& op(X& x, const Expression<Y>& y) { \
	ifn(x,x,X(y)); return x; } \

#define ARIADNE_UNARY_OPERATOR(op,Fn) \
  template<class X> inline Expression< Unary<Fn,X> > op(const Scalar<X>& x) { \
	return make_expression(Fn(),x.promote()); } \

#define ARIADNE_BINARY_OPERATOR(op,Fn) \
  template<class X,class Y> inline Expression< Binary<Fn,X,Y> > op(const Scalar<X>& x,const Scalar<Y>& y) {  \
	return make_expression(Fn(),x.promote(),y.promote()); } \
  template<class X> inline Expression< Binary<Fn,X,unsigned short int> > op(const Scalar<X>& x,const unsigned short int& y) {  \
	return make_expression(Fn(),x.promote(),y); } \
  template<class X> inline Expression< Binary<Fn,X,unsigned int> > op(const Scalar<X>& x,const unsigned int& y) {  \
	return make_expression(Fn(),x.promote(),y); } \
  template<class X> inline Expression< Binary<Fn,X,int> > op(const Scalar<X>& x,const int& y) {  \
	return make_expression(Fn(),x.promote(),y); } \
  template<class X> inline Expression< Binary<Fn,X,double> > op(const Scalar<X>& x,const double& y) {  \
	return make_expression(Fn(),x.promote(),y); } \
  template<class X> inline Expression< Binary<Fn,unsigned short int,X> > op(const unsigned short int& x,const Scalar<X>& y) {  \
	return make_expression(Fn(),x,y.promote()); } \
  template<class X> inline Expression< Binary<Fn,unsigned int,X> > op(const unsigned int& x,const Scalar<X>& y) {  \
	return make_expression(Fn(),x,y.promote()); } \
  template<class X> inline Expression< Binary<Fn,int,X> > op(const int& x,const Scalar<X>& y) {  \
	return make_expression(Fn(),x,y.promote()); } \
  template<class X> inline Expression< Binary<Fn,double,X> > op(const double& x,const Scalar<X>& y) {  \
	return make_expression(Fn(),x,y.promote()); } \


#define ARIADNE_COMPARISON(R,X,Y) \
  inline R operator==(const X& x, const Y& y) { return cmp(x,y)==0; } \
  inline R operator!=(const X& x, const Y& y) { return cmp(x,y)!=0; } \
  inline R operator<=(const X& x, const Y& y) { return cmp(x,y)<=0; } \
  inline R operator>=(const X& x, const Y& y) { return cmp(x,y)>=0; } \
  inline R operator< (const X& x, const Y& y) { return cmp(x,y)< 0; } \
  inline R operator> (const X& x, const Y& y) { return cmp(x,y)> 0; } \

#define ARIADNE_MIXED_FUNCTION_COMPARISON(R,X,Y) \
  inline R operator==(const X& x, const Y& y) { return cmp(x,y)==0; } \
  inline R operator!=(const X& x, const Y& y) { return cmp(x,y)!=0; } \
  inline R operator<=(const X& x, const Y& y) { return cmp(x,y)<=0; } \
  inline R operator>=(const X& x, const Y& y) { return cmp(x,y)>=0; } \
  inline R operator< (const X& x, const Y& y) { return cmp(x,y)< 0; } \
  inline R operator> (const X& x, const Y& y) { return cmp(x,y)> 0; } \
  inline R operator==(const Y& y, const X& x) { return cmp(x,y)==0; } \
  inline R operator!=(const Y& y, const X& x) { return cmp(x,y)!=0; } \
  inline R operator<=(const Y& y, const X& x) { return cmp(x,y)>=0; } \
  inline R operator>=(const Y& y, const X& x) { return cmp(x,y)<=0; } \
  inline R operator< (const Y& y, const X& x) { return cmp(x,y)> 0; } \
  inline R operator> (const Y& y, const X& x) { return cmp(x,y)< 0; } \

namespace Ariadne {
  namespace Numeric {
    

    ARIADNE_NULLARY_FUNCTION(NaN,nan,nan_,"nan");
    ARIADNE_NULLARY_FUNCTION(Inf,inf,inf_,"inf");
    ARIADNE_NULLARY_FUNCTION(Eps,eps,eps_,"eps");

    ARIADNE_BINARY_FUNCTION(Max,max,max_,"max");
    ARIADNE_BINARY_FUNCTION(Min,min,min_,"min");
    ARIADNE_UNARY_FUNCTION(Pos,pos,pos_,"pos");
    ARIADNE_UNARY_FUNCTION(Neg,neg,neg_,"neg");
    ARIADNE_UNARY_FUNCTION(Abs,abs,abs_,"abs");
    ARIADNE_BINARY_FUNCTION(Add,add,add_,"add");
    ARIADNE_BINARY_FUNCTION(Sub,sub,sub_,"sub");
    ARIADNE_BINARY_FUNCTION(Mul,mul,mul_,"mul");
    ARIADNE_BINARY_FUNCTION(Div,div,div_,"div");
    ARIADNE_BINARY_FUNCTION(Pow,pow,pow_,"pow");
    ARIADNE_BINARY_FUNCTION(Mod,mod,mod_,"mod");
    ARIADNE_BINARY_FUNCTION(Quot,quot,quot_,"quot");
    ARIADNE_BINARY_FUNCTION(Rem,rem,rem_,"rem");
  
    ARIADNE_UNARY_FUNCTION(Factorial,fac,fac_,"fac");
    ARIADNE_BINARY_FUNCTION(Choose,bin,bin_,"bin");
    ARIADNE_BINARY_FUNCTION(LCM,lcm,lcm_,"lcm");
    ARIADNE_BINARY_FUNCTION(GCD,gcd,gcd_,"gcd");

    ARIADNE_UNARY_FUNCTION(Floor,floor,floor_,"floor");
    ARIADNE_UNARY_FUNCTION(Ceil,ceil,ceil_,"ceil");

    ARIADNE_UNARY_FUNCTION(Sqrt,sqrt,sqrt_,"sqrt");
    ARIADNE_BINARY_FUNCTION(Hypot,hypot,hypot_,"hypot");
    ARIADNE_UNARY_FUNCTION(Exp,exp,exp_,"exp");
    ARIADNE_UNARY_FUNCTION(Log,log,log_,"log");

    ARIADNE_NULLARY_FUNCTION(Pi,pi,pi_,"pi");
    ARIADNE_UNARY_FUNCTION(Sin,sin,sin_,"sin");
    ARIADNE_UNARY_FUNCTION(Cos,cos,cos_,"cos");
    ARIADNE_UNARY_FUNCTION(Tan,tan,tan_,"tan");
    ARIADNE_UNARY_FUNCTION(Asin,asin,asin_,"asin");
    ARIADNE_UNARY_FUNCTION(Acos,acos,acos_,"acos");
    ARIADNE_UNARY_FUNCTION(Atan,atan,atan_,"atan");

    ARIADNE_UNARY_OPERATOR(operator+,Pos);
    ARIADNE_UNARY_OPERATOR(operator-,Neg);
    ARIADNE_BINARY_OPERATOR(operator+,Add);
    ARIADNE_BINARY_OPERATOR(operator-,Sub);
    ARIADNE_BINARY_OPERATOR(operator*,Mul);
    ARIADNE_BINARY_OPERATOR(operator/,Div);
	  
    ARIADNE_INPLACE_OPERATOR(operator+=,add_);
    ARIADNE_INPLACE_OPERATOR(operator-=,sub_);
    ARIADNE_INPLACE_OPERATOR(operator*=,mul_);
    ARIADNE_INPLACE_OPERATOR(operator/=,div_);
	  
   
    // Rounded operators
    template<class Rnd, class X>
    X add(const X& x, const X& y) {
      X r; add_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X sub(const X& x, const X& y) {
      X r; sub_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X mul(const X& x, const X& y) {
      X r; mul_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X div(const X& x, const int& y) {
      X r; div_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X div(const X& x, const long int& y) {
      X r; div_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X div(const X& x, const X& y) {
      X r; div_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X med(const X& x, const X& y) {
      X r; med_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X rad(const X& x, const X& y) {
      X r; rad_(r,x,y,Rnd()); return r; }

    template<class Rnd, class X>
    X exp(const X& x) {
      X r; exp_(r,x,Rnd()); return r; }

    template<class R, class X> inline R approx(const X& x) {
      R r; set_(r,x,round_approx); return r; }




  }   
}
  

#endif /* ARIADNE_NUMERIC_OPERATORS_H */
