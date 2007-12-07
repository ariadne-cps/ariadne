/***************************************************************************
 *            arithmetic.h
 *
 *  Wed 18 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file arithmetic.h
 *  \brief Simple arithmetic functions using expression templates.
 */
 
#ifndef ARIADNE_ARITHMETIC_H
#define ARIADNE_ARITHMETIC_H

#include "numeric/expression.h"

namespace Ariadne {
  namespace Numeric {

#ifdef DOXYGEN

    //! \name Exact arithmetical operations. 
    //@{
    //! \ingroup Numeric
    /*! \brief Minimum. */
    template<class R> inline R min_(const R& x1, const R& x2);
    template<class R> inline Interval<R> min_(const Interval<R>& x1, 
                                                const Interval<R>& x2);

    /*! \brief Maximum. */
    template<class R> inline R max_(const R& x1, const R& x2);;
    template<class R> inline Interval<R> max_(const Interval<R>& x1, 
                                                const Interval<R>& x2);

    /*! \brief The med_ian (average) of two values. */
    template<class R> inline R med_(const R& x1, const R& x2);

    /*! \brief Absolute value. */
    template<class R> inline R abs_(const R& x);  
    template<class R> inline Interval<R> abs_(const Interval<R>& x);

    /*! \brief Unary neg_ation. */
    template<class R> inline R neg_(const R& x);
    
    /*! \brief Reciprocal. */
    template<class R> inline R rec(const R& x);
    
    /*! \brief Addition. */
    template<class R> inline R add_(const R& x1,const R& x2);
    
    /*! \brief Subtraction. */
    template<class R> inline R sub_(const R& x1,const R& x2);
    
    /*! \brief Multiplication. */
    template<class R> inline R mul_(const R& x1,const R& x2);
    
    /*! \brief Division. */
    template<class R> inline R div_(const R& x1, const R& x2);
    
    /*! \brief The pow_er of a real number type by an integer. */
    template<class R, class N> inline R pow_(const R& x, const N& n);
    
    /*! \brief The integer part of \a x, rounding towards 0. */
    template<class R> inline R floor(const R& x);

    /*! \brief The integer ceiling of \a x, rounding away from 0. */
    template<class R> inline R ceil(const R& x);
    
    //@}
#endif // DOXYGEN

/*
    template<class A1, class A2> BinaryExpression<A1,A2,Min> inline min_(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Min()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Max> inline max_(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Max()); }
    template<class A> UnaryExpression<A,Pos> inline pos_(const A& x) {
      return make_expression(x,Pos()); }
    template<class A> UnaryExpression<A,Neg> inline neg_(const A& x) {
      return make_expression(x,Neg()); }
    template<class A> UnaryExpression<A,Abs> inline abs_(const A& x) {
      return make_expression(x,Abs()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Add> inline add_(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Add()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Sub> inline sub_(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Sub()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Mul> inline mul_(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Mul()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Div> inline div_(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Div()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Pow> inline pow_(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Pow()); }

    template<class A> UnaryExpression<A,Factorial> inline factorial(const A& x) {
      return make_expression(x,Factorial()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Choose> inline factorial(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Choose()); }
    template<class A1, class A2> BinaryExpression<A1,A2,LCM> inline lcm(const A1& x1, const A2& x2) {
       return make_expression(x1,x2,LCM()); }
    template<class A1, class A2> BinaryExpression<A1,A2,GCD> inline gcd(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,GCD()); }

    template<class A> UnaryExpression<A,Floor> inline floor(const A& x) {
      return make_expression(x,Floor()); }
    template<class A> UnaryExpression<A,Floor> inline ceil(const A& x) {
      return make_expression(x,Ceil()); }

    template<class A> UnaryExpression<A,Pos> inline operator+(const A& x) {
      return make_expression(x,Pos()); }
    template<class A> UnaryExpression<A,Neg> inline operator-(const A& x) {
      return make_expression(x,Neg()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Add> inline operator+(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Add()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Sub> inline operator-(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Sub()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Mul> inline operator*(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Mul()); }
    template<class A1, class A2> BinaryExpression<A1,A2,Div> inline operator/(const A1& x1, const A2& x2) {
      return make_expression(x1,x2,Div()); }

    template<class T, class A> inline T& operator+=(T& t, const A& a) { add_(t,t,a); return t; }
    template<class T, class A> inline T& operator-=(T& t, const A& a) { sub_(t,t,a); return t; }
    template<class T, class A> inline T& operator*=(T& t, const A& a) { mul_(t,t,a); return t; }
    template<class T, class A> inline T& operator/=(T& t, const A& a) { div_(t,t,a); return t; }
*/
  }
}



#endif /* ARIADNE_ARITHMETIC_H */
