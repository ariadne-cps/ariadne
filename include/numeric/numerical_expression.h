/***************************************************************************
 *            numerical_traits.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file numerical_expression.h
 *  \brief Numerical expression templates for deferred evaluation.
 */

#ifndef _ARIADNE_NUMERICAL_EXPRESSION_H
#define _ARIADNE_NUMERICAL_EXPRESSION_H

#include <string>

#include "approximation.h"

namespace Ariadne {
  namespace Numeric {
    
    template<class R>
    class ConstantFunction
    {
     public:
      typedef R result_type;
    };

    template<class A, class R>
    class UnaryFunction
    {
     public:
      typedef A argument_type;
      typedef R result_type;
    };

    template<class A1, class A2, class R>
    class BinaryFunction
    {
     public:
      typedef A1 first_argument_type;
      typedef A2 second_argument_type;
      typedef R result_type;
    };

    template<class R, class E> 
    class Expression 
    { 
     public:
      typedef R result_type; 
      void operator()(R& r) const { static_cast<const E*>(this)->operator()(r); }
      operator R () const { R r; this->operator()(r); return r; }
      E& promote() { return static_cast<E&>(*this); }
      const E& promote() const { return static_cast<const E&>(*this); }
    };
    
    template<typename Op> 
    class ConstantExpression
      : public Expression<typename Op::result_type, ConstantExpression<Op> > 
    {
     public:
      ConstantExpression(const Op& o) : op(o) { }
      void operator() (typename Op::result_type& r) const {
        op(r);
      }
      Op op; 
    };
  
    template<typename A, typename Op> 
    class UnaryExpression
      : public Expression<typename Op::result_type, UnaryExpression<A,Op> > 
    {
     public:
      UnaryExpression(const A& a, const Op& o) : arg(a), op(o) { }
      void operator() (typename Op::result_type& r) const {
        typename A::result_type tmp;
        arg(tmp);
        op(r,tmp);
      }
      A arg; Op op; 
    };
  
    template<class Op> 
    class UnaryExpression<typename Op::argument_type,Op>
      : public Expression< typename Op::result_type, UnaryExpression<typename Op::argument_type,Op> > 
    {
      typedef typename Op::argument_type A;
     public:
      UnaryExpression(const A& a, const Op& o) : arg(a), op(o) { }
      void operator() (typename Op::result_type& r) const {
        op(r,arg);
      }
      A arg; Op op; 
    };
  

    template<typename A1, class A2, typename Op> 
    class BinaryExpression
      : public Expression<typename Op::result_type, BinaryExpression<A1,A2,Op> > 
    {
     public:
      BinaryExpression(const A1& a1, const A2& a2,const Op& o) : arg1(a1), arg2(a2), op(o) { }
      void operator()(typename Op::result_type& r) const {
        typename A1::result_type tmp1;
        arg1(tmp1);
        typename A2::result_type tmp2;
        arg2(tmp2);
        op(r,tmp1,tmp2);
      }
      A1 arg1; A2 arg2; Op op; 
    };
  
    template<typename A1, typename Op> 
    class BinaryExpression<A1,typename Op::second_argument_type,Op>
      : public Expression<typename Op::result_type, BinaryExpression<A1,typename Op::second_argument_type,Op> >
    {
      typedef typename Op::second_argument_type A2;
     public:
      BinaryExpression(const A1& a1, const A2& a2,const Op& o) : arg1(a1), arg2(a2), op(o) { }
      void operator()(typename Op::result_type& res) const {
        typename A1::result_type tmp1;
        arg1(tmp1);
        op(res,tmp1,arg2);
      }
      A1 arg1; A2 arg2; Op op; 
    };
  
    template<typename A2, typename Op> 
    class BinaryExpression<typename Op::first_argument_type,A2,Op>
      : public Expression<typename Op::result_type, BinaryExpression<typename Op::first_argument_type,A2,Op> >
    {
      typedef typename Op::first_argument_type A1;
     public:
      BinaryExpression(const A1& a1, const A2& a2,const Op& o) : arg1(a1), arg2(a2), op(o) { }
      void operator()(typename Op::result_type& res) const {
        typename A2::result_type tmp2;
        arg2(tmp2);
        op(res,arg1,tmp2);
      }
      A1 arg1; A2 arg2; Op op; 
    };
  
    template<class Op> 
    class BinaryExpression<typename Op::first_argument_type,typename Op::second_argument_type,Op>
      : public Expression< typename Op::result_type, BinaryExpression<typename Op::first_argument_type,typename Op::second_argument_type,Op> >
    {
      typedef typename Op::first_argument_type A1;
      typedef typename Op::second_argument_type A2;
     public:
      BinaryExpression(const A1& a1, const A2& a2, const Op& o) : arg1(a1), arg2(a2), op(o) { }
      void operator()(typename Op::result_type& res) const {
        op(res,arg1,arg2);
      }
      A1 arg1; A2 arg2; Op op; 
    };
  
  }   
}
  

#endif /* _ARIADNE_NUMERICAL_TRAITS_H */
