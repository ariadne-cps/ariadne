/***************************************************************************
 *            numeric/expression.h
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
 
/*! \file numeric/expression.h
 *  \brief Numerical expression templates for deferred evaluation.
 */

#ifndef ARIADNE_NUMERIC_EXPRESSION_H
#define ARIADNE_NUMERIC_EXPRESSION_H

#include <iostream>

#ifdef DEBUG 
#warning "Debug mode"
#define ARIADNE_DEBUG_PRINT(expr) { std::cerr << expr << std::endl; }
#else
#define ARIADNE_DEBUG_PRINT(expr) { }
#endif 

namespace Ariadne {
  
    
    template<class X> class Scalar; 
    template<class E> class Expression; 
    template<class T> class Value; 
    template<class Op> class Nullary;
    template<class Op, class Arg> class Unary;
    template<class Op, class Arg1, class Arg2> class Binary;

    template<class Op> inline Expression< Nullary<Op> > make_expression(const Op& op);
    template<class Op, class Arg> inline Expression< Unary<Op,Arg> > make_expression(const Op& op, const Arg& a);
    template<class Op, class Arg1, class Arg2> inline Expression< Binary<Op,Arg1,Arg2> > make_expression(const Op&, const Arg1& a1, const Arg2& a2);

    template<class Op> inline std::ostream& operator<<(std::ostream& os, const Expression< Nullary<Op> >& e);
    template<class Op, class Arg> inline std::ostream& operator<<(std::ostream& os, const Expression< Unary<Op,Arg> >& e);
    template<class Op, class Arg1, class Arg2> inline std::ostream& operator<<(std::ostream& os, const Expression< Binary<Op,Arg1,Arg2> >& e);


    template<class X> 
    class Scalar
    {
     public:
      X& promote() { return static_cast<X&>(*this); }
      const X& promote() const { return static_cast<const X&>(*this); }
    };


    template<class V> 
    class Value 
      : public Scalar<V>
    {
    };


    template<class Op>
    class Expression< Nullary<Op> >
      : public Scalar< Expression< Nullary<Op> > >
    { 
      Op op;
     public:
      Expression(const Op& o) : op(o) { }
      template<class R> void assign_to(R& r) const { op(r); }
      template<class R> R evaluate() const { R r; this->assign_to(r); return r; }
      template<class R> operator R () const { R r; this->assign_to(r); return r; }
    };


    template<class Op, class Arg> 
    class Expression< Unary< Op, Arg > >
      : public Scalar< Expression< Unary< Op, Arg > > >
	{
      Op op; const Arg& arg;
     public:
      Expression(const Op& o, const Arg& a) : op(o), arg(a) { }
      template<class R> void assign_to(R& r) const { op(r,arg); }
      template<class R> R evaluate() const { R r; this->assign_to(r); return r; }
      template<class R> operator R () const { R r; this->assign_to(r); return r; }
    };

    template<class Op, class E> 
    class Expression< Unary< Op, Expression<E> > >
      : public Scalar< Expression< Unary< Op, Expression<E> > > >
    {
      Op op; const Expression<E>& arg; 
     public:
      Expression(const Op& o, const Expression<E>& e) : op(o), arg(e) { }
      template<class R> void assign_to (R& r) const { op(r,arg.evaluate<R>()); }
      template<class R> R evaluate() const { R r; this->assign_to(r); return r; }
      template<class R> operator R () const { R r; this->assign_to(r); return r; }
    };


    template<class Op, class Arg1, class Arg2> 
    class Expression< Binary< Op, Arg1, Arg2 > >
      : public Scalar< Expression< Binary< Op, Arg1, Arg2 > > >
    {
      Op op; const Arg1& arg1; const Arg2& arg2; 
     public:
      Expression(const Op& o, const Arg1& a1, const Arg2& a2) : op(o), arg1(a1), arg2(a2){ }
      template<class R> void assign_to(R& r) const { op(r,arg1,arg2); }
      template<class R> R evaluate() const { R r; this->assign_to(r); return r; }
      template<class R> operator R () const { R r; this->assign_to(r); return r; }
    };

    template<class Op, class Arg1, class E2> 
    class Expression< Binary< Op, Arg1, Expression<E2> > >
      : public Scalar< Expression< Binary< Op, Arg1, Expression<E2> > > >
    {
      Op op; const Arg1& arg1; const Expression<E2>& arg2; 
     public:
      Expression(const Op& o, const Arg1& a1, const Expression<E2>& a2) : op(o), arg1(a1), arg2(a2) { }
      template<class R> void assign_to(R& r) const { op(r,arg1,arg2.evaluate<R>()); }
      template<class R> R evaluate() const { R r; this->assign_to(r); return r; }
      template<class R> operator R () const { R r; this->assign_to(r); return r; }
    };

    template<class Op, class E1, class Arg2> 
    class Expression< Binary< Op, Expression<E1>, Arg2 > >
      : public Scalar< Expression< Binary< Op, Expression<E1>, Arg2 > > >
    {
      Op op; const Expression<E1>& arg1; const Arg2& arg2; 
     public:
      Expression(const Op& o, const Expression<E1>& a1, const Arg2& a2) : op(o), arg1(a1), arg2(a2) { }
      template<class R> void assign_to(R& r) const { op(r,arg1.evaluate<R>(),arg2); }
      template<class R> R evaluate() const { R r; this->assign_to(r); return r; }
      template<class R> operator R () const { R r; this->assign_to(r); return r; }
    };

    template<class Op, class E1, class E2> 
    class Expression< Binary< Op, Expression<E1>, Expression<E2> > >
      : public Scalar< Expression< Binary< Op, Expression<E1>, Expression<E2> > > >
    {
      Op op; const Expression<E1>& arg1; const Expression<E2>& arg2;
     public:
      Expression(const Op& o, const Expression<E1>& a1, const Expression<E2>& a2) : op(o), arg1(a1), arg2(a2) { }
      template<class R> void assign_to(R& r) const { op(r,arg1.evaluate<R>(),arg2.evaluate<R>()); }
      template<class R> R evaluate() const { R r; this->assign_to(r); return r; }
      template<class R> operator R () const { R r; this->assign_to(r); return r; }
    };


    template<class Op> inline 
    Expression< Nullary<Op> > make_expression(const Op& op) { 
      return Expression< Nullary<Op> >(op); }
    
    template<class Op, class Arg> inline 
    Expression< Unary<Op,Arg> > make_expression(const Op& op, const Arg& arg) { 
      return Expression< Unary<Op,Arg> >(op,arg); }

    template<class Op, class Arg1, class Arg2> inline 
    Expression< Binary<Op,Arg1,Arg2> > make_expression(const Op& op, const Arg1& a1, const Arg2& a2) { 
      return Expression< Binary<Op,Arg1,Arg2> >(op,a1,a2); }


    template<class Op> inline 
    std::ostream& operator<<(std::ostream& os, const Expression< Nullary<Op> >& e) { 
      return os << e.op << '(' << ')';
    }

    template<class Arg, class Op> inline 
    std::ostream& operator<<(std::ostream& os, const Expression< Unary<Op,Arg> >& e) { 
      return os << e.op << '(' << e.arg << ')';
    }

    template<class Arg1, class Arg2, class Op> inline 
    std::ostream& operator<<(std::ostream& os, const Expression< Binary<Op,Arg1,Arg2> >& e) { 
      return os << e.op << '(' << e.arg1 <<',' << e.arg2 << ')';
    }



  }   

  

#endif /* ARIADNE_NUMERIC_EXPRESSION_H */
