/***************************************************************************
 *            vector_expression.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file vector_expressions.h
 *  \brief Vector expression templates.
 */

#ifndef ARIADNE_VECTOR_EXPRESSION_H
#define ARIADNE_VECTOR_EXPRESSION_H 

#include <iosfwd>
#include <algorithm>

#include "../base/types.h"
#include "../base/array.h"
#include "../numeric/integer.h"
#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/exceptions.h"

namespace Ariadne {
  namespace LinearAlgebra {
    
    struct plus {
      template<class R1,class R2> 
      typename Numeric::traits<R1,R2>::arithmetic_type 
      operator() (const R1& x1, const R2& x2) const { return x1+x2; }
    };

    struct minus {
      template<class R1,class R2> 
      typename Numeric::traits<R1,R2>::arithmetic_type 
      operator() (const R1& x1, const R2& x2) const { return x1-x2; }
    };

    struct times {
      template<class R1,class R2> 
      typename Numeric::traits<R1,R2>::arithmetic_type 
      operator() (const R1& x1, const R2& x2) const { return x1*x2; }
    };

    struct divides {
      template<class R1,class R2> 
      typename Numeric::traits<R1,R2>::arithmetic_type 
      operator() (const R1& x1, const R2& x2) const { return x1/x2; }
    };

    
    
    /*!\brief %Base class for all vector expressions. */
    template<class E>
    class VectorExpression 
    {
     public:
      /*!\brief Convert \a *this to a reference to E. */
      E& operator() () { return static_cast<E&>(*this); }
      /*!\brief Convert \a *this to a constant reference to E. */
      const E& operator() () const { return static_cast<const E&>(*this); }
    };


    
  template<class Op, class VE1, class VE2>
    class BinaryVectorVectorExpression :
      public VectorExpression< BinaryVectorVectorExpression<Op,VE1,VE2> >
    {
     public:
      //      typedef typename Numeric::traits<typename VE1::value_type, typename VE2::value_type>::arithmetic_type value_type;
      typedef Numeric::Expression< Numeric::Binary<Op, typename VE1::value_type,typename VE1::value_type> > value_type;
      BinaryVectorVectorExpression(const Op& o, const VE1& v1, const VE2& v2)
        : _ve1(v1), _ve2(v2), _op(o) { }
      size_type size() const { return _ve1.size(); }
      value_type operator()(const size_type& i) const { return value_type(_op,_ve1(i),_ve2(i)); }
      value_type operator[](const size_type& i) const { return value_type(_op,_ve1(i),_ve2(i)); }
      //      value_type operator()(const size_type& i) const { return _op(_ve1(i),_ve2(i)); }
      //      value_type operator[](const size_type& i) const { return _op(_ve1(i),_ve2(i)); }
     private:
      const VE1& _ve1; const VE2& _ve2; Op _op;
    };
    
    
    template<class Op, class VE, class SE>
    class BinaryVectorScalarExpression
      : public VectorExpression< BinaryVectorScalarExpression<Op,VE,SE> >
    {
      typedef typename VE::value_type vector_value_type;
      //typedef typename Numeric::traits<SE>::closure_type scalar_closure_type;
     public:
      //      typedef typename Numeric::traits<vector_value_type,scalar_closure_type>::arithmetic_type value_type;
      typedef Numeric::Expression< Numeric::Binary<Op, typename VE::value_type,SE> > value_type;
      BinaryVectorScalarExpression(const Op& o, const VE& ve, const SE& se) 
        : _ve(ve), _se(se), _op(o) { }
      size_type size() const { return _ve.size(); }
      value_type operator()(const size_type& i) const { return value_type(_op,_ve(i),_se); }
      value_type operator[](const size_type& i) const { return value_type(_op,_ve(i),_se); }
     private:
      const VE& _ve; 
      //const scalar_closure_type _se; 
      const SE& _se;
      Op _op;
    };
    

  }
}


#endif /* ARIADNE_VECTOR_EXPRESSION_H */
