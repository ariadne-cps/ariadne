/***************************************************************************
 *            function.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
/*! \file system/function.h
 *  \brief General function interface.
 */
 
#ifndef ARIADNE_SYSTEM_FUNCTION_H
#define ARIADNE_SYSTEM_FUNCTION_H

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

#include "../base/types.h"
#include "../numeric/numerical_traits.h"
#include "../linear_algebra/declarations.h"

namespace Ariadne {
  namespace System {

    class SyntaxError : public std::runtime_error { public: SyntaxError(const std::string& what) : std::runtime_error(what) { } };

    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief Abstract base class for (differentiable) functionss.
     * 
     * The function is specified by the method operator()(const LinearAlgebra::Vector< Interval<R> >& A) const,
     * This method should compute an interval vector \f$\overline{f}(A)\f$ with the
     * following properties:
     *   -# \f$f(A)\subset\overline{f}(A)\f$,
     *   -# If \f$A_1\subset A_0\f$, then \f$\overline{f}(A_1)\subset 
     *       \overline{f}(A_0)\f$, and
     *   -# If \f$\bigcap_{n\in\mathbb{N}}A_n=\{x\}\f$, then 
     *       \f$\bigcap_{n\in\mathbb{N}}\overline{f}(A_n)=\{f(x)\}\f$.
     *
     * More succinctly, we say that \f$\overline{f}(A_n)\f$ converges monotonically 
     * as \f$A_n\f$ tends to a point.
     *
     * Additional accuracy can be obtained be using derivatives.
     * The method derivative(const LinearAlgebra::Vector< Interval<R> >& A) const computes the \a i th component of the derivative over the set \a A 
     * with respect to the variables in the multi-index \a j.
     */
    template<class R>
    class Function {
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      
      /*! \brief Virtual destructor. */
      virtual ~Function() { }
     
      /*! \brief Make a copy (clone) of the vector field. */
      virtual Function<R>* clone() const = 0;
     
      /*! \brief A bound for the function over a set of vectors. */
      LinearAlgebra::Vector<I> operator() (const LinearAlgebra::Vector<I>& v) const { return this->image(v); }

      /*! \brief A bound for the vector field over aa set of vectors. */
      virtual LinearAlgebra::Vector<I> image(const LinearAlgebra::Vector<I>& v) const = 0;

      /*! \brief A bound for the vector field over a set of vectors. */
      virtual I derivative(const LinearAlgebra::Vector<I>& v, const size_type& i, const LinearAlgebra::MultiIndex& j) const = 0;

    
      /*! \brief The degree of differentiability of the function. */
      virtual size_type smoothness() const = 0;
      /*! \brief The dimension of the function argument. */
      virtual dimension_type argument_dimension() const = 0;
      /*! \brief The dimension of the function result. */
      virtual dimension_type result_dimension() const = 0;
    };
  
    enum Operation { PUSH, PULL, CONST, NEG, ADD, SUB, MUL, DIV };
    enum TokenType { STR, NUM, INT, REAL, END, ASSIGN, LP, RP, LSP, RSP, PLUS, MINUS, TIMES, DIVIDES };

    union ByteCode { Operation op; int var; int ind; int val; };
    struct Variable { std::string name; uint num; uint ind; };

    struct Instruction { Operation op; Variable var; };

    struct Token { 
      TokenType type; 
      std::string str;
      int val; 
      double rval; 
      Token() { }
      Token(std::string s) : type(STR), str(s) { }
      Token(int n) : type(INT), val(n) { }
      Token(double x) : type(REAL), rval(x) { }
      Token(TokenType t) : type(t) { }
    };
  



    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief Concrete class for functionss.
     */
    template<class R>
    class GeneralFunction {
      typedef typename Numeric::traits<R>::arithmetic_type A; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      
      /*! \brief Virtual destructor. */
      virtual ~GeneralFunction() { }
     
      /*! \brief Construct from a string literal. */
      GeneralFunction(const std::string& str);
     
      /*! \brief Make a copy (clone) of the vector field. */
      virtual GeneralFunction<R>* clone() const;
     
      /*! \brief A bound for the vector field over aa set of vectors. */
      virtual LinearAlgebra::Vector<A> image(const LinearAlgebra::Vector<A>& x) const;
 
      /*! \brief A bound for the vector field over a set of vectors. */
      virtual A derivative(const LinearAlgebra::Vector<A>& x, const size_type& i, const LinearAlgebra::MultiIndex& j) const;

      /*! \brief A bound for the vector field over a set of vectors. */
      virtual LinearAlgebra::Matrix<A> jacobian(const LinearAlgebra::Vector<A>& x) const;

      /*! \brief The degree of differentiability of the function. */
      virtual size_type smoothness() const;
      /*! \brief The dimension of the function argument. */
      virtual dimension_type argument_dimension() const;
      /*! \brief The dimension of the function result. */
      virtual dimension_type result_dimension() const;

      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
   
      void _evaluate(A** args) const;

      uint variable_number(const std::string& name) const;
     private:

      Token get_token(std::istream& is);
      Token peek_token(std::istream& is);
      Variable get_variable(std::istream& is);

      std::string read_string(std::istream& is);
      int read_int(std::istream& is);
      void read_statement(std::istream& is);
      void read_assignment(std::istream& is);
      void read_expression(std::istream& is);
      void read_factor(std::istream& is);
      void read_term(std::istream& is);
      void read_primary(std::istream& is);
      std::vector<R> _parameters;
      std::vector<R> _constants;
      std::vector<ByteCode> _operations;
      std::vector<std::string> _variable_names;
    };
   
    template<class R> inline
    std::istream& operator>>(std::istream& is, GeneralFunction<R>& f) {
      return f.read(is);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const  GeneralFunction<R>& f) {
      return f.write(os);
    };
  }
}

#endif /* ARIADNE_SYSTEM_FUNCTION_H */
