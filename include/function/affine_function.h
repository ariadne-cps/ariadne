/***************************************************************************
 *            affine_function.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file function/affine_function.h
 */

#ifndef ARIADNE_AFFINE_FUNCTION_H
#define ARIADNE_AFFINE_FUNCTION_H

#include "base/types.h"

#include "numeric/declarations.h"
#include "numeric/traits.h"

#include "linear_algebra/declarations.h"
#include "linear_algebra/vector.h"

#include "function/function_interface.h"


namespace Ariadne {


    /*! \brief An affine function \f$f(x)=Ax+b\f$ on Euclidean space. 
     *  \ingroup FunctionTypes
     */
    template<class R>
    class AffineFunction
      : public FunctionInterface<R> 
    {
      typedef typename traits<R>::approximate_arithmetic_type AA;
      typedef typename traits<R>::arithmetic_type F;
      typedef typename traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real number used to describe the system. */
      typedef R real_type;
      
      /*! \brief Default constructor constructs a function on a zero-dimensional space. */
      explicit AffineFunction() {}
      /*! \brief Construct from the matrix \f$A\f$ and the vector \f$b\f$. */
      explicit AffineFunction(const Matrix<F>& A, const Vector<F>& b)
        : _a(A), _b(b) { }
      /*! \brief Construct from the vector \f$b\f$ and the matrix \f$A\f$. */
      explicit AffineFunction(const Vector<F>& b, const Matrix<F>& A)
        : _a(A), _b(b) { }
      /*! \brief Construct a linear function from the matrix \f$A\f$. */
      explicit AffineFunction(const Matrix<F>& A)
        : _a(A), _b(A.number_of_rows()) { }
      /*! \brief Construct a translation from the vector \f$b\f$. */
      explicit AffineFunction(const Vector<F>& b)
        : _a(Matrix<R>::identity(b.size())), _b(b) { }
      
      /*! \brief Copy constructor. */
      AffineFunction(const AffineFunction<R>& f)
        : _a(f._a), _b(f._b) { }
      /*! \brief Assignment operator. */
      AffineFunction<R>& operator=(const AffineFunction<R>& f) {
        this->_a=f._a; this->_b=f._b; return *this; }
      /*! \brief Returns a pointer to a dynamically-allocated copy of the function. */
      AffineFunction<R>* clone() const { return new AffineFunction<R>(*this); }

      
      /*! \brief  An approximation to the image of an approximate point. */
      Vector<F> evaluate(const Vector<F>& x) const;
      /*! \brief The Jacobian derivative matrix at a point. */
      Matrix<F> jacobian(const Vector<F>& x) const;
      /*! \brief All the derivative values up to degree \a s. */
      DifferentialVector<F> derivative(const Vector<F>& x, const smoothness_type& s) const;
      /*! \brief All the derivative values up to degree \a s, computed using approximate arithmetic. */
      SparseDifferentialVector<AA> expansion(const Vector<AA>& x, const smoothness_type& s) const;


           
      /*! \brief  The linear transformation of the function. */
      const Matrix<F>& A() const { return _a; }
      /*! \brief  The offset vector of the function. */
      const Vector<F>& b() const { return _b; }
      
      /*! \brief The size of the result. */
      virtual size_type result_size() const {
        return _b.size(); }

      /*! \brief  The size of the argument. */
      virtual size_type argument_size() const {
        return _a.number_of_columns(); }
      
      /*! \brief The smoothness of the result. */
      virtual smoothness_type smoothness() const { 
        return std::numeric_limits<smoothness_type>::max(); }
      

      /*! \brief  The name of the system. */
      std::string name() const { return "AffineFunction"; }
      
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     protected:
      Matrix<F> _a;
      Vector<F> _b;
    };


  
} // namespace Ariadne


#endif /* ARIADNE_AFFINE_FUNCTION_H */
