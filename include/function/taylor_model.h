/***************************************************************************
 *            taylor_model.h
 *
 *  Copyright  2007 Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file taylor_model.h
 *  \brief TaylorModels.
 */

#ifndef ARIADNE_TAYLOR_MODEL_H
#define ARIADNE_TAYLOR_MODEL_H

#include <iosfwd>
#include <string>
#include <sstream>

#include "base/types.h"
#include "base/array.h"
#include "linear_algebra/declarations.h"


namespace Ariadne {
  
  namespace Geometry { template<class R> class Rectangle; }
  namespace Output { class latexstream; }

  namespace Function {
  
    class MultiIndex;
    template<class R> class TaylorModel;

    template<class R0, class R1, class R2> void add(TaylorModel<R0>&, const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R0, class R1, class R2> void sub(TaylorModel<R0>&, const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R0, class R1, class R2> void mul(TaylorModel<R0>&, const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R0, class R1> void pow(TaylorModel<R0>&, const TaylorModel<R1>&, const unsigned int&);

    template<class R0, class R1> void scale(TaylorModel<R0>&, const R1&);

    template<class R0, class R1> void derivative(TaylorModel<R0>&, const TaylorModel<R1>&, const size_type&);
    template<class R0, class R1, class R2> void compose(TaylorModel<R0>&, const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R0, class R1, class R2> void inverse(TaylorModel<R0>&, const TaylorModel<R1>&, const LinearAlgebra::Vector<R2>&);
    template<class R0, class R1, class R2> void implicit(TaylorModel<R0>&, const TaylorModel<R1>&, const LinearAlgebra::Vector<R2>&);
    template<class R0, class R1, class R2> void combine(TaylorModel<R0>&, const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R0, class R1, class R2> void join(TaylorModel<R0>&, const TaylorModel<R1>&, const TaylorModel<R2>&);

    template<class R0, class R1, class R2> void truncate(TaylorModel<R0>&, const TaylorModel<R1>&, const Geometry::Rectangle<R1>&, const size_type&, const size_type&);


    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator+(const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator-(const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const TaylorModel<R1>&, const TaylorModel<R2>&);

    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const R1&, const TaylorModel<R2>&);
    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const TaylorModel<R1>&, const R2&);
    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator/(const TaylorModel<R1>&, const R2&);

    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> compose(const TaylorModel<R1>&, const TaylorModel<R2>&);
    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1>::arithmetic_type> inverse(const TaylorModel<R1>&, const LinearAlgebra::Vector<R2>&);
    template<class R1, class R2> TaylorModel<typename Numeric::traits<R1>::arithmetic_type> implicit(const TaylorModel<R1>&, const LinearAlgebra::Vector<R2>&);
    template<class R> TaylorModel<typename Numeric::traits<R>::arithmetic_type> pow(const TaylorModel<R>& p, const unsigned int& n);
    template<class R> TaylorModel<typename Numeric::traits<R>::arithmetic_type> derivative(const TaylorModel<R>&, const size_type& k);
    template<class R> TaylorModel<typename Numeric::traits<R>::interval_type> truncate(const TaylorModel<R>&, const Geometry::Rectangle<typename Numeric::traits<R>::number_type>&, const size_type&, const size_type&);

    template<class R> std::ostream& operator<<(std::ostream&, const TaylorModel<R>&);
  

    /*! \brief A taylor_model with multivalued output, using a den.
     *  \ingroup FunctionModel
     */
    template<class R>
    class TaylorModel {
      typedef typename Numeric::traits<R>::number_type Flt;
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
     public:
      /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
      TaylorModel();
      /*! \brief Construct from a string literal. */
      TaylorModel(const std::string& s);
      /*! \brief The zero Taylor model in \a as variables with size \a rs image, order \a d and smoothness \a d. */
      TaylorModel(const size_type& rs, const size_type& as, const size_type& d);
      /*! \brief The zero Taylor model in \a as variables with size \a rs image, order \a d and smoothness \a s. */
      TaylorModel(const size_type& rs, const size_type& as, const size_type& d, const size_type& s);
      /*! \brief Construct a Taylor model in \a as variables with size \a rs image, order \a d and smoothness \a s from the data array \a a. */
      template<class RR> TaylorModel(const size_type& rs, const size_type& as, const size_type& d, const size_type& s, const RR* a);
      /*! \brief Copy constructor. */
      TaylorModel(const TaylorModel<R>& p);
      /*! \brief Copy assignment. */
      TaylorModel<R>& operator=(const TaylorModel<R>& p);
      /*! \brief Conversion constructor. */
      template<class RR> explicit TaylorModel(const TaylorModel<RR>& p);
      /*! \brief Conversion assignment. */
      template<class RR> TaylorModel<R>& operator=(const TaylorModel<RR>& p);
        
      /*! \brief Equality operator. */
      bool operator==(const TaylorModel<R>& p) const;
      /*! \brief Inequality operator. */
      bool operator!=(const TaylorModel<R>& p) const;

      /*! \brief The size of the argument. */
      size_type argument_size() const;
      /*! \brief The size of the result. */
      size_type result_size() const;
      /*! \brief The order of the Taylor model. */
      size_type order() const;
      /*! \brief The smoothness of the function. */
      size_type smoothness() const;
      /*! \brief The data used to define the Taylor model. */
      const array<R>& data() const;
      
      /*! \brief Set the \a j th value of the \a i th component to \a x. */
      void set(const size_type& i, const MultiIndex& j, const R& x); 
      /*! \brief A reference to \a j th value of the i th component. */
      R& at(const size_type& i, const MultiIndex& j);
      /*! \brief The \a j th value of the i th component. */
      const R& get(const size_type& i, const MultiIndex& j) const;
      /*! \brief Resize to a Taylor model in \a as variables with size \a rs image, order \a d and smoothness \a s. */
      void resize(const size_type& rs, const size_type& as, const size_type& d, const size_type& s);

      /*! \brief The \a i th component Taylor model. */
      TaylorModel<R> component(const size_type& i) const;

      /*! \brief Truncate the Taylor model to a Taylor model of order \a d and smoothness \a s within the domain \a domain. */
      TaylorModel<I> truncate(const size_type& order, const size_type& smoothness, const Geometry::Rectangle<R>& domain) const;

      /*! \brief Evaluate the Taylor model at the point \a x. */
      LinearAlgebra::Vector<F> evaluate(const LinearAlgebra::Vector<F>& x) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      LinearAlgebra::Matrix<F> jacobian(const LinearAlgebra::Vector<F>& s) const;

      /*! \brief The derivative of the model with respect to the \a k<sup>th</sup> variable. */
      TaylorModel<F> derivative(const size_type& k) const;

      /*! \brief The domain of validity of the Taylor model. */
      Geometry::Rectangle<Flt> domain() const;

      /*! \brief The zero Taylor model with result size \a rs and argument size \a as. */
      static TaylorModel<R> zero(const size_type& rs, const size_type& as);
      /*! \brief The unit Taylor model with result size 1 and argument size \a as. */
      static TaylorModel<R> one(const size_type& as);
      /*! \brief The constant Taylor model with result size 1 and argument size \a as. */
      static TaylorModel<R> constant(const size_type& as, const R& c);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);

#ifdef DOXYGEN
      /*! \brief Addition. */
      friend template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator+(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Subtraction. */
      friend template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator-(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Multiplication. At least one argument must be scalar-valued. */
      friend template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const TaylorModel<R1>&, const TaylorModel<R2>&);

      /*! \brief Multiplication by a scalar. */
      friend template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const R1&, const TaylorModel<R2>&);
      /*! \brief Multiplication by a scalar. */
      friend template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const TaylorModel<R1>&, const R2&);
      /*! \brief Division by a scalar. */
      friend template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> operator/(const TaylorModel<R1>&, const R2&);

      /*! \brief Composition \f$p\circ q(x)=p(q(x))\f$. */
      friend template<class R1, class R2> TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> compose(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Power of a scalar Taylor model. */
      friend template<class R> TaylorModel<typename Numeric::traits<R>::arithmetic_type> pow(const TaylorModel<R>& p, const unsigned int& n);
      /*! \brief Derivative with respect to variable \a k. */
      friend template<class R> TaylorModel<typename Numeric::traits<R>::arithmetic_type> derivative(const TaylorModel<R>&, const size_type& k);
      /*! \brief Truncate within \a r to a Taylor model of order at most \a d, putting the error into terms of order \a s. */
      friend template<class R> TaylorModel<typename Numeric::traits<R>::arithmetic_type> truncate(const TaylorModel<R>& p, const Rectangle<R>& bb, const size_type& d, const size_type& s);
#endif
     private:
      static void instantiate();
      template<class RR> array< array<F> > _powers(const LinearAlgebra::Vector<RR>&) const;
      void _compute_jacobian() const;
      void _set_argument_size(const size_type& n);
      size_type _compute_maximum_component_size() const;
     private:
      template<class R0,class R1,class R2> friend void add(TaylorModel<R0>&,const TaylorModel<R1>&,const TaylorModel<R2>&);
      template<class R0,class R1,class R2> friend void sub(TaylorModel<R0>&,const TaylorModel<R1>&,const TaylorModel<R2>&);
      template<class R0,class R1,class R2> friend void mul(TaylorModel<R0>&,const TaylorModel<R1>&,const TaylorModel<R2>&);
      template<class R0,class R1,class R2> friend void div(TaylorModel<R0>&,const TaylorModel<R1>&,const TaylorModel<R2>&);
      template<class R0,class R1,class R2> friend void compose(TaylorModel<R0>&,const TaylorModel<R1>&,const TaylorModel<R2>&);
      template<class R0,class R1> friend void scale(TaylorModel<R0>&,const R1&);
     private:
      /* Components of the map. */
      size_type _result_size; 
      size_type _argument_size;
      size_type _order; 
      size_type _smoothness; 
      array<R> _data;
    };
    
  }
}

#include "taylor_model.inline.h"

#endif /* ARIADNE_TAYLOR_MODEL_H */
