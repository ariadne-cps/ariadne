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

#include "function/taylor_derivative.h"

namespace Ariadne {
  
  namespace Geometry { template<class R> class Point; template<class R> class Box; }

  namespace Function {
  
    class MultiIndex;
    template<class R> class TaylorVariable;
    template<class R> class TaylorDerivative;
    template<class R> class TaylorModel;

 
    /*! \brief A taylor_model with multivalued output, using a den.
     *  \ingroup FunctionModel
     */
    template<class R>
    class TaylorModel {
      typedef typename Numeric::Interval<R> I;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
     public:
      /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
      TaylorModel();
      /*! \brief Construct from a domain, centre, and derivative expansion. */
      TaylorModel(const Geometry::Box<R>& domain, const LinearAlgebra::Vector<R>& centre, const TaylorDerivative<I>& td);
      /*! \brief Construct from a domain, centre, smoothness and derivative expansion. */
      TaylorModel(const Geometry::Box<R>& domain, const LinearAlgebra::Vector<R>& centre, const smoothness_type& s, const TaylorDerivative<I>& td);
      /*! \brief The zero Taylor model in \a as variables with size \a rs image, order \a d and smoothness \a d, defined on the whole space with centre at the origin. */
      TaylorModel(const size_type& rs, const size_type& as, const size_type& d, const size_type& s);
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
      const TaylorDerivative<I>& derivatives() const;

      /*! \brief The centre of the derivative expansion. */
      LinearAlgebra::Vector<R> centre() const;
     
      /*! \brief Set the \a j th value of the \a i th component to \a x. */
      void set(const size_type& i, const MultiIndex& j, const I& x); 
      void set(const size_type& i, const MultiIndex& j, const R& x); 
      /*! \brief The \a j th value of the i th component. */
      const I& get(const size_type& i, const MultiIndex& j) const;
      /*! \brief Resize to a Taylor model in \a as variables with size \a rs image, order \a d and smoothness \a s. */
      void resize(const size_type& rs, const size_type& as, const size_type& d, const size_type& s);

      /*! \brief The \a i th component Taylor model. */
      TaylorModel<R> component(const size_type& i) const;

      /*! \brief Truncate the Taylor model to a Taylor model of order \a d and smoothness \a s within the domain \a domain. */
      TaylorModel<R> truncate(const size_type& order, const size_type& smoothness, const Geometry::Box<R>& domain) const;

      /*! \brief Evaluate the Taylor model at the point \a x. */
      LinearAlgebra::Vector<I> evaluate(const LinearAlgebra::Vector<R>& x) const;
      
      /*! \brief Evaluate the Taylor model at the point \a x. */
      LinearAlgebra::Vector<I> evaluate(const LinearAlgebra::Vector<I>& x) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      LinearAlgebra::Matrix<I> jacobian(const LinearAlgebra::Vector<I>& s) const;

       /*! \brief The domain of validity of the Taylor model. */
      Geometry::Box<R> domain() const;

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
      friend template<class R> TaylorModel<R> operator+(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Subtraction. */
      friend template<class R> TaylorModel<R> operator-(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Multiplication. At least one argument must be scalar-valued. */
      friend template<class R> TaylorModel<R> operator*(const TaylorModel<R1>&, const TaylorModel<R2>&);

      /*! \brief Multiplication by a scalar. */
      friend template<class R> TaylorModel<R> operator*(const R1&, const TaylorModel<R2>&);
      /*! \brief Multiplication by a scalar. */
      friend template<class R> TaylorModel<R> operator*(const TaylorModel<R1>&, const R2&);
      /*! \brief Division by a scalar. */
      friend template<class R> TaylorModel<R> operator/(const TaylorModel<R1>&, const R2&);

      /*! \brief Composition \f$p\circ q(x)=p(q(x))\f$. */
      friend template<class R> TaylorModel<R> compose(const TaylorModel<R1>&, const TaylorModel<R2>&);
      /*! \brief Power of a scalar Taylor model. */
      friend template<class R> TaylorModel<R> pow(const TaylorModel<R>& p, const unsigned int& n);
      /*! \brief Derivative with respect to variable \a k. */
      friend template<class R> TaylorModel<R> derivative(const TaylorModel<R>&, const size_type& k);
      /*! \brief Truncate within \a r to a Taylor model of order at most \a d, putting the error into terms of order \a s. */
      friend template<class R> TaylorModel<R> truncate(const TaylorModel<R>& p, const Rectangle<R>& bb, const size_type& d, const size_type& s);
#endif
     private:
      static void instantiate();
      array< array<I> > _powers(const LinearAlgebra::Vector<I>&) const;
      void _compute_jacobian() const;
      void _set_argument_size(const size_type& n);
      size_type _compute_maximum_component_size() const;
     private:
      //template<class R> friend void add(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void sub(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void mul(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void div(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void compose(TaylorModel<R>&,const TaylorModel<R>&,const TaylorModel<R>&);
      //template<class R> friend void scale(TaylorModel<R>&,const R&);
     private:
      /* Domain of definition. */
      Geometry::Box<R> _domain;
      /* The centre of the derivative expansion. */
      LinearAlgebra::Vector<R> _centre;
      size_type _smoothness; 
      Function::TaylorDerivative<I> _derivatives;
    };
    
    template<class R> TaylorModel<R> recentre(const TaylorModel<R>, const Geometry::Box<R>& bx, const Geometry::Point<R>& pt);

    template<class R> TaylorModel<R> add(const TaylorModel<R>&, const TaylorModel<R>&);
    template<class R> TaylorModel<R> sub(const TaylorModel<R>&, const TaylorModel<R>&);

  //template<class R> TaylorModel<R> mul(const TaylorModel<R>&, const R&);
  //template<class R> TaylorModel<R> div(const TaylorModel<R>&, const R&);

    template<class R> TaylorModel<R> derivative(const TaylorModel<R>&, const size_type&);
    template<class R> TaylorModel<R> compose(const TaylorModel<R>&, const TaylorModel<R>&);
    template<class R> TaylorModel<R> inverse(const TaylorModel<R>&, const LinearAlgebra::Vector<R>&);
    template<class R> TaylorModel<R> implicit(const TaylorModel<R>&, const LinearAlgebra::Vector<R>&);
    template<class R> TaylorModel<R> combine(const TaylorModel<R>&, const TaylorModel<R>&);
    template<class R> TaylorModel<R> join(const TaylorModel<R>&, const TaylorModel<R>&);

    template<class R> TaylorModel<R> truncate(const TaylorModel<R>&, const Geometry::Box<R>&, const size_type&, const size_type&);

    template<class R> TaylorModel<R> compose(const TaylorModel<R>&, const TaylorModel<R>&);
    template<class R> TaylorModel<R> inverse(const TaylorModel<R>&, const LinearAlgebra::Vector<R>&);
    template<class R> TaylorModel<R> implicit(const TaylorModel<R>&, const LinearAlgebra::Vector<R>&);
    template<class R> TaylorModel<R> derivative(const TaylorModel<R>&, const size_type& k);
    template<class R> TaylorModel<R> truncate(const TaylorModel<R>&, const Geometry::Box<R>&, const size_type&, const size_type&);

    template<class R> std::ostream& operator<<(std::ostream&, const TaylorModel<R>&);
  



  }
}

#include "taylor_model.inline.h"

#endif /* ARIADNE_TAYLOR_MODEL_H */
