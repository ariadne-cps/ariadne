/***************************************************************************
 *            affine_model.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
/*! \file system/affine_model.h
 *  \brief An affine approximation to a function.
 */
 
#ifndef ARIADNE_AFFINE_MODEL_H
#define ARIADNE_AFFINE_MODEL_H

#include "base/array.h"

#include "numeric/traits.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/function_model_concept.h"

namespace Ariadne {
    
    template<class R> class FunctionInterface;
    template<class R> class AffineVariable;
  
    template<class R> class AffineModel;
    template<class R> AffineModel<R> add(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> sub(const AffineModel<R>&, const AffineModel<R>&);

    template<class R> AffineModel<R> recentre(const AffineModel<R>&, const Vector< Interval<R> >& bx, const Vector<R>& pt);
    template<class R> AffineModel<R> restrict(const AffineModel<R>&, const Vector< Interval<R> >& bx);
    template<class R> AffineModel<R> reduce(const AffineModel<R>&, size_type);
    template<class R> AffineModel<R> translate(const AffineModel<R>&, const Vector<R>& c);

    template<class R> AffineModel<R> compose(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> inverse(const AffineModel<R>&);
    template<class R> AffineModel<R> implicit(const AffineModel<R>&);


    /*!\ingroup FunctionModel
     * \brief Concrete class for functions.
     */
    template<class R>
    class AffineModel
    {
      typedef typename traits<R>::arithmetic_type A; 
      typedef typename traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;

      /*! \brief Destructor. */
      ~AffineModel() { };
     
      /*! \brief Constructor. */
      AffineModel(const Vector<I>& d,
                  const Vector<R>& c, 
                  const Vector<I>& v, 
                  const Matrix<I>& j);
     
      /*! \brief Constructor. */
      AffineModel(const Vector<I>& d,
                  const Vector<R>& c, 
                  const Vector< AffineVariable<I> >& av);

      /*! \brief Constructor. */
      AffineModel(const Vector<I>& d,
                  const Vector<R>& c,
                  const FunctionInterface<R>& f);

      // Data access
      /*! \brief The image of the centre of the domain. */
      Vector<I> const& value() const { return this->_value; }
      /*! \brief The jacobian derivative of the model. */
      Matrix<I> const& jacobian() const {  return this->_jacobian; }


      /*! \brief The size of the result. */
      size_type result_size() const { return this->_value.size(); }
      /*! \brief The size of the argument. */
      size_type argument_size() const { return this->_centre.size(); }
      /*! \brief The order of the function model. Returns the constant 1. */
      smoothness_type order() const { return 1u; }
      /*! \brief The smoothness of the function. */
      smoothness_type smoothness() const { return 1u; }

      /*! \brief The domain of validity of the model. */
      Vector<I> domain() const { return this->_domain; }
      /*! \brief The centre of the  of validity of the model. */
      Vector<R> centre() const { return this->_centre; }
      /*! \brief The image of the domain. */
      Vector<I> range() const;

      /*! \brief Evaluate at a point. */
      Vector<I> evaluate(const Vector<I>& pt) const;
      Vector<I> evaluate(const Vector<R>& pt) const {
        return this->evaluate(Vector<I>(pt)); }
      /*! \brief The Jacobian derivative at a point. */
      Vector<I> jacobian(const Vector<I>& pt) const;
      Vector<I> jacobian(const Vector<R>& pt) const {
        return this->evaluate(Vector<I>(pt)); }

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      static void instantiate();
     private:
      Vector<I> _domain;
      Vector<R> _centre;
      Vector<I> _value;
      Matrix<I> _jacobian;
     private:
      // Doesn't conform to FunctionModelConcept since no derivative/antiderivative
      // BOOST_CONCEPT_ASSERT((FunctionModelConcept< AffineModel<R> >));
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const AffineModel<R>& am) {
      return am.write(os);
    }



} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
