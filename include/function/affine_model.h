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
#include "geometry/point.h"
#include "geometry/box.h"

namespace Ariadne {
    
    template<class R> class FunctionInterface;
    template<class R> class AffineVariable;
  
    template<class R> class AffineModel;
    template<class R> AffineModel<R> add(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> sub(const AffineModel<R>&, const AffineModel<R>&);

    template<class R> AffineModel<R> recentre(const AffineModel<R>&, const Box<R>& bx, const Point<R>& pt);
    template<class R> AffineModel<R> restrict(const AffineModel<R>&, const Box<R>& bx);
    template<class R> AffineModel<R> reduce(const AffineModel<R>&, size_type);
    template<class R> AffineModel<R> translate(const AffineModel<R>&, const Point<R>& c);

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
      AffineModel(const Box<R>& d,
                  const Point<R>& c, 
                  const Point<I>& v, 
                  const Matrix<I>& j);
     
      /*! \brief Constructor. */
      AffineModel(const Box<R>& d,
                  const Point<R>& c, 
                  const array< AffineVariable<I> >& av);

      /*! \brief Constructor. */
      AffineModel(const Box<R>& d,
                  const Point<R>& c,
                  const FunctionInterface<R>& f);

      // Data access
      /*! \brief The image of the centre of the domain. */
      Point<I> const& value() const { return this->_value; }
      /*! \brief The jacobian derivative of the model. */
      Matrix<I> const& jacobian() const {  return this->_jacobian; }


      /*! \brief The size of the result. */
      size_type result_size() const { return this->_value.dimension(); }
      /*! \brief The size of the argument. */
      size_type argument_size() const { return this->_centre.dimension(); }
      /*! \brief The order of the function model. Returns the constant 1. */
      smoothness_type order() const { return 1u; }
      /*! \brief The smoothness of the function. */
      smoothness_type smoothness() const { return 1u; }

      /*! \brief The domain of validity of the model. */
      Box<R> domain() const { return this->_domain; }
      /*! \brief The centre of the  of validity of the model. */
      Point<R> centre() const { return this->_centre; }
      /*! \brief The image of the domain. */
      Box<R> range() const;

      /*! \brief Evaluate at a point. */
      Point<I> evaluate(const Point<I>& pt) const;
      Point<I> evaluate(const Point<R>& pt) const {
        return this->evaluate(Point<I>(pt)); }
      /*! \brief The Jacobian derivative at a point. */
      Point<I> jacobian(const Point<I>& pt) const;
      Point<I> jacobian(const Point<R>& pt) const {
        return this->evaluate(Point<I>(pt)); }

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      static void instantiate();
     private:
      Box<R> _domain;
      Point<R> _centre;
      Point<I> _value;
      Matrix<I> _jacobian;
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const AffineModel<R>& am) {
      return am.write(os);
    }

  
} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_H */
