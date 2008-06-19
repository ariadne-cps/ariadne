/***************************************************************************
 *            image_set.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file image_set.h
 *  \brief Sets defined by function and images.
  */

#ifndef ARIADNE_IMAGE_SET_H
#define ARIADNE_IMAGE_SET_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "base/tribool.h"

#include "function/function_interface.h"
#include "geometry/set_interface.h"

namespace Ariadne {
  
    
    //! \ingroup ExactSet
    /*! \brief A set defined by the conditions \f$f(x)\geq0\f$ for some function \f$f\f$. 
     *   Satisfies the conditions of the RegularSetInterface, 
     *   which means that only superset(Box) and disjoint(Box) need be meaningfully defined, 
     *   and only outer- and inner-approximations can be computed.
     */
    template<class R>
    class ImageSet
      : public SetInterface<R>
    {
      typedef typename traits<R>::arithmetic_type A;
     public:
      /*! \brief Construct the set \f$f(x)\in B\f$ from the function \f$f\f$ and the box \a D. */
      ImageSet(const Box<R> D, const FunctionInterface<R>& f)
        : _domain(D), _function_ptr(f.clone());

      /*! \brief Return a new dynamically-allocated copy of the set. */
      virtual ImageSet<R>* clone() const { return new ImageSet<R>(*this); }
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const { return this->_function_ptr->argument_size(); }
      /*! \brief Test if the set contains a point. */
      virtual tribool contains(const Point<R>& pt) const { return indeterminate; }
      /*! \brief */
        virtual tribool superset(const Box<R>& r) const { return indeterminate; }
      /*! \brief */
      virtual tribool intersects(const Box<R>& r) const { return indeterminate; }
      /*! \brief */
      virtual tribool disjoint(const Box<R>& r) const { return disjoint(this->_range(),r) or indeterminate; }
      /*! \brief */
      virtual tribool subset(const Box<R>& r) const { return subset(this->_range(),r) or indeterminate; }
      /*! \brief */
      virtual tribool bounded() const { return this->_range().bounded(); }
      /*! \brief */
      virtual Box<R> bounding_box() const { return this->_range(); }
      /*! \brief */
      virtual std::ostream& write(std::ostream& os) const;

      /*! \brief The codomain given the allowable values of the image function. */
      const Box<R>& domain() const { return this->_domain; }
      /*! \brief The function describing the images. */
      const FunctionInterface<R>& function() const { return *this->_function_ptr; }
     private:
      inline Box<R> _range() const { 
        return Box<R>(this->_function_ptr->evaluate(this->_domain.position_vectors())); }
     private:
      Box<R> _domain;
      boost::shared_ptr< const FunctionInterface<R> > _function_ptr;
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const ImageSet<R>& cset) {
      return cset.write(os);
    }
    




  }
}

#endif /* ARIADNE_IMAGE_SET_H */
