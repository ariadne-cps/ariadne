/***************************************************************************
 *            halfspace.h
 *
 *  Copyright  2008 Pieter Collins
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

/*! \file halfspace.h
 *  \brief A halfspace in Euclidean space.
 */
 
#ifndef ARIADNE_HALFSPACE_H
#define ARIADNE_HALFSPACE_H

#include <iosfwd>
#include <vector>

#include "base/tribool.h"
#include "base/iterator.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"

namespace Ariadne {  
  

    /*! \ingroup BasicSet
     *  \brief A halfspace described by a linear inequalities.
     *
     *  The set is described as
     *  \f$ x\in\mathbb{R}^d \mid a^x\leq b , \f$
     *  where \f$A\f$ is a \f$d\f$-vector and \f$b\f$ is a scalar.
     */ 
    template<class X>
    class Halfspace {
      typedef typename traits<X>::number_type R;
      typedef typename traits<X>::arithmetic_type F;
     private:
      array<X> _data;
    public:
      //@{ 
      //! \name Constructors
      /*! \brief Construct the degenerate halfspace \f$\sum_{i=0}^{d-1} 0 x_i \geq 0 \f$ in \f$\R^d\f$. */
      explicit Halfspace(dimension_type d) : _data(d+1u) { }
      /*! \brief Construct the halfspace \f$\sum_{i=0}^{d-1} c_i x_i + c_d \geq 0 \f$. */
      template<class XX> explicit Halfspace(dimension_type d, const XX* c) : _data(c,c+d+1u) { }
      /*! \brief Construct the halfspace \f$\sum_{i=0}^{d-1} a_i x_i \leq b \f$. */
      template<class XX> explicit Halfspace(dimension_type d, const XX* a, const XX& b) : _data(d+1u) { 
        for(size_type i=0; i!=d; ++i) { this->_data[i]=-a[i]; } this->_data[d]=b; }
      //@}

      //@{
      //! \name Data access
      /*! A constant reference to the array of real data. */
      const array<X>& data() const { return this->_data; };
      //@}

      //@{
      //! \name Geometric operations
      /*! \brief The dimension the halfspace lies in. */
      dimension_type dimension() const { return this->_data.size()-1; }
      //@}
    };
    

    template<class X1, class X2>
    tribool 
    satisfies(const Point<X1>& pt, const Halfspace<X2>& hs)
    {
      ARIADNE_ASSERT(pt.dimension()==hs.dimension());
      typedef typename traits<X1,X2>::arithmetic_type F;
      F dot=hs.data()[hs.dimension()];
      for(dimension_type i=0; i!=hs.dimension(); ++i) {
        dot+=pt[i]*hs.data()[i];
      }
      return (dot>0 ? tribool(true) : dot<0 ? tribool(false) : indeterminate);
    }


    template<class X>
    std::ostream& operator<<(std::ostream& os, const Halfspace<X>& hsp) {
      for(dimension_type i=0; i!=hsp.dimension(); ++i) {
        os << (i==0?"[":",") << X(-hsp.data()[i]); }
      return os << ";" << hsp.data()[hsp.dimension()] << "]";
    }

  
} // namespace Ariadne

#endif /* ARIADNE_HALFSPACE_H */
