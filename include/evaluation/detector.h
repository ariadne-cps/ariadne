/***************************************************************************
 *            detector.h
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
 
/*! \file detector.h
 *  \brief Methods for detecting crossings with constraints.
 */

#ifndef ARIADNE_DETECTOR_H
#define ARIADNE_DETECTOR_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "geometry/constraint_interface.h"
#include "evaluation/detector_interface.h"

namespace Ariadne {
  

    template<class R> class TimeModel;

    /*! \ingroup Solvers
     *  \brief A class for detecting constraint crossings. 
     */
    template<class R>
    class Detector
      : public DetectorInterface<R>
    {
      typedef Interval<R> I;
     public:
      /*! \brief Virtual destructor. */
      virtual ~Detector();

      /*! \brief Default constructor. */
      Detector();

      /*! \brief Copy constructor. */
      Detector(const Detector<R>& det);

      /*! \brief Make a dynamically-allocated copy. */
      virtual Detector<R>* clone() const;

      /*! \brief Test if a set entirely satisfies the constraint. */
      template<class BS>
      tribool satisfies(const BS& bs, const ConstraintInterface<R>& c) const;


      /*! \brief Compute the value of a constraint over a rectangle. */
      Interval<R> value(const ConstraintInterface<R>& c, 
                                 const Box<R>& r) const;

      /*! \brief Compute the value of a constraint over a zonotope. */
      Interval<R> value(const ConstraintInterface<R>& c, 
                                 const Zonotope<R>& z) const;

      /*! \brief Determine whether constraint \a c1 forces constraint \a c2 within \a dom.
       */
      virtual tribool forces(const ConstraintInterface<R>& c1,
                             const ConstraintInterface<R>& c2,
                             const Box<R>& dom) const;

      /*! \brief Compute the normal derivative to of the vector field \a vf to the constraint \a c at the point \a pt.
       */
      virtual Interval<R> normal_derivative(const VectorField<R>& vf, 
                                                     const ConstraintInterface<R>& c, 
                                                     const Point<I>& pt) const;

      /*! \brief Estimate the time needed for the point \a pt to reach constraint \a c under vector field \a vf,
       *  assuming that the flow remains in \a bb. 
       */
      virtual Interval<R> crossing_time(const VectorField<R>& vf, 
                                                 const ConstraintInterface<R>& c, 
                                                 const Point<I>& pt, 
                                                 const Box<R>& bb) const;

      /*! \brief Compute the time needed for points in the domain rectangle \a dom to reach constraint \a c under vector field \a vf,
       *  assuming that the flow remains in \a bb. The integrator \a i is used to integrate the flow.
       */
      virtual TimeModel<R> crossing_time(const VectorField<R>& vf, 
                                                     const ConstraintInterface<R>& c, 
                                                     const Box<R>& dom, 
                                                     const Box<R>& bb) const;

    };



template<class R> template<class BS> inline
tribool 
Detector<R>::satisfies(const BS& bs,
                                   const ConstraintInterface<R>& c) const
{
  Interval<R> ivl=this->value(c,bs);
  Comparison cmp=c.comparison();
  if(ivl.upper()<0) {
    return (cmp==less);
  } else if(ivl.lower()>0) {
    return (cmp==greater);
  } else {
    return indeterminate;
  }
}

} // namespace Ariadne

#endif /* ARIADNE_DETECTOR_H */
