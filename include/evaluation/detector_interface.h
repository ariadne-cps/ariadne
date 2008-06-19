/***************************************************************************
 *            detector_interface.h
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
 
/*! \file detector_interface.h
 *  \brief Methods for detecting crossings with constraints.
 */

#ifndef ARIADNE_DETECTOR_INTERFACE_H
#define ARIADNE_DETECTOR_INTERFACE_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "evaluation/time_model.h"

namespace Ariadne {
  

   
    /*! \ingroup EvaluatorInterfaces \ingroup Solvers
     *  \brief Interface for constraint crossing and detection schemes. 
     */
    template<class R>
    class DetectorInterface 
    {
      typedef Interval<R> I;
      typedef Box<R> BS;
     public:
      /*! \brief Virtual destructor. */
      virtual ~DetectorInterface() { };

      /*! \brief Make a dynamically-allocated copy. */
      virtual DetectorInterface<R>* clone() const = 0;

      /*! \brief Compute the value of a constraint over a set. */
      virtual Interval<R> value(const ConstraintInterface<R>& c, 
                                         const Box<R>& r) const = 0;

      /*! \brief Determine whether constraint \a c1 forces constraint \a c2 within \a dom.
       */
      virtual tribool forces(const ConstraintInterface<R>& c1,
                             const ConstraintInterface<R>& c2,
                             const Box<R>& dom) const = 0;

      /*! \brief Compute the normal derivative to of the vector field \a vf to the constraint \a c at the point \a pt.
       */
      virtual Interval<R> normal_derivative(const VectorField<R>& vf, 
                                                     const ConstraintInterface<R>& c, 
                                                     const Point<I>& pt) const = 0;

      /*! \brief Estimate the time needed to cross a constraint. */
      virtual Interval<R> crossing_time(const VectorField<R>& vf, 
                                                 const ConstraintInterface<R>& c, 
                                                 const Point<I>& pt, 
                                                 const Box<R>& b) const = 0;

      /*! \brief Compute the value of the crossing time over a set. */
      virtual TimeModel<R> crossing_time(const VectorField<R>& vf, 
                                                     const ConstraintInterface<R>& c, 
                                                     const Box<R>& d, 
                                                     const Box<R>& b) const = 0;

    };

  
} // namespace Ariadne

#endif /* ARIADNE_DETECTOR_INTERFACE_H */
