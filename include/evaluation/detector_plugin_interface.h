/***************************************************************************
 *            detector_plugin_interface.h
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
 
/*! \file detector_plugin_interface.h
 *  \brief Methods for detecting crossings with constraints.
 */

#ifndef ARIADNE_DETECTOR_PLUGIN_INTERFACE_H
#define ARIADNE_DETECTOR_PLUGIN_INTERFACE_H

#include "../base/types.h"
#include "../base/declarations.h"
#include "../numeric/declarations.h"
#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"

#include "../evaluation/time_model.h"

namespace Ariadne {
  namespace Evaluation {

   
    /*! \brief %Base class for constraint crossing and detection schemes. 
     *  \ingroup Detection
     */
    template<class R>
    class DetectorPluginInterface 
    {
      typedef Numeric::Interval<R> I;
      typedef typename Geometry::Rectangle<R> BS;
     public:
      /*! \brief Virtual destructor. */
      virtual ~DetectorPluginInterface() { };

      /*! \brief Make a dynamically-allocated copy. */
      virtual DetectorPluginInterface<R>* clone() const = 0;

      /*! \brief Compute the value of a constraint over a set. */
      virtual Numeric::Interval<R> value(const Geometry::ConstraintInterface<R>& c, 
                                         const Geometry::Rectangle<R>& r) const = 0;

      /*! \brief Determine whether constraint \a c1 forces constraint \a c2 within \a dom.
       */
      virtual tribool forces(const Geometry::ConstraintInterface<R>& c1,
                             const Geometry::ConstraintInterface<R>& c2,
                             const Geometry::Rectangle<R>& dom) const = 0;

      /*! \brief Compute the normal derivative to of the vector field \a vf to the constraint \a c at the point \a pt.
       */
      virtual Numeric::Interval<R> normal_derivative(const System::VectorFieldInterface<R>& vf, 
                                                     const Geometry::DifferentiableConstraintInterface<R>& c, 
                                                     const Geometry::Point<I>& pt) const = 0;

      /*! \brief Estimate the time needed to cross a constraint. */
      virtual Numeric::Interval<R> crossing_time(const System::VectorFieldInterface<R>& vf, 
                                                 const Geometry::ConstraintInterface<R>& c, 
                                                 const Geometry::Point<I>& pt, 
                                                 const Geometry::Rectangle<R>& b) const = 0;

      /*! \brief Compute the value of the crossing time over a set. */
      virtual Evaluation::TimeModel<R> crossing_time(const System::VectorFieldInterface<R>& vf, 
                                                     const Geometry::ConstraintInterface<R>& c, 
                                                     const Geometry::Rectangle<R>& d, 
                                                     const Geometry::Rectangle<R>& b) const = 0;

    };

  }
}

#endif /* ARIADNE_DETECTOR_PLUGIN_INTERFACE_H */
