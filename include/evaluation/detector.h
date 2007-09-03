/***************************************************************************
 *            detector.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
    class Detector 
    {
     public:
      /*! \brief Virtual destructor. */
      virtual ~Detector();

      /*! \brief Copy constructor. */
      Detector(const Detector<R>& det);

      /*! \brief Make a dynamically-allocated copy. */
      virtual Detector<R>* clone() const;

      /*! \brief Compute the value of a constraint over a set. */
      virtual  Numeric::Interval<R> value(const Geometry::ConstraintInterface<R>& c, 
                                          const Geometry::BasicSetInterface<R>& bs);

      /*! \brief Compute the value of a constraint over a set. */
      virtual TimeModel<R> crossing_time(const System::VectorFieldInterface<R> vf, 
                                         const Geometry::ConstraintInterface<R>& c, 
                                         const Geometry::BasicSetInterface<R>& bs);

    };

  }
}

#endif /* ARIADNE_DETECTOR_H */
