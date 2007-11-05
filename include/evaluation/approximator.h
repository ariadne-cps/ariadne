/***************************************************************************
 *            approximator.h
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
 
/*! \file approximator.h
 *  \brief Methods for approximating basic sets
 */

#ifndef ARIADNE_APPROXIMATOR_H
#define ARIADNE_APPROXIMATOR_H

#include "../base/types.h"
#include "../base/declarations.h"
#include "../numeric/declarations.h"
#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"

#include "approximator_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Geomerical approximation schemes.
     *  \ingroup Approximation
     */
    template<class BS>
    class Approximator
      : public ApproximatorInterface<BS>
    {
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Virtual destructor. */
      virtual ~Approximator();

      /*! \brief Default constructor. */
      Approximator();

      /*! \brief Copy constructor. */
      Approximator(const Approximator<BS>& approx);

      /*! \brief Make a dynamically-allocated copy. */
      virtual Approximator<BS>* clone() const;

      /*! \brief Computets and over-approximation of a set from a rectangle. */
      virtual BS over_approximation(const Geometry::Rectangle<R>& r) const;

      /*! \brief Computets and over-approximation of a set from a rectangle. */
      Geometry::GridCellListSet<R> outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const;

    };


    // Need a specialisation for Rectangle since outer_approximation returns a GridBlock
    template<class R>
    class Approximator< Geometry::Rectangle<R> >
      : public ApproximatorInterface< Geometry::Rectangle<R> >
    {
      typedef Numeric::Interval<R> I;
      typedef Geometry::Rectangle<R> BS;
     public:
      virtual Approximator<BS>* clone() const;
      virtual BS over_approximation(const Geometry::Rectangle<R>& r) const;
      Geometry::GridCellListSet<R> outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const;
    };


    /*! \brief Geomerical approximation schemes.
     *  \ingroup Approximation
     */
    template<class BS>
    class FastApproximator
      : public ApproximatorInterface<BS>
    {
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Virtual destructor. */
      virtual ~FastApproximator();

      /*! \brief Default constructor. */
      FastApproximator();

      /*! \brief Copy constructor. */
      FastApproximator(const FastApproximator<BS>& approx);

      /*! \brief Make a dynamically-allocated copy. */
      virtual FastApproximator<BS>* clone() const;

      /*! \brief Computets and over-approximation of a set from a rectangle. */
      virtual BS over_approximation(const Geometry::Rectangle<R>& r) const;

      /*! \brief Computets and over-approximation of a set from a rectangle. */
      Geometry::GridCellListSet<R> outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const;

    };

  }
}

#endif /* ARIADNE_APPROXIMATOR_H */
