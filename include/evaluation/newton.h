/***************************************************************************
 *            newton.h
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
 
/*! \file newton.h
 *  \brief Newton and interval Newton methods.
 */

#ifndef _ARIADNE_NEWTON_H
#define _ARIADNE_NEWTON_H

#include <exception>
#include <stdexcept>
#include <string>

#include "../declarations.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "system/map.h"
#include "system/vector_field.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*! \brief %Base class for exceptions in the Evaluation module. */
    class EvaluationException
      : public std::exception 
    {
     public:
       EvaluationException(const std::string& s) : _what(s) { }
      ~EvaluationException() throw () { }
      const char* what() const throw () { return this->_what.c_str(); }
     private:
      std::string _what;
    };
    
    /*! \brief Interval Newton solver.
     *  \ingroup Solve
     */
    template<typename R>
    Geometry::Rectangle<R>
    interval_newton(const System::VectorField<R>& f, 
                    const Geometry::Rectangle<R>& r, 
                    const R& e,
                    uint max_steps=64);

  }
}


namespace Ariadne {
  namespace System {

    /*! \brief A class representing the difference of a Map and the identity.
     *
     * Useful for computing fixed points.
     */
    template<typename R> class DifferenceMap
      : public VectorField<R>
    {
     public:
      /*!\brief Construct from a map \a f, which must have the same argument dimension as result dimension. */
      DifferenceMap(const Map<R>& f) : _base(f) { assert(f.argument_dimension()==f.result_dimension()); }
      /*!\brief The dimension of the space the map acts on. */
      virtual size_type smoothness() const { return _base.smoothness(); }
      /*!\brief The dimension of the space the map acts on. */
      virtual dimension_type dimension() const { return _base.argument_dimension(); }
      /*!\brief Evaluate the function \f$f(x)-x\f$, where \f$f\f$ is the map used to construct the difference map. */
      LinearAlgebra::Vector< Interval<R> > image(const Geometry::Rectangle<R>& r) const {
        return _base.image(r)-r; }
      /*!\brief Evaluate the derivative of function \f$f(x)-x\f$, which is \f$Df(x)-I\f$. */
      virtual LinearAlgebra::Matrix< Interval<R> > jacobian(const Geometry::Rectangle<R>& r) const {
        LinearAlgebra::Matrix< Interval<R> > d=_base.jacobian(r);
        LinearAlgebra::Matrix< Interval<R> > i=LinearAlgebra::Matrix< Interval<R> >::identity(this->dimension());
        return d-i; }
      /*!\brief The name of the class. */
      virtual std::string name() const { return "DifferenceMap"; }
     private:
      const Map<R>& _base;
    };
  }
}


#endif /* _ARIADNE_NEWTON_H */
