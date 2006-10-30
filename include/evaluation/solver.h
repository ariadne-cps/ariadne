/***************************************************************************
 *            solver.h
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
 
/*! \file solver.h
 *  \brief Base class for solving equations.
 */

#ifndef _ARIADNE_SOLVER_H
#define _ARIADNE_SOLVER_H

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
   
    
    template<class R>
    class Solver
    {
     public:
      Solver(R max_error, uint max_steps)
        : _max_error(max_error), _max_steps(max_steps) { }
      
      virtual ~Solver();
        
      const R& maximum_error() const { return this->_max_error; }
      const uint& maximum_number_of_steps() const { return this->_max_steps; }
      
      void set_maximum_error(R max_error) { this->_max_error=max_error; }
      void maximum_number_of_steps(uint max_steps) { this->_max_steps=max_steps; }
      
      virtual Geometry::Rectangle<R>solve(const System::VectorField<R>& f, 
                                          const Geometry::Rectangle<R>& r) = 0;
     private:
      R _max_error;
      uint _max_steps;
    };
    
    template<class R> Solver<R>::~Solver() { }
  }
}


namespace Ariadne {
  namespace System {

    /*! \brief A class representing the difference of a Map and the identity.
     *
     * Useful for computing fixed points.
     */
    template<class R> class DifferenceMap
      : public VectorField<R>
    {
     public:
      /*!\brief Construct from a map \a f, which must have the same argument dimension as result dimension. */
      DifferenceMap(const Map<R>& f) : _base(f) { 
        if(f.argument_dimension()!=f.result_dimension()) { 
          throw IncompatibleDimensions("DifferenceMap<R>::DifferenceMap(Map<R>): "
                                       "The argument and result dimensions must be equal"); } }
      /*! \brief Make a copy (clone) of the vector field. */
      DifferenceMap<R>* clone() const { return new DifferenceMap<R>(this->_base); }
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


#endif /* _ARIADNE_SOLVER_H */
