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

#ifndef ARIADNE_SOLVER_H
#define ARIADNE_SOLVER_H

#include <exception>
#include <stdexcept>
#include <string>

#include "../declarations.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../system/vector_field.h"
#include "../system/map.h"

namespace Ariadne {
  namespace System {

    /*!
     * \brief A class representing the difference of a Map and the identity.
     *
     * Useful for computing fixed points.
     */
    template<class R> class DifferenceMap
      : public VectorField<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*!\brief Construct from a map \a f, which must have the same argument dimension as result dimension. */
      DifferenceMap(const Map<R>& f) : _base(f) { 
        if(f.argument_dimension()!=f.result_dimension()) { 
          throw Geometry::IncompatibleDimensions(__PRETTY_FUNCTION__); } }
      /*! \brief Make a copy (clone) of the vector field. */
      DifferenceMap<R>* clone() const { return new DifferenceMap<R>(this->_base); }
      /*!\brief The dimension of the space the map acts on. */
      virtual size_type smoothness() const { return _base.smoothness(); }
      /*!\brief The dimension of the space the map acts on. */
      virtual dimension_type dimension() const { return _base.argument_dimension(); }
      /*!\brief Evaluate the function \f$f(x)-x\f$, where \f$f\f$ is the map used to construct the difference map. */
      virtual LinearAlgebra::Vector<F> image(const Geometry::Point<F>& p) const {
        return _base.image(p)-p; }
      /*!\brief Evaluate the derivative of function \f$f(x)-x\f$, which is \f$Df(x)-I\f$. */
      virtual LinearAlgebra::Matrix< Interval<R> > jacobian(const Geometry::Point<F>& p) const {
        LinearAlgebra::Matrix<F> d=_base.jacobian(p);
        LinearAlgebra::Matrix<F> i=LinearAlgebra::Matrix< Interval<R> >::identity(this->dimension());
        return d-i; }
      /*!\brief The name of the class. */
      virtual std::string name() const { return "DifferenceMap"; }
     private:
      const Map<R>& _base;
    };
  }
}



namespace Ariadne {
  namespace Evaluation {
    
    
    
    /*!\ingroup Solve
     * \brief %Base class for solving (nonlinear) equations. 
     */
    template<class R>
    class Solver
    {
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief Constructor. */
      Solver(R max_error, uint max_steps)
        : _max_error(max_error), _max_steps(max_steps) { }
      
      /*! \brief Virtual destructor. */
      virtual ~Solver();
        
      /*! \brief The maximum permissible error of the solution. */
      const R& maximum_error() const { return this->_max_error; }
      /*! \brief The maximum number of steps allowed before the method must quit. */
      const uint& maximum_number_of_steps() const { return this->_max_steps; }
      
      /*! \brief Set the maximum error. */
      void set_maximum_error(R max_error) { this->_max_error=max_error; }
      /*! \brief Set the maximum number of steps. */
      void set_maximum_number_of_steps(uint max_steps) { this->_max_steps=max_steps; }
      
      /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
      virtual Geometry::Point<I>
      solve(const System::VectorField<R>& f,const Geometry::Point<I>& pt) = 0;
      /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
      virtual Geometry::Point<I> 
      fixed_point(const System::Map<R>& f,const Geometry::Point<I>& pt) {
        return this->solve(System::DifferenceMap<R>(f),pt); }
      
     private:
      R _max_error;
      uint _max_steps;
    };
    
    template<class R> Solver<R>::~Solver() { }
  }
}



#endif /* ARIADNE_SOLVER_H */
