/***************************************************************************
 *            integrator_interface.h
 *
 *  Copyright  2006-10  Pieter Collins
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

/*! \file integrator_interface.h
 *  \brief Interface for solver classes for differential equations.
 */

#ifndef ARIADNE_INTEGRATOR_INTERFACE_H
#define ARIADNE_INTEGRATOR_INTERFACE_H

#include <string>

namespace Ariadne {

template<class T1,class T2> class Pair;
class Real;
class Float;
class Interval;
template<class X> class Vector;
typedef Vector<Interval> IntervalVector;

template<class X> class VectorFunction;
typedef VectorFunction<Real> RealVectorFunction;
class VectorTaylorFunction;

struct FlowBoundsException : public std::runtime_error {
    FlowBoundsException(const std::string& what) : std::runtime_error(what) { }
};

struct FlowTimeStepException : public std::runtime_error {
    FlowTimeStepException(const std::string& what) : std::runtime_error(what) { }
};

//! \ingroup SolverModule EvaluationModule
//! \brief Interface for integrating differential equations of the form \f$\dot{x}=f(x)\f$.
class IntegratorInterface
{
  public:
    //! \brief Virtual destructor.
    virtual ~IntegratorInterface() { };

    /*! \brief Make a dynamically-allocated copy. */
    virtual IntegratorInterface* clone() const = 0;

    /*! \brief Write to an output stream. */
    virtual void write(std::ostream& os) const = 0;

    /*! \brief Set the maximum allowable error in the flow. */
    virtual void set_maximum_error(double) = 0;

    //! \brief Compute a pair \a (h,B) consisting of a bound \a B for the flow
    //! starting in the \a state_domain for time step \a h.
    virtual Pair<Float,IntervalVector>
    flow_bounds(const RealVectorFunction& vector_field,
                const IntervalVector& state_domain,
                const Float& maximum_time_step) const = 0;

    //! \brief Solve \f$\dot{\phi}(x,t)=f(\phi(x,t))\f$ for \f$t\in[0,h]\f$ where \f$h\f$ is a time step based on \a suggested_time_step.
    //! The value of \a suggested_time_step is overwritten with the actual time step used.
    virtual VectorTaylorFunction
    flow_step(const RealVectorFunction& vector_field,
              const IntervalVector& state_domain,
              Float& suggested_time_step) const = 0;

    //! \brief Solve \f$\dot{\phi}(x,t)=f(\phi(x,t))\f$ for \f$t\in[0,h]\f$ where \f$h\f$ is the \a time_step used,
    //! and \a state_bounding_box is a bound for the trajectories.
    //! Throws a FlowTimeStepException if the flow cannot be computed sufficiently accurately for the given time step.
    virtual VectorTaylorFunction
    flow_step(const RealVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& state_bounding_box) const = 0;

    //! \brief Solve \f$\dot{\phi}(x,t)=f(\phi(x,t))\f$.
    virtual VectorTaylorFunction
    flow(const RealVectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Real& time) const = 0;

    //! \brief Solve \f$\dot{\phi}(x,t)=f(\phi(x,t))\f$.
    virtual VectorTaylorFunction
    flow(const RealVectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Interval& time_domain) const = 0;

};

inline std::ostream& operator<<(std::ostream& os, const IntegratorInterface& integrator) {
    integrator.write(os); return os;
}


} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_INTERFACE_H */
