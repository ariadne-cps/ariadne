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
#include <iosfwd>

#include "utility/declarations.h"

namespace Ariadne {

struct FlowBoundsException : public std::runtime_error {
    FlowBoundsException(const StringType& what) : std::runtime_error(what) { }
};

//! \ingroup SolverModule EvaluationModule
//! \brief A solution to a differential equation could not be computed within the requested tolerances.
struct FlowTimeStepException : public std::runtime_error {
    FlowTimeStepException(const StringType& what) : std::runtime_error(what) { }
};

//! \ingroup SolverModule EvaluationModule
//! \brief Interface for integrating differential equations of the form \f$\dt{x}=f(x)\f$.
class IntegratorInterface
{
  public:
    //! \brief Virtual destructor.
    virtual ~IntegratorInterface() { };

    //! \brief Make a dynamically-allocated copy.
    virtual IntegratorInterface* clone() const = 0;

    //! \brief Write to an output stream.
    virtual Void write(OutputStream& os) const = 0;

    //! \brief Get the maximum allowable error in the flow.
    virtual double maximum_error() const = 0;
    //! \brief Set the maximum allowable error in the flow.
    virtual Void set_maximum_error(double) = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x)\f$ starting in \f$D\f$  for time step \f$h\leq h_{\max}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain and \f$h_{\max}\f$ is the \a maximum_time_step.
    virtual Pair<Float64Value,UpperBoxType>
    flow_bounds(const ValidatedVectorFunction& vector_field,
                const ExactBoxType& state_domain,
                const RawFloat64& maximum_time_step) const = 0;

    //! \brief Compute a validated version \f$\hat{\phi}\f$ of the flow \f$\phi(x,t)\f$ satisfying \f$\dt{\phi}(x,t)=f(\phi(x,t))\f$ for \f$x\in D\f$ and \f$t\in[0,h]\f$, where \f$h\f$ is a time step which is taken to be equal to \f$h_\mathrm{sug}\f$ if possible. The value of \f$h_\mathrm{sug}\f$ is overwritten with \f$h\f$, the actual time step used.
    //! <br>
    //! Returns: A validated version \f$\hat{\phi}\f$ of the flow over a short step represented as a single function over a box.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain, and \f$h_\mathrm{sug}\f$ is the \a suggested_time_step.
    virtual ValidatedVectorFunctionModel64
    flow_step(const ValidatedVectorFunction& vector_field,
              const ExactBoxType& state_domain,
              RawFloat64& suggested_time_step) const = 0;

    //! \brief Solve \f$\dt{\phi}(x,t)=f(\phi(x,t))\f$ for \f$x\in D\f$ and \f$t\in[0,h]\f$, assuming that the flow remains in \f$B\f$.
    //! If the flow does not remain in \f$B\f$, then \f$\hat{\phi}(x,t)\f$ may not be a bound for \f$\phi(x,t)\f$ if \f$\exists \tau\in[0,t],\ \phi(x,\tau)\not\in B\f$.
    //! If \f$(h,B)\f$ have been computed using the flow_bounds() method, then the flow is guarenteed to be correct for all \f$x\in D\f$ and \f$t\in[0,h]\f$.
    //! <br>
    //! Returns: A validated version \f$\hat{\phi}\f$ of the flow over a short step represented as a single function over a box.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain, \f$h\f$ is the \a time_step, and \f$B\f$ is the \a state_bounding_box.
    //! <br>
    //! Throws: A FlowTimeStepException if the flow cannot be computed sufficiently accurately for the given time step.
    virtual ValidatedVectorFunctionModel64
    flow_step(const ValidatedVectorFunction& vector_field,
              const ExactBoxType& state_domain,
              const Float64Value& time_step,
              const UpperBoxType& state_bounding_box) const = 0;

    //! \brief Solve \f$\dt{\phi}(x,t)=f(\phi(x,t))\f$ for initial conditions in \f$x\in D\f$ over the interval \f$[0,t_f]\f$.
    //! <br>
    //! Returns: A validated version \f$\hat{\phi}\f$ of the flow represented as a single function over a box.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain, and \f$t_f\f$ is the \a final_time.
    //! <br>
    //! Throws: A FlowTimeStepException if the flow cannot be represented as a single function to sufficiently accurately for the given time interval.
    virtual ValidatedVectorFunctionModel64
    flow_to(const ValidatedVectorFunction& vector_field,
            const ExactBoxType& state_domain,
            const Real& final_time) const = 0;

    //! \brief Solve \f$\dt{\phi}(x,t)=f(\phi(x,t))\f$ for initial conditions in \f$x\in D\f$ over the interval \f$[t_b,t_f]\f$.
    //! <br>
    //! Returns: A validated version of the flow represented as a list of functions whose spacial domains are all \f$D\f$ and whose time domains have union \f$[t_b,t_f]\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain, \f$t_b\f$ is the \a beginning_time, and \f$t_f\f$ is the \a final_time.
    virtual List<ValidatedVectorFunctionModel64>
    flow(const ValidatedVectorFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& beginning_time,
         const Real& final_time) const = 0;

    //! \brief Solve \f$\dt{\phi}(x,t)=f(\phi(x,t))\f$ for initial conditions in \f$x\in D\f$ over the interval \f$[0,t_f]\f$..
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain,  and \f$t_f\f$ is the \a final_time.
    virtual List<ValidatedVectorFunctionModel64>
    flow(const ValidatedVectorFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& final_time) const = 0;

    //! \brief Write to an output stream.
    friend inline OutputStream& operator<<(OutputStream& os, const IntegratorInterface& integrator) {
        integrator.write(os); return os;
}
};



} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_INTERFACE_H */
