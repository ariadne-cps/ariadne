/***************************************************************************
 *            integrator.h
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

/*! \file integrator.h
 *  \brief Solver classes for differential equations.
 */

#ifndef ARIADNE_INTEGRATOR_H
#define ARIADNE_INTEGRATOR_H

#include <exception>
#include <stdexcept>
#include <string>

#include "integrator_interface.h"
#include "function_interface.h"

#include "logging.h"
#include "pointer.h"

namespace Ariadne {

class Real;
class Interval;
template<class X> class Vector;
template<class X> class Differential;
template<class X> class Procedure;
template<class X> class Polynomial;
typedef Vector< Procedure<Real> > RealVectorProcedure;
typedef Vector< Differential<Interval> > IntervalDifferentialVector;
class TaylorFunctionFactory;
typedef shared_ptr<const TaylorFunctionFactory> FunctionFactoryPointer;

class IntegratorBase
    : public IntegratorInterface
    , public Loggable
{
  public:
    IntegratorBase(double e);
    //! \brief A threshold for the error estimate of the approximation.
    virtual void set_maximum_error(double e) { assert(e>0.0); this->_maximum_error=e; }
    //! \brief The fraction L(f)*h used for a time step.
    //! The convergence of the Picard iteration is approximately Lf*h.
    void set_lipschitz_tolerance(double lt) { _lipschitz_tolerance = lt; }
    double maximum_error() const  { return this->_maximum_error; }

    //! \brief The class which constructs functions for representing the flow.
    const TaylorFunctionFactory& function_factory() const;
    //! \brief Set the class which constructs functions for representing the flow.
    void set_function_factory(const TaylorFunctionFactory& factory);


    virtual Pair<Float,IntervalVector>
    flow_bounds(const RealVectorFunction& vector_field,
                const IntervalVector& state_domain,
                const Float& suggested_time_step) const;

    virtual VectorTaylorFunction
    flow_step(const RealVectorFunction& vector_field,
              const IntervalVector& state_domain,
              Float& suggested_time_step) const;

    virtual VectorTaylorFunction
    flow(const RealVectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Real& time) const;

    virtual VectorTaylorFunction
    flow(const RealVectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Interval& time_domain) const;

    virtual VectorTaylorFunction
    flow_step(const RealVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& suggested_time_step,
              const IntervalVector& bounding_box) const = 0;

  public:
    double _maximum_error;
    double _lipschitz_tolerance;
    FunctionFactoryPointer _function_factory_ptr;
};


//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator
    : public IntegratorBase
{
    uint _maximum_temporal_order;
    double _sweep_threshold;
  public:
    TaylorPicardIntegrator() : IntegratorBase(1e-4), _sweep_threshold(1e-8) { }
    TaylorPicardIntegrator(uint maxto, double e, double sw=0.0) : IntegratorBase(e), _maximum_temporal_order(maxto), _sweep_threshold(sw) {
        if(_sweep_threshold==0.0) { _sweep_threshold=e; } }
    // Set the sweep threshold of the Taylor model.
    // Note: This method may be removed or modified in the future.
    //! \brief The order of the method in time.
    void set_maximum_temporal_order(uint m) { this->_maximum_temporal_order=m; }
    uint maximum_temporal_order() const { return this->_maximum_temporal_order; }
    void set_sweep_threshold(double lt) { _sweep_threshold = lt; }
    virtual TaylorPicardIntegrator* clone() const { return new TaylorPicardIntegrator(*this); }

    virtual VectorTaylorFunction
    flow_step(const RealVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    using IntegratorBase::flow_step;
};

//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class TaylorSeriesIntegrator
    : public IntegratorBase
{
    uint _maximum_temporal_order;
    uint _spacial_order;
    double _sweep_threshold;
  public:
    TaylorSeriesIntegrator() : IntegratorBase(1e-4), _maximum_temporal_order(6), _spacial_order(2), _sweep_threshold(1e-8) { }
    TaylorSeriesIntegrator(uint so, uint maxto, double e, double sw=0.0) : IntegratorBase(e), _maximum_temporal_order(maxto), _spacial_order(so), _sweep_threshold(sw) {
        if(_sweep_threshold==0.0) { _sweep_threshold=e; } }
    // Set the sweep threshold of the Taylor model.
    // Note: This method may be removed or modified in the future.
    //! \brief The order of the method in time.
    void set_spacial_order(uint m) { this->_spacial_order=m; }
    void set_maximum_temporal_order(uint n) { this->_maximum_temporal_order=n; }
    uint maximum_temporal_order() const { return this->_maximum_temporal_order; }
    void set_sweep_threshold(double lt) { _sweep_threshold = lt; }
    virtual TaylorSeriesIntegrator* clone() const { return new TaylorSeriesIntegrator(*this); }

    virtual VectorTaylorFunction
    flow_step(const RealVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    virtual Vector< Differential<Interval> >
    flow_diff(const RealVectorProcedure& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    using IntegratorBase::flow_step;
};


//! \brief An integrator computes a approximation to the flow which is affine in space.
//! \internal This code is written to allow higher-spacial order approximations.
class AffineIntegrator
    : public IntegratorBase
{
    uint _spacial_order;
    uint _temporal_order;
  public:
    AffineIntegrator(double e)
        : IntegratorBase(e), _spacial_order(1u), _temporal_order(4u) { }
    AffineIntegrator(double e, uint temporal_order)
        : IntegratorBase(e), _spacial_order(1u), _temporal_order(temporal_order) { }
    AffineIntegrator(double e, uint spacial_order, uint temporal_order)
        : IntegratorBase(e), _spacial_order(spacial_order), _temporal_order(temporal_order) { }
    //! \brief The order of the method in space.
    uint spacial_order() const { return this->_spacial_order; }
    //! \brief The order of the method in time.
    uint temporal_order() const { return this->_temporal_order; }
    virtual AffineIntegrator* clone() const { return new AffineIntegrator(*this); }

    virtual VectorTaylorFunction
    flow_step(const RealVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    using IntegratorBase::flow_step;

    //! \brief Compute the derivative of the flow of f at time zero within \a dom.
    IntervalDifferentialVector
    flow_derivative(const IntervalVectorFunction& f,
                    const IntervalVector& dom) const;
};


} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_H */
