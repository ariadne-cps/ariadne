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



class IntegratorBase
    : public IntegratorInterface
    , public Loggable
{
  public:
    IntegratorBase(uint to, double e) : _temporal_order(to), _maximum_error(e) { assert(e>0.0); }
    virtual void set_temporal_order(uint m) { this->_temporal_order=m; }
    virtual void set_maximum_error(double e) { assert(e>0.0); this->_maximum_error=e; }
    uint temporal_order() const { return this->_temporal_order; }
    double maximum_error() const  { return this->_maximum_error; }

    virtual Pair<Float,IntervalVector>
    flow_bounds(const VectorFunction& vector_field,
                const IntervalVector& state_domain,
                const Float& suggested_time_step) const;

    virtual VectorTaylorFunction
    flow_step(const VectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& suggested_time_step) const;

    virtual VectorTaylorFunction
    flow(const VectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Real& time) const;

    virtual VectorTaylorFunction
    flow(const VectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Interval& time_domain) const;

    virtual VectorTaylorFunction
    flow_step(const VectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& suggested_time_step,
              const IntervalVector& bounding_box) const = 0;

  public:
    uint _temporal_order;
    double _maximum_error;
};


class TaylorIntegrator
    : public IntegratorBase
{
    double _sweep_threshold;
  public:
    TaylorIntegrator() : IntegratorBase(4,1e-4), _sweep_threshold(1e-8) { }
    TaylorIntegrator(uint to, double e, double sw=0.0) : IntegratorBase(to,e), _sweep_threshold(sw) {
        if(_sweep_threshold==0.0) { _sweep_threshold=e; } }
    virtual TaylorIntegrator* clone() const { return new TaylorIntegrator(*this); }

    virtual VectorTaylorFunction
    flow_step(const VectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    using IntegratorBase::flow_step;
};




} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_H */
