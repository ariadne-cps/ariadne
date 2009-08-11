/***************************************************************************
 *            integrator.h
 *
 *  Copyright  2006-9  Pieter Collins
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
    virtual Pair<Float,IVector> flow_bounds(const FunctionInterface& vector_field,
                                            const IVector& parameter_domain,
                                            const IVector& state_domain,
                                            const Float& suggested_time_step) const;

    virtual TaylorFunction time_step(const FunctionInterface& vector_field,
                                     const IVector& parameter_domain,
                                     const IVector& state_domain,
                                     const Float& suggested_time_step) const;

    virtual TaylorFunction flow(const FunctionInterface& vector_field,
                                const Vector<Interval>& state_domain,
                                const Float& suggested_time_step) const;

    using IntegratorInterface::flow;
};


class TaylorIntegrator
    : public IntegratorBase
{
  public:
    TaylorIntegrator(uint to) : _spacial_order(1), _temporal_order(to), _error(1e-16) { }

    uint temporal_order() const { return this->_temporal_order; }

    virtual TaylorFunction flow(const FunctionInterface& vector_field,
                                const Vector<Interval>& parameter_domain,
                                const Vector<Interval>& state_domain,
                                const Float& suggested_time_step) const;

    using IntegratorInterface::flow;

  private:
    uint _spacial_order;
    uint _temporal_order;
    double _error;
};




} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_H */
