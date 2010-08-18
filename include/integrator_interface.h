/***************************************************************************
 *            integrator_interface.h
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

/*! \file integrator_interface.h
 *  \brief Interface for solver classes for differential equations.
 */

#ifndef ARIADNE_INTEGRATOR_INTERFACE_H
#define ARIADNE_INTEGRATOR_INTERFACE_H


namespace Ariadne {

template<class T1,class T2> class Pair;
class Float;
class Interval;
template<class X> class Vector;
class VectorFunction;
class VectorTaylorFunction;

//! \ingroup SolverModule EvaluationModule
//! \brief Interface for integrating differential equations of the form \f$\dot{x}=f(x)\f$.
class IntegratorInterface
{
  protected:
    typedef Vector<Interval> IVector;
  public:
    //! \brief Virtual destructor.
    virtual ~IntegratorInterface() { };

    /*! \brief Make a dynamically-allocated copy. */
    virtual IntegratorInterface* clone() const = 0;

    /*! \brief Set the maximum allowable error in the flow. */
    virtual void set_maximum_error(double) = 0;

    /*! \brief Set the temporal order (if appropriate). */
    virtual void set_temporal_order(unsigned int) = 0;

    //! \brief Compute a pair \a (h,B) consisting of a bound \a B for the flow 
    //! starting in the \a state_domain for time step \a h.
    virtual Pair<Float,IVector> flow_bounds(const VectorFunction& vector_field,
                                            const IVector& parameter_domain,
                                            const IVector& state_domain,
                                            const Float& suggested_time_step) const = 0;

    //! \brief Solve \f$\dot{\phi}(x,t)=f(\phi(x,t))\f$.
    virtual VectorTaylorFunction flow(const VectorFunction& vector_field,
                                const IVector& state_domain,
                                const Float& time_domain) const = 0;


    //! \brief Solve \f$\dot{\phi}(a,x,t)=f(a,\phi(a,x,t))\f$.
    virtual VectorTaylorFunction flow(const VectorFunction& vector_field,
                                              const IVector& parameter_domain,
                                              const IVector& state_domain,
                                              const Float& time_domain) const = 0;

    //! \brief Compute \f$\phi(a,x,h)\f$, where \f$\dot{\phi}(a,x,t)=f(a,\phi(a,x,t))\f$.
    virtual VectorTaylorFunction time_step(const VectorFunction& vector_field,
                                     const IVector& parameter_domain,
                                     const IVector& state_domain,
                                     const Float& time_domain) const = 0;

};



} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_INTERFACE_H */
