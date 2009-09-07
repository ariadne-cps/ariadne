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
typedef double Float;
class Interval;
template<class X> class Vector;
class VectorFunctionInterface;
class VectorTaylorFunction;

/*! \ingroup \ingroup Solvers
 *  \brief %Common functionality for solving (nonlinear) equations.
 */
class IntegratorInterface
{
  protected:
    typedef Vector<Interval> IVector;
  public:
    //! \brief Virtual destructor.
    virtual ~IntegratorInterface() { };


    //! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt.
    virtual Pair<Float,IVector> flow_bounds(const VectorFunctionInterface& vector_field,
                                            const IVector& parameter_domain,
                                            const IVector& state_domain,
                                            const Float& suggested_time_step) const = 0;

    //! \brief Solve \f$\dot{\phi}(x,t)=f(\phi(x,t))\f$.
    virtual VectorTaylorFunction flow(const VectorFunctionInterface& vector_field,
                                const IVector& state_domain,
                                const Float& time_domain) const = 0;


    //! \brief Solve \f$\dot{\phi}(a,x,t)=f(a,\phi(a,x,t))\f$.
    virtual VectorTaylorFunction flow(const VectorFunctionInterface& vector_field,
                                              const IVector& parameter_domain,
                                              const IVector& state_domain,
                                              const Float& time_domain) const = 0;

    //! \brief Compute \f$\phi(a,x,h)\f$, where \f$\dot{\phi}(a,x,t)=f(a,\phi(a,x,t))\f$.
    virtual VectorTaylorFunction time_step(const VectorFunctionInterface& vector_field,
                                     const IVector& parameter_domain,
                                     const IVector& state_domain,
                                     const Float& time_domain) const = 0;

};



} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_INTERFACE_H */
