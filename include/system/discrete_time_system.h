/***************************************************************************
 *            discrete_time_system.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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
 
/*! \file discrete_time_system.h
 *  \brief Discrete-time systems of the \f$x_{n+1}=f(x_n,u_n,v_n)\f$.
 */

#ifndef ARIADNE_DISCRETE_TIME_SYSTEM_H
#define ARIADNE_DISCRETE_TIME_SYSTEM_H

#include <string>

#include "geometry/declarations.h"

namespace Ariadne {
  

    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief Abstract base class for noisy control systems operating in discrete time.
     * 
     * A discrete-time control system is defined by a continuous function
     * \f$f:X\times U\times V\rightarrow X\f$ where \f$X\f$ is the state space,
     * \f$U\f$ is the control input space and \f$V\f$ is the noise.
     * The set of admissible controls \f$U(x)\f$ may depend on \f$x\f$.
     * 
     * To compute with a discrete-time control system, the multivalued mapping
     * \f$X\Rightarrow X\times U,\ x\mapsto\{x\}\times U(x)\f$ must be 
     * lower-semicomputable, and \f$X\times U\Rightarrow X\f$, \f$(x,u)\mapsto f(x,u,v)\f$
     * must be upper-semicomputable.
     *
     * Currently, a discrete-time control system is defined by operator()(const Point<F>&, const Point<F>&, const Point<F>&) const 
     * acting on fuzzy points in \f$X\times U\times V\f$, and by derivatives (see the Map class for details).
     * In a future version, restrictions on \f$U(x)\f$ will also be allowed.
     */
    template<class R>
    class DiscreteTimeControlSystemInterface
    {
      typedef typename traits<R>::arithmetic_type F;
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Point<R> state_type;
      
      /*! \brief Virtual destructor.  */
      virtual ~DiscreteTimeControlSystemInterface();
      
      /*! \brief  The image of a point, computed approximately.  */
      Point<F> 
      operator() (const Point<F>& x,
                  const Point<F>& u,
                  const Point<F>& v) const
      {
        return this->image(x,u,v);
      }
      
      /*! \brief  The image of a point, computed approximately.  */
      virtual 
      Point<F> 
      image(const Point<F>& x,
            const Point<F>& u,
            const Point<F>& v) const = 0;
      
      /*! \brief  The jacobian derivate of a point with respect to the state, computed approximately.  */
      virtual
      Matrix<F> 
      jacobian(const Point<F>& x,
               const Point<F>& u,
               const Point<F>& v) const;
      
      /*! \brief  The dimension of the state space. */
      virtual dimension_type state_space_dimension() const = 0;
      
      /*! \brief  The dimension of the control space. */
      virtual dimension_type control_space_dimension() const = 0;
      
      /*! \brief  The dimension of the noise space. */
      virtual dimension_type noise_space_dimension() const = 0;
      
      /*! \brief  The name of the system. */
      virtual std::string name() const { return "DiscreteTimeSystem"; }
    };

  
} // namespace Ariadne


#endif /* ARIADNE_DISCRETE_TIME_SYSTEM_H */
