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

#ifndef _ARIADNE_DISCRETE_TIME_SYSTEM_H
#define _ARIADNE_DISCRETE_TIME_SYSTEM_H

#include <string>

#include "../declarations.h"

namespace Ariadne {
  namespace System {

    /*! \brief Abstract base class for noisy control systems operating in discrete time.
     *  \ingroup System
     *  \ingroup DiscreteTime
     */
    template<class R>
    class DiscreteTimeSystem
    {
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief  An the image of a point. Only available if the image can be 
       *  computed exactly. */
      virtual Geometry::Point<R> 
      operator() (const Geometry::Point<R>& x,
                  const Geometry::Point<R>& u,
                  const Geometry::Point<R>& v) const;
      
      /*! \brief  The system operating on rectangular sets. */
      virtual Geometry::Rectangle<R> 
      operator() (const Geometry::Rectangle<R>& x,
                  const Geometry::Rectangle<R>& u,
                  const Geometry::Rectangle<R>& v) const;

      /*! \brief  The system acting on zonotopic sets. */
      virtual Geometry::Zonotope<R> 
      operator() (const Geometry::Zonotope<R>& x,
                  const Geometry::Zonotope<R>& u,
                  const Geometry::Zonotope<R>& v) const;

      /*! \brief  The dimension of the state space. */
      virtual dimension_type state_space_dimension() const;
      
      /*! \brief  The dimension of the control space. */
      virtual dimension_type control_space_dimension() const;
      
      /*! \brief  The dimension of the noise space. */
      virtual dimension_type noise_space_dimension() const;
      
      /*! \brief  The name of the system. */
      virtual std::string name() const { return "DiscreteTimeSystem"; }
    };

  }
}


#endif /* _ARIADNE_DISCRETE_TIME_SYSTEM_H */
