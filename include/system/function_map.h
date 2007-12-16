/***************************************************************************
 *            function_map.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it
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
 
/*! \file function_map.h
 *  \briefMapInterfaces described by function objects.
 */

#ifndef ARIADNE_FUNCTION_MAP_H
#define ARIADNE_FUNCTION_MAP_H

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "base/types.h"
#include "numeric/declarations.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"
#include "system/map.h"

#include "function/function_interface.h"

namespace Ariadne {
  namespace System {

    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief A map (discrete-time dynamical system) described by an object satisfying the FunctionInterface.
     */
    template<class R>
    class FunctionMap
      : public MapInterface<R>
    {
     protected:
      typedef typename Numeric::traits<R>::arithmetic_type A; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      /*! \brief The type obtained by applying the map to a state. */
      typedef Geometry::Point<A> result_type;
      
      /*! \brief Construct from a function interface and parameters. */
      FunctionMap(const Function::DifferentiableFunctionInterface<R>& f, const Geometry::Point<A>& param);
      /*! \brief Construct from a function and parameters. */
      FunctionMap(const Function::InterpretedFunction<R>& f, const Geometry::Point<A>& param);

      /*! \brief Make a copy (clone) of the map. */
      virtual FunctionMap<R>* clone() const;

      /*! \brief An over-approximation to the image of a point. */
      virtual Geometry::Point<A> image(const Geometry::Point<A>& pt) const;
      /*! \brief The derivative of the \a i th component with respect to the multi-index j. */
      virtual A derivative(const Geometry::Point<A>& x, const size_type& i, const Function::MultiIndex& j) const;
      /*! \brief The Jacobian derivative matrix over a rectangle. */
      virtual LinearAlgebra::Matrix<A> jacobian(const Geometry::Point<A>& x) const;
        
      /*! \brief The dimension of the range space. */
      virtual dimension_type number_of_parameters() const;
      /*! \brief Set the parameters of the function. */
      virtual void set_parameters(const Geometry::Point<A>&) ;
      /*! \brief The parameters of the function. */
      virtual const Geometry::Point<A>& parameters() const;

      /*! \brief The degree of differentiability of the map. */
      virtual smoothness_type smoothness() const;
      /*! \brief The dimension of the domain space. */
      virtual dimension_type argument_dimension() const;
      /*! \brief The dimension of the range space. */
      virtual dimension_type result_dimension() const;
    
      /*! \brief The name of the map. */
      virtual std::string name() const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      boost::shared_ptr< Function::DifferentiableFunctionInterface<R> > _function_ptr;
      Geometry::Point<A> _parameters;
    };
   
  
  }
}

#endif /* ARIADNE_FUNCTION_MAP_H */
