/***************************************************************************
 *            impact_system.h
 *
 *  Copyright  2008, Pieter Collins
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
 
/*! \file impact_system.h
 *  \brief Vector_type field interface.
 */
 
#ifndef ARIADNE_IMPACT_SYSTEM_H
#define ARIADNE_IMPACT_SYSTEM_H

#include <limits>

#include "base/types.h"
#include "numeric/declarations.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"

#include <boost/shared_ptr.hpp>
#include "function/function_interface.h"
#include "geometry/euclidean_space.h"


namespace Ariadne {



    /*!\ingroup System
     * \ingroup ContinuousTime
     * \brief A class representing a system in which continuous evolution is interrupted 
     * by discrete resets in a single mode.
     *
     * An impact system is a simple type of hybrid system. The state space is 
     * Euclidean space i.e. a single mode, and there is one guard set with an
     * associated reset map.
     *
     */
    template<class R>
    class ImpactSystem {
     protected:
      typedef typename traits<R>::arithmetic_type F; 
      typedef typename traits<R>::interval_type I; 
     public:
      /*! \brief The type used to represent time. */
      typedef Rational time_type;
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type used to describe the state space. */
      typedef EuclideanSpace state_space_type;
      
      /*! \brief Destructor. */
      ~ImpactSystem();

      /*! \brief Construct from functions giving the vector field \a vf, the guard set \a g and the impact (reset) map \a h. 
       *  The guard condition is a function \f$\R^n\rightarrow\R\f$ which is positive when continuous evolution is allowed.
       */
      ImpactSystem(const FunctionInterface<R>& vf, const FunctionInterface<R>& g, const FunctionInterface<R>& h);

      /*! \brief Make a copy (clone) of the vector field. */
      VectorField<R>* clone() const { return new VectorField<R>(*this); }
     
      /*! \brief The function defining the vector field. */
      FunctionInterface<R>& vector_field() const { return *this->_vector_field_ptr; }

      /*! \brief The function defining the vector field. */
      FunctionInterface<R>& guard_condition() const { return *this->_guard_ptr; }

      /*! \brief The function defining the vector field. */
      FunctionInterface<R>& impact_map() const { return *this->_impact_ptr; }

      /*! \brief The state space of the vector field. */
      EuclideanSpace state_space() const { return EuclideanSpace(this->_function_ptr->result_size()); }

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      boost::shared_ptr< FunctionInterface<R> > _vector_field_ptr;
      boost::shared_ptr< FunctionInterface<R> > _guard_ptr;
      boost::shared_ptr< FunctionInterface<R> > _impact_ptr;
    };
   
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ImpactSystem<R>& sys) {
      return sys.write(os);
    }
    
  
} // namespace Ariadne

#endif /* ARIADNE_VECTOR_FIELD_H */
