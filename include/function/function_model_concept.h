/***************************************************************************
 *            function_model_concept.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file function_model_concept.h
 *  \brief Concept check class for function models.
 */
 
#ifndef ARIADNE_FUNCTION_MODEL_CONCEPT_H
#define ARIADNE_FUNCTION_MODEL_CONCEPT_H

#include "boost/concept_check.hpp"
//#include "boost/concept_archetype.hpp"
//#include "boost/concept/assert.hpp"
//#include "boost/concept/usage.hpp"

namespace Ariadne {

    template<class T> void same_type(T const&, T const&);
    template<class T, class X> void exact_type(X const& x) { T* t=0; same_type(t,x); }

    template<class R> class Interval;
    template<class X> class Vector;
    template<class X> class Matrix;
    template<class R> class FunctionInterface;

    /*!\ingroup Function
     * \brief Concept checking class for function models.
     * 
     * A function model is an approximation to a function on a 
     * restricted domain that can be manipulated numerically. 
     * Function model classes are the main workhorses of the 
     * Ariadne numerical engine. They contain methods for
     * evaluation, differentiation, solving implicit function
     * systems and integration.
     */
    template<class X>
    class FunctionModelConcept {
     public:
      /*! \brief The real number type. */
      typedef typename X::real_type real_type;
     private:
      typedef typename X::real_type R;
      typedef Interval<R> I;
     public:

      BOOST_CONCEPT_USAGE(FunctionModelConcept)
      {   
        X* model_ptr=0;
        FunctionInterface<R>* function_ptr=0;
        
        X& model=*model_ptr;
        FunctionInterface<R>& function=*function_ptr;
        Vector<R> point;
        Vector<I> domain;
        Vector<I> float_vector;
        Vector<I> interval_vector;
        Matrix<I> interval_matrix;
        ushort smoothness;
        uint accuracy;
        uint variable;
        R time;      
        
        //std::ostream* os_ptr=0;
        //std::ostream& os_ref=*os_ptr;

        X copy_model(model); /// Require copy construction.
        X new_model(domain,function,accuracy,smoothness); /// Require construction from a domain, a function and accuracy and smoothness parameters.
        domain = model.domain(); /// Require a domain of validity.
        interval_vector = model.evaluate(float_vector); /// Require evaluation at a point.
        interval_vector = model.evaluate(interval_vector); /// Require evaluation over a (small) box.
        interval_matrix = model.jacobian(float_vector); /// Require the Jacobian derivative at a point.
        interval_matrix = model.jacobian(interval_vector); /// Require the Jacobian derivative over a (small) box.

        // model = model + model; /// Require addition.
        // model = model - model; /// Require subtraction.

        model = restrict(model,domain); /// Require restriction to a subdomain.
        
        model = join(model,model); /// Require the direct sum of two models.
        model = derivative(model,variable); /// Require differentiation with respect to a variable.
        model = antiderivative(model,variable); /// Require antidifferentiation with respect to a variable.
        model = compose(model,model); /// Require composition of two models.
        model = inverse(model); /// Require inverse function two models.
        model = implicit(model); /// Require implicit function of two models.
        //model = flow(model); /// Require flow of an autonomous vector field.
        //model = hitting(model,model); /// Require hitting of a flow with a hypersurface.

        // os_ref = operator<<(os_ref, model);
      }

    };
  
} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_MODEL_CONCEPT_H */
