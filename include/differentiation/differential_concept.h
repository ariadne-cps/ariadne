/***************************************************************************
 *            differential_concept.h
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
 
/*! \file _concept.h
 *  \brief Concept check class for automatic differentiation classes.
 */
 
#ifndef ARIADNE_DIFFERENTIAL_CONCEPT_H
#define ARIADNE_DIFFERENTIAL_CONCEPT_H

#include "boost/concept_check.hpp"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;

/*!\ingroup Differentiation
 * \brief Concept checking class for differntial variables.
 * 
 * A differential variable is a quantity depending on
 * many arguments that follows the rules of a 
 * differential algebra. 
 */
template<class D>
class DifferentialConcept {
 private:
  typedef typename D::value_type X;
  D& differential_ref;
  const D& differential;
  const X scalar;
  const double dbl;
 public:
  BOOST_CONCEPT_USAGE(DifferentialConcept)
  {   
    /*! \brief The scalar type used for values. */
    typedef typename D::value_type value_type;
    /*! \brief The type used for the index. */
    typedef typename D::index_type index_type;

    
    differential_ref = D::constant(0u,0u,scalar);
    differential_ref = D::variable(0u,0u,scalar,0u);
    
    differential_ref = scalar;
    
    differential_ref = differential_ref += differential;
    differential_ref = differential_ref -= differential;
    
    differential_ref = differential_ref += scalar;
    differential_ref = differential_ref -= scalar;
    differential_ref = differential_ref *= scalar;
    differential_ref = differential_ref /= scalar;
    
    differential_ref = differential_ref += dbl;
    differential_ref = differential_ref -= dbl;
    differential_ref = differential_ref *= dbl;
    differential_ref = differential_ref /= dbl;
    
    differential_ref = +differential;
    differential_ref = -differential;
    differential_ref = differential + differential;
    differential_ref = differential - differential;
    differential_ref = differential * differential;
    differential_ref = differential / differential;
    
    differential_ref = differential + scalar;
    differential_ref = differential - scalar;
    differential_ref = differential * scalar;
    differential_ref = differential / scalar;
    differential_ref = scalar + differential;
    differential_ref = scalar - differential;
    differential_ref = scalar * differential;
    differential_ref = scalar / differential;
  }
};
  
/*!\ingroup Differentiation
 * \brief Concept checking class for differntial variables.
 * 
 * A differential variable is a quantity depending on
 * many arguments that follows the rules of a 
 * differential algebra. 
 */
template<class DV>
class DifferentialVectorConcept {
 private:
  typedef typename DV::value_type D;
  typedef typename D::value_type X;
    DV& differential_vector_ref;
    const DV& differential_vector;
    D& differential_ref;
    const D& differential;
    Vector<X>& vector_ref;
    const Vector<X>& vector;
    Matrix<X>& matrix_ref;
    
    const X scalar;
    const double dbl;
    const Slice slice;
    const uint index;

 public:
  BOOST_CONCEPT_USAGE(DifferentialVectorConcept)
  {   
    typedef typename DV::value_type D;
    typedef typename D::value_type X;

    differential_vector_ref = project(differential_vector,slice);
    project(differential_vector_ref,slice) = differential_vector;
    differential_vector_ref = join(differential_vector,differential);
    differential_vector_ref = join(differential_vector,differential_vector);
    
    vector_ref = differential_vector.value();
    matrix_ref = differential_vector.jacobian();

    vector_ref = evaluate(differential_vector,vector);
    differential_vector_ref = translate(differential_vector,vector);
    differential_vector_ref = compose(differential_vector,differential_vector);
    differential_vector_ref = derivative(differential_vector,index);
    differential_vector_ref = antiderivative(differential_vector,index);
    differential_vector_ref = inverse(differential_vector);
    differential_vector_ref = implicit(differential_vector);
    differential_vector_ref = flow(differential_vector);
    differential_vector_ref = hitting(differential_vector,differential_vector);
  }
  
};

} // namespace Ariadne

#endif /* ARIADNE_DIFFERENTIAL_CONCEPT_H */

