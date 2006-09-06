/***************************************************************************
 *            ppl_polyhedron.h
 *
 *  Thu Jan 27 10:26:36 2005
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

/*! \file ppl_polyhedron.h
 *  \brief Wrapper for Parmal Polyhedral Library Polyhedron.
 */
 
#ifndef _ARIADNE_PPL_POLYHEDRON_H
#define _ARIADNE_PPL_POLYHEDRON_H

#include <iosfwd>
#include <vector>

#include "../declarations.h"

// For forward declaration only
#include "../numeric/rational.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"

namespace Parma_Polyhedra_Library {
  class Polyhedron;
  class C_Polyhedron;
  class NNC_Polyhedron;
  class Constraint_System;
  class Generator_System;
}

namespace Ariadne {  
  namespace Geometry {

    Parma_Polyhedra_Library::NNC_Polyhedron 
    ppl_open_polyhedron(const LinearAlgebra::Matrix<Rational>& A, 
                        const LinearAlgebra::Vector<Rational>& b);
    
    Parma_Polyhedra_Library::NNC_Polyhedron 
    ppl_nnc_polyhedron(const LinearAlgebra::Matrix<Rational>& G);
    
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::Matrix<Rational>& A, 
                   const LinearAlgebra::Vector<Rational>& b);
    
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::IntervalVector<Rational>& c, 
                   const LinearAlgebra::Matrix<Rational>& A);
    
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::Vector<Rational>& c, 
                   const LinearAlgebra::Matrix<Rational>& A);
    
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::Matrix<Rational>& G);
    
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::IntervalVector<Rational>& r);
    
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::Vector<Rational>& p);
  
    
    LinearAlgebra::Matrix<Rational> 
    generators(const Parma_Polyhedra_Library::C_Polyhedron& A);
    
    LinearAlgebra::Matrix<Rational> 
    constraints_matrix(const Parma_Polyhedra_Library::C_Polyhedron& A);
    
    LinearAlgebra::Vector<Rational> 
    constraints_vector(const Parma_Polyhedra_Library::C_Polyhedron& A);
    
   
    dimension_type dimension(const Parma_Polyhedra_Library::C_Polyhedron& A);
                  
    bool empty(const Parma_Polyhedra_Library::C_Polyhedron& A);
                  
    bool empty_interior(const Parma_Polyhedra_Library::C_Polyhedron& A);
                  
    bool equal(const Parma_Polyhedra_Library::C_Polyhedron& A,
               const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    bool disjoint(const Parma_Polyhedra_Library::C_Polyhedron& A,
                  const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    bool interiors_intersect(const Parma_Polyhedra_Library::C_Polyhedron& A,
                             const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    bool inner_subset(const Parma_Polyhedra_Library::C_Polyhedron& A,
                      const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    bool subset(const Parma_Polyhedra_Library::C_Polyhedron& A,
                const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    
    
    Parma_Polyhedra_Library::NNC_Polyhedron
    interior(const Parma_Polyhedra_Library::C_Polyhedron& A);
    
    Parma_Polyhedra_Library::NNC_Polyhedron
    closure(const Parma_Polyhedra_Library::C_Polyhedron& A);
    
    
                  
    Parma_Polyhedra_Library::C_Polyhedron
    intersection(const Parma_Polyhedra_Library::C_Polyhedron& A,
                 const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    Parma_Polyhedra_Library::C_Polyhedron
    regular_intersection(const Parma_Polyhedra_Library::C_Polyhedron& A,
                         const Parma_Polyhedra_Library::C_Polyhedron& B);
                
                
    Parma_Polyhedra_Library::C_Polyhedron
    convex_hull(const Parma_Polyhedra_Library::C_Polyhedron& A,
                const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    Parma_Polyhedra_Library::C_Polyhedron
    minkowski_sum(const Parma_Polyhedra_Library::C_Polyhedron& A,
                  const Parma_Polyhedra_Library::C_Polyhedron& B);
                  
    Parma_Polyhedra_Library::C_Polyhedron
    minkowski_difference(const Parma_Polyhedra_Library::C_Polyhedron& A,
                         const Parma_Polyhedra_Library::C_Polyhedron& B);
                         
                         
  }
}
                   
#endif /* _ARIADNE_PPL_POLYHEDRON_H */
