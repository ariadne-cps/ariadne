/***************************************************************************
 *            polytope.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it   Pieter.Collins@cwi.nl
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

/*! \file polytope.h
 *  \brief Polytopes.
 */
 
#ifndef _ARIADNE_POLYTOPE_H
#define _ARIADNE_POLYTOPE_H

#include <iosfwd>
#include <vector>

#include "../declarations.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/ppl_polyhedron.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"

namespace Ariadne {  
  namespace Geometry {


    /*! \ingroup BasicSet
     *  \brief A polytope (bounded polyhedral set) described by its vertices.
     */ 
    template<typename R>
    class Polytope {
      public:
       /*! \brief The type of denotable real numbers used to describe the convex hull. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the convex hull. */
      typedef Point<R> state_type;
      /*! \brief The type of vector. */
      typedef Ariadne::LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix. */
      typedef Ariadne::LinearAlgebra::Matrix<R> matrix_type;
      /*! \brief The type of a list of points. */
      typedef PointList<R>  state_list_type;
     public:
       /*! \brief Construct an empty polytope in dimension \a n. */
      Polytope(dimension_type d=0);
     
     /*! \brief Construct from a matrix whose columns give the vertices. */
      explicit Polytope(const LinearAlgebra::Matrix<R>& A);
     
      /*! \brief Construct from a list of vertices. */
      explicit Polytope(const PointList<R>& v);
     
      /*! \brief Construct from a rectangle. */
      explicit Polytope<R>(const Rectangle<R>& rect);
            
      /*! \brief Copy constructor. 
       */
      Polytope<R>(const Polytope<R>& original);
          
      /*! \brief Copy assignment operator. 
       *
       * \param original is the original polytope.
       * \return A reference to the current object.
       */
      Polytope<R>& operator=(const Polytope<R>& original);
      //@}
 
      /*! \brief The number of vertices of the convex set. */
      size_type number_of_vertices() const;

      /*! \brief The vertices of the convex set. */
      PointList<R> vertices() const;
      
      //@{
      //! \name Conversion operations
      /*! \brief Convert to a Parma Polyhedral Library closed polyhedron. */
      operator Parma_Polyhedra_Library::C_Polyhedron () const;
      //@}


      //@{
      //! \name Geometric operations
      /*! \brief Returns the polytope's space dimension.
       */
      dimension_type dimension() const;
      
      /*! \brief Checks for emptyness.
       */
      bool empty() const;
      
      /*! \brief Checks for emptyness of the interior.
       */
      bool empty_interior() const;
   
#ifdef DOXYGEN
      /*! \brief The vertices of the polytope. */
      state_list_type vertices() const; 
#endif

      /*! \brief Tests if a point is an element of the polytope.
       */
      bool contains(const state_type& point) const;

      /*! \brief Tests if a point is an element of the interior of the polytope.
       */
      bool interior_contains(const state_type& point) const;
           
      /*! \brief A rectangle containing the polytope. */
      Rectangle<R> bounding_box() const;
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality. */
      friend bool Geometry::equal<>(const Polytope<R>& A, 
                                    const Polytope<R>& B);
      /*! \brief Tests disjointness. */
      friend bool Geometry::disjoint<>(const Polytope<R>& A, 
                                       const Polytope<R>& B);
    
      /*! \brief Tests intersection of interiors. */
      friend bool Geometry::interiors_intersect<>(const Polytope<R>& A, 
                                                  const Polytope<R>& B);
      
      /*! \brief Tests inclusion of \a A in interior of \a B. */
      friend bool Geometry::inner_subset<>(const Polytope<R>& A, 
                                           const Polytope<R>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend bool Geometry::subset<>(const Polytope<R>& A, 
                                     const Polytope<R>& B);
    
      //@}
      

      //@{
      //! \name Geometric binary operations
      /*! \brief The intersection of two polyhedra. */
      friend Polytope<R> intersection<>(const Polytope<R>& A, 
                                          const Polytope<R>& B);
    
      /*! \brief The closure of the intersection of the interiors of two polyhedra. */
      friend Polytope<R> regular_intersection<>(const Polytope<R>& A, 
                                                  const Polytope<R>& B);
    
      /*! \brief The convex hull of two polyhedra. */
      friend Polytope<R> convex_hull<>(const Polytope<R>& A, 
                                         const Polytope<R>& B);
      
      /*! \brief The Minkowski (pointwise) sum of two polyhedra. */
      friend Polytope<R> minkowski_sum<>(const Polytope<R>& A, 
                                           const Polytope<R>& B);
   
      /*! \brief The Minkowski (pointwise) difference of two polyhedra. */
      friend Polytope<R> minkowski_difference<>(const Polytope<R>& A, 
                                                  const Polytope<R>& B);
      //@}
 #else
      static bool equal(const Polytope<R>& A, const Polytope<R>& B);
      static bool disjoint(const Polytope<R>& A, const Polytope<R>& B);
      static bool interiors_intersect(const Polytope<R>& A, const Polytope<R>& B);
      static bool subset(const Polytope<R>& A, const Polytope<R>& B);
      static bool inner_subset(const Polytope<R>& A, const Polytope<R>& B);
      
      static Polytope<R> convex_hull(const Polytope<R>& A, const Polytope<R>& B);
#endif
     private:    
      LinearAlgebra::Matrix<R> _vertices;
    };
   
    template<typename R> std::ostream& operator<<(std::ostream& os, const Polytope<R>& p);
    template<typename R> std::istream& operator>>(std::istream& os, Polytope<R>& p);
   
    
    template<typename R> inline 
    bool equal(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::equal(A,B);
    }
   
    template<typename R> inline 
    bool disjoint(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::disjoint(A,B);
    }
    
    template<typename R> inline 
    bool 
    disjoint(const Polytope<R>& A, const Rectangle<R>& B) {
      return Polytope<R>::disjoint(A,Polytope<R>(B));
    }
    
    template<typename R> inline 
    bool 
    disjoint(const Rectangle<R>& A, const Polytope<R>& B) {
      return Polytope<R>::disjoint(Polytope<R>(A),B);
    }
    
    
    template<typename R> inline 
    bool
    interiors_intersect(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::interiors_intersect(A,B);
    }
    
    template<typename R> inline 
    bool
    interiors_intersect(const Polytope<R>& A, const Rectangle<R>& B) {
      return Polytope<R>::interiors_intersect(A,Polytope<R>(B));
    }
    
    template<typename R> inline 
    bool
    interiors_intersect(const Rectangle<R>& A, const Polytope<R>& B) {
      return Polytope<R>::interiors_intersect(Polytope<R>(A),B);
    }
    
   
    template<typename R> inline 
    bool 
    inner_subset(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::inner_subset(A,B);
    }
    
    template<typename R> inline 
    bool 
    inner_subset(const Polytope<R>& A, const Rectangle<R>& B) {
      return Polytope<R>::inner_subset(A,Polytope<R>(B));
    }

    template<typename R> inline 
    bool 
    inner_subset(const Rectangle<R>& A, const Polytope<R>& B) {
      return Polytope<R>::inner_subset(Polytope<R>(A),B);
    }
      
      
    template<typename R> inline 
    bool 
    subset(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::subset(A,B);
    }

    template<typename R> inline 
    bool 
    subset(const Polytope<R>& A, const Rectangle<R>& B) {
      return Polytope<R>::subset(A,Polytope<R>(B));
    }

    template<typename R> inline 
    bool 
    subset(const Rectangle<R>& A, const Polytope<R>& B) {
      return Polytope<R>::subset(Polytope<R>(A),B);
    }


    template<typename R> inline 
    Polytope<R> 
    convex_hull(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::convex_hull(A,B);
    }
  
    
    
    

  }
}

#endif /* _ARIADNE_POLYTOPE_H */
