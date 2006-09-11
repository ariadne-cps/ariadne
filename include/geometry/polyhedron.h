/***************************************************************************
 *            polyhedron.h
 *
 *  Thu Jan 27 10:26:36 2005
 *  Copyright  2005  Alberto Casagrande
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

/*! \file polyhedron.h
 *  \brief Polyhedra.
 */
 
#ifndef _ARIADNE_POLYHEDRON_H
#define _ARIADNE_POLYHEDRON_H

#include <iosfwd>
#include <vector>

#include "../declarations.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"

#include "../geometry/ppl_polyhedron.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"

namespace Ariadne {  
  namespace Geometry {


    /*! \ingroup BasicSet
     *  \brief A polyhedron (bounded polyhedral set) described by its vertices.
     */ 
    template<typename R>
    class Polyhedron {
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
       /*! \brief Construct an empty polyhedron in dimension \a n. */
      Polyhedron(dimension_type n=0);
     
     /*! \brief Construct from a matrix whose columns give the vertices. */
      explicit Polyhedron(const LinearAlgebra::Matrix<R>& A);
     
      /*! \brief Convert from a list of vertices. */
      Polyhedron(const PointList<R>& v);
     
      /*! \brief Convert from a rectangle. */
      Polyhedron<R>(const Rectangle<R>& rect);
            
      /*! \brief Copy constructor. 
       */
      Polyhedron<R>(const Polyhedron<R>& original);
          
      /*! \brief Copy assignment operator. 
       *
       * \param original is the original polyhedron.
       * \return A reference to the current object.
       */
      Polyhedron<R>& operator=(const Polyhedron<R>& original);
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
      /*! \brief Returns the polyhedron's space dimension.
       */
      dimension_type dimension() const;
      
      /*! \brief Checks for emptyness.
       */
      bool empty() const;
      
      /*! \brief Checks for emptyness of the interior.
       */
      bool empty_interior() const;
   
#ifdef DOXYGEN
      /*! \brief The vertices of the polyhedron. */
      state_list_type vertices() const; 
#endif

      /*! \brief Tests if a point is an element of the polyhedron.
       */
      bool contains(const state_type& point) const;

      /*! \brief Tests if a point is an element of the interior of the polyhedron.
       */
      bool interior_contains(const state_type& point) const;
           
      /*! \brief A rectangle containing the polyhedron. */
      Rectangle<R> bounding_box() const;
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality. */
      friend bool Geometry::equal<>(const Polyhedron<R>& A, 
                                    const Polyhedron<R>& B);
      /*! \brief Tests disjointness. */
      friend bool Geometry::disjoint<>(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
    
      /*! \brief Tests intersection of interiors. */
      friend bool Geometry::interiors_intersect<>(const Polyhedron<R>& A, 
                                                  const Polyhedron<R>& B);
      
      /*! \brief Tests inclusion of \a A in interior of \a B. */
      friend bool Geometry::inner_subset<>(const Polyhedron<R>& A, 
                                           const Polyhedron<R>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend bool Geometry::subset<>(const Polyhedron<R>& A, 
                                     const Polyhedron<R>& B);
    
      //@}
      

      //@{
      //! \name Geometric binary operations
      /*! \brief The intersection of two polyhedra. */
      friend Polyhedron<R> intersection<>(const Polyhedron<R>& A, 
                                          const Polyhedron<R>& B);
    
      /*! \brief The closure of the intersection of the interiors of two polyhedra. */
      friend Polyhedron<R> regular_intersection<>(const Polyhedron<R>& A, 
                                                  const Polyhedron<R>& B);
    
      /*! \brief The convex hull of two polyhedra. */
      friend Polyhedron<R> convex_hull<>(const Polyhedron<R>& A, 
                                         const Polyhedron<R>& B);
      
      /*! \brief The Minkowski (pointwise) sum of two polyhedra. */
      friend Polyhedron<R> minkowski_sum<>(const Polyhedron<R>& A, 
                                           const Polyhedron<R>& B);
   
      /*! \brief The Minkowski (pointwise) difference of two polyhedra. */
      friend Polyhedron<R> minkowski_difference<>(const Polyhedron<R>& A, 
                                                  const Polyhedron<R>& B);
      //@}
 #else
      static bool equal(const Polyhedron<R>& A, const Polyhedron<R>& B);
      static bool disjoint(const Polyhedron<R>& A, const Polyhedron<R>& B);
      static bool interiors_intersect(const Polyhedron<R>& A, const Polyhedron<R>& B);
      static bool subset(const Polyhedron<R>& A, const Polyhedron<R>& B);
      static bool inner_subset(const Polyhedron<R>& A, const Polyhedron<R>& B);
      
      static Polyhedron<R> convex_hull(const Polyhedron<R>& A, const Polyhedron<R>& B);
#endif
     private:
      LinearAlgebra::Matrix<R> _vertices;
    };
   
    template<typename R> std::ostream& operator<<(std::ostream& os, const Polyhedron<R>& p);
    template<typename R> std::istream& operator>>(std::istream& os, Polyhedron<R>& p);
   
    
    template<typename R> inline 
    bool equal(const Polyhedron<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::equal(A,B);
    }
   
    template<typename R> inline 
    bool disjoint(const Polyhedron<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::disjoint(A,B);
    }
    
    template<typename R> inline 
    bool 
    disjoint(const Polyhedron<R>& A, const Rectangle<R>& B) {
      return Polyhedron<R>::disjoint(A,Polyhedron<R>(B));
    }
    
    template<typename R> inline 
    bool 
    disjoint(const Rectangle<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::disjoint(Polyhedron<R>(A),B);
    }
    
    
    template<typename R> inline 
    bool
    interiors_intersect(const Polyhedron<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::interiors_intersect(A,B);
    }
    
    template<typename R> inline 
    bool
    interiors_intersect(const Polyhedron<R>& A, const Rectangle<R>& B) {
      return Polyhedron<R>::interiors_intersect(A,Polyhedron<R>(B));
    }
    
    template<typename R> inline 
    bool
    interiors_intersect(const Rectangle<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::interiors_intersect(Polyhedron<R>(A),B);
    }
    
   
    template<typename R> inline 
    bool 
    inner_subset(const Polyhedron<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::inner_subset(A,B);
    }
    
    template<typename R> inline 
    bool 
    inner_subset(const Polyhedron<R>& A, const Rectangle<R>& B) {
      return Polyhedron<R>::inner_subset(A,Polyhedron<R>(B));
    }

    template<typename R> inline 
    bool 
    inner_subset(const Rectangle<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::inner_subset(Polyhedron<R>(A),B);
    }
      
      
    template<typename R> inline 
    bool 
    subset(const Polyhedron<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::subset(A,B);
    }

    template<typename R> inline 
    bool 
    subset(const Polyhedron<R>& A, const Rectangle<R>& B) {
      return Polyhedron<R>::subset(A,Polyhedron<R>(B));
    }

    template<typename R> inline 
    bool 
    subset(const Rectangle<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::subset(Polyhedron<R>(A),B);
    }


    template<typename R> inline 
    Polyhedron<R> 
    convex_hull(const Polyhedron<R>& A, const Polyhedron<R>& B) {
      return Polyhedron<R>::convex_hull(A,B);
    }
  
    
    
    
    /*! \ingroup BasicSet
     *  \brief A polytope (not necessarily bounded polyhedral set) described by a system of linear inequalities.
     */ 
    template<typename R>
    class Polytope {
     public:
      /*! \brief The type of denotable real numbers used to describe the polytope. */
      typedef Rational real_type;
      /*! \brief The type of denotable point contained by the polytope. */
      typedef Point<R> state_type;
      /*! \brief The type of vector. */
      typedef Ariadne::LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix. */
      typedef Ariadne::LinearAlgebra::Matrix<R> matrix_type;
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Construct full Euclidean space of dimension \a n.
       */
      explicit Polytope<R>(dimension_type n=0);
     
      /*! \brief Construct the polytope defined by the matrix equations \f$Ax\leq b\f$.
       */
      explicit Polytope<R>(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b);
            
      /*! \brief Convert from a rectangle. */
      Polytope<R>(const Rectangle<R>& rect);
            
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
      
      
      //@{
      //! \name Conversion operations
      /*! \brief Convert to a Parma Polyhedral Library closed polyhedron. */
      operator Parma_Polyhedra_Library::C_Polyhedron () const;
      //@}

      //@{
      //! \name Data access
      /*! \brief The matrix \f$A\f$ in the inequalities \f$Ax\leq b\f$. */
      LinearAlgebra::Matrix<R> A() const { return _A; }
      /*! \brief The vector \f$b\f$ in the inequalities \f$Ax\leq b\f$. */
      LinearAlgebra::Vector<R> b() const { return _b; }
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
   
      /*! \brief The vertices of the polytope. */
      PointList<Rational>  vertices() const; 
           
      /*! \brief Tests if a point is an element of the polytope.
       */
      bool contains(const state_type& point) const;

      /*! \brief Tests if a point is an element of the interior of the polytope.
       */
      bool interior_contains(const state_type& point) const;
           
      /*! \brief A rectangle containing the polytope. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief An over-approximation of the polytope governed by the parameter \a delta.
       *
       * WARNING: The metric error of the approximation may be larger than \a delta.
       */
      Polytope<R> over_approximation(const R& delta) const;
      //@}
      

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
      static bool equal(const Polytope& A, const Polytope& B);
      static bool disjoint(const Polytope& A, const Polytope& B);
      static bool interiors_intersect(const Polytope& A, const Polytope& B);
      static bool subset(const Polytope& A, const Polytope& B);
      static bool inner_subset(const Polytope& A, const Polytope& B);
      
      static Polytope<R> intersection(const Polytope& A, const Polytope& B);
      static Polytope<R> regular_intersection(const Polytope& A, const Polytope& B);
#endif
     private:
       LinearAlgebra::Matrix<R> _A;
       LinearAlgebra::Vector<R> _b;
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
    regular_intersection(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::regular_intersection(A,B);
    }

    template<typename R> inline 
    Polytope<R> 
    intersection(const Polytope<R>& A, const Polytope<R>& B) {
      return Polytope<R>::intersection(A,B);
    }



  }
}

#endif /* _ARIADNE_POLYHEDRON_H */
