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

#include "../base/tribool.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"

namespace Ariadne {  
  namespace Geometry {

    template<class R> class PolytopeVerticesIterator;
      
    /*! \ingroup BasicSet
     *  \brief A polytope (bounded polyhedral set) described by its vertices.
     *
     *  The vertices are stored in a matrix of \em generators. The columns of
     *  this matrix are the position vectors of the vertices, padded by 1's.
     *  Hence the generator matrix has size \f$(d+1)\times nv\f$, where nv is
     *  the number of vertices.
     *
     *  The polytope is described as 
     *     \f$ \{ x \in \mathbb{R}^d \mid \exists s\in\mathbb{R}^n,\ x=Vs, s\geq 0,\ 1^Ts=1 \} . \f$
     *
     *  \internal It would be preferable to store the generator matrix in 
     *  column major format.
     */ 
    template<class R>
    class Polytope {
      typedef typename Numeric::traits<R>::arithmetic_type F;
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
      /*! \brief The type of a list of points. */
      typedef PolytopeVerticesIterator<R> vertices_iterator;
      typedef PolytopeVerticesIterator<R> vertices_const_iterator;
     public:
      /*! \brief Construct an empty polytope in dimension \a n. */
      Polytope(dimension_type d=0);
     
      /*! \brief Construct a polytope of dimension \a d with \a nv constraints from the data in the
       *  array beginning at \a data. The ith element of the jth vertex is stored in position j*(d+1)+i. 
       *  The j(d+1)+dth data element is padded with a 1.
       */
      explicit Polytope<R>(dimension_type d, size_type nv, const R* data);

     /*! \brief Construct from a matrix whose columns give the vertices. */
      explicit Polytope(const LinearAlgebra::Matrix<R>& A);
     
      /*! \brief Construct from a list of vertices. */
      explicit Polytope(const PointList<R>& v);
     
      /*! \brief Construct from a rectangle. */
      explicit Polytope(const Rectangle<R>& rect);
            
      /*! \brief Construct from a polyhedron. */
      explicit Polytope(const Polyhedron<R>& p);
            
      /*! \brief Copy constructor. 
       */
      template<class OR> inline Polytope(const Polytope<OR>& original);
          
      /*! \brief Copy assignment operator. 
       *
       * \param original is the original polytope.
       * \return A reference to the current object.
       */
      Polytope<R>& operator=(const Polytope<R>& original);
      //@}
 
      /*! \brief The matrix of generators. */
      const LinearAlgebra::Matrix<R>& generators() const;
      /*! \brief The number of vertices of the convex set. */
      size_type number_of_vertices() const;
      /*! \brief The vertices of the convex set. */
      PointList<R> vertices() const;
      /*! \brief The \a i th vertex. */
      Point<R> vertex(const size_type& i) const;
      /*! \brief An iterator to the beginning of the vertices. */
      vertices_const_iterator vertices_begin() const;
      /*! \brief An iterator to the end of the vertices. */
      vertices_const_iterator vertices_end() const;
      
      //@{
      //! \name Conversion operations
      //@}


      //@{
      //! \name Geometric operations
      /*! \brief Returns the polytope's space dimension.
       */
      dimension_type dimension() const;
      
      /*! \brief Checks for emptyness.
       */
      tribool empty() const;
         
#ifdef DOXYGEN
      /*! \brief The vertices of the polytope. */
      state_list_type vertices() const; 
      /*! \brief An iterator to the vertices of the polytope. */
      vertices_const_iterator vertices_begin() const; 
      /*! \brief An iterator to the vertices of the polytope. */
      vertices_const_iterator vertices_end() const; 
#endif

      /*! \brief Tests if a point is an element of the polytope.
       */
      tribool contains(const state_type& point) const;
           
      /*! \brief A rectangle containing the polytope. */
      Rectangle<R> bounding_box() const;
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality. */
      friend tribool Geometry::equal<>(const Polytope<R>& A, 
                                       const Polytope<R>& B);
      /*! \brief Tests disjointness. */
      friend tribool Geometry::disjoint<>(const Polytope<R>& A, 
                                          const Polytope<R>& B);
        
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polytope<R>& A, 
                                        const Rectangle<R>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polytope<R>& A, 
                                        const Polyhedron<R>& B);
      //@}
      

      //@{
      //! \name Geometric binary operations
      /*! \brief The intersection of two polyhedra. */
      friend Polytope<R> closed_intersection<>(const Polytope<R>& A, 
                                               const Polytope<R>& B);
    
      /*! \brief The closure of the intersection of the interiors of two polyhedra. */
      friend Polytope<R> open_intersection<>(const Polytope<R>& A, 
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
     private:
      static void _instantiate_geometry_operators();
     public:
  #endif

      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     private:    
      LinearAlgebra::Matrix<R> _generators;
    };
   
    template<class R>
    class Polytope< Interval<R> >
    {
      typedef Interval<R> I;
     public:
      template<class Rl> Polytope(const PointList<Rl>& v);
      template<class Rl> explicit Polytope(const Rectangle<Rl>& r);
      dimension_type dimension() const { return _vertices.dimension(); }
      dimension_type number_of_vertices() const { return _vertices.size(); }
      const PointList< Interval<R> >& vertices() const;
     private:
      PointList<I> _vertices;
    };

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Polytope<R>& p) {
      return p.write(os); }
    template<class R> inline
    std::istream& operator>>(std::istream& os, Polytope<R>& p) {
      return p.read(os); }
   

    template<class R>
    class PolytopeVerticesIterator
      : public boost::iterator_facade<PolytopeVerticesIterator<R>,
                                      Point<R>,
                                      boost::forward_traversal_tag,
                                      Point<R> const&,
                                      Point<R> const*
                                     >
    {
     public:
      PolytopeVerticesIterator(const Polytope<R>& plytp, const size_type& j) 
        : _p(&plytp), _j(j) { }
      bool equal(const PolytopeVerticesIterator<R>& other) const {
        return this->_j==other._j && this->_p ==other._p; }
      void increment() {
        ++this->_j; }
      const Point<R>& dereference() const {
        this->_v=this->_p->vertex(this->_j); return this->_v; }
     private:
      const Polytope<R>* _p; size_type _j; mutable Point<R> _v;
    };
    
    template<class R> 
    tribool equal(const Polytope<R>& A, const Polytope<R>& B);
   
    template<class R>  
    tribool disjoint(const Polytope<R>& A, const Polytope<R>& B);
    
    template<class R>
    tribool disjoint(const Polytope<R>& A, const Rectangle<R>& B);
    
    template<class R> 
    tribool disjoint(const Rectangle<R>& A, const Polytope<R>& B);
    
      
      
    template<class R> 
    tribool subset(const Polytope<R>& A, const Polytope<R>& B);

    template<class R> 
    tribool subset(const Polytope<R>& A, const Rectangle<R>& B);

    template<class R> 
    tribool subset(const Rectangle<R>& A, const Polytope<R>& B);


    template<class R> 
    Polytope<R> 
    convex_hull(const Polytope<R>& A, const Polytope<R>& B);
  
    
    
    template<class R> template<class Rl> inline 
    Polytope<R>::Polytope(const Polytope<Rl>& original)
      : _generators(original.generators()) { }


  }
}

#endif /* _ARIADNE_POLYTOPE_H */
