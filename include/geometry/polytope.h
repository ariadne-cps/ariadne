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
 
#ifndef ARIADNE_POLYTOPE_H
#define ARIADNE_POLYTOPE_H

#include <iosfwd>
#include <vector>

#include "base/tribool.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/rectangle.h"

namespace Ariadne {  
  namespace Geometry {

    class basic_set_tag;
    template<class R> class Polyhedron;
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
    template<class X>
    class Polytope {
      typedef typename Numeric::traits<X>::number_type R;
      typedef typename Numeric::traits<X>::arithmetic_type A;
      typedef typename Numeric::traits<X>::interval_type I;
     private:    
      dimension_type _dimension;
      size_type _number_of_vertices;
      array<X> _data;
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
       /*! \brief The type of denotable real numbers used to describe the convex hull. */
      typedef R value_type;
       /*! \brief The type of real numbers used to describe the state space. */
      typedef typename Numeric::traits<X>::number_type real_type;
      /*! \brief The type of denotable point contained by the convex hull. */
      typedef Point<X> state_type;
      /*! \brief An iterator through the vertices of the polytope. */
      typedef PolytopeVerticesIterator<X> vertices_const_iterator;
     public:
      /*! \brief Construct an empty polytope in dimension \a n. */
      Polytope(dimension_type d=0);
     
      /*! \brief Construct the polyhedron defined by a string literal.
       */
      explicit Polytope<X>(const std::string& str);
            
      /*! \brief Construct a polytope of dimension \a d with \a nv constraints from the data in the
       *  array beginning at \a data. The ith element of the jth vertex is stored in position j*(d+1)+i. 
       *  The j(d+1)+dth data element is padded with a 1.
       */
      explicit Polytope<X>(dimension_type d, size_type nv, const X* data);

      /*! \brief Construct from a matrix whose columns give the vertices. */
      explicit Polytope(const LinearAlgebra::Matrix<X>& A);
     
      /*! \brief Construct from a list of vertices. */
      explicit Polytope(const PointList<X>& v);
     
      /*! \brief Construct from a rectangle. */
      explicit Polytope(const Rectangle<X>& rect);
            
      /*! \brief Construct from a polyhedron. */
      template<class XX> explicit Polytope(const Polyhedron<XX>& p);
            
      /*! \brief Copy constructor. */
      template<class XX> Polytope(const Polytope<XX>& original);
          
      /*! \brief Copy assignment operator. */
      template<class XX> Polytope<X>& operator=(const Polytope<XX>& original);
      //@}
 
      /*! \brief The matrix of generators. */
      const LinearAlgebra::Matrix<X> generators() const;
      /*! \brief The number of vertices of the convex set. */
      size_type number_of_vertices() const;
      /*! \brief The vertices of the convex set. */
      PointList<X> vertices() const;
      /*! \brief The \a i th vertex. */
      Point<X> vertex(const size_type& i) const;
      /*! \brief An iterator to the beginning of the vertices. */
      vertices_const_iterator vertices_begin() const;
      /*! \brief An iterator to the end of the vertices. */
      vertices_const_iterator vertices_end() const;
      
      const array<X>& data() const { 
        return this->_data; 
      }
      void resize(dimension_type d, size_type nv) { 
        this->_dimension=d; 
        this->_number_of_vertices=nv;
        this->_data.resize((d+1)*nv);
      }
      
      LinearAlgebra::MatrixSlice<X> _generators_() { 
        return LinearAlgebra::MatrixSlice<X>(this->_dimension+1,this->_number_of_vertices,this->_data.begin(),1,this->_dimension+1); 
      }
      //@{
      //! \name Conversion operations
      //@}


      //@{
      //! \name Geometric operations
      /*! \brief Returns the polytope's space dimension.
       */
      dimension_type dimension() const;
      
      /*! \brief Checks for emptyness. */
      tribool empty() const;
         
      /*! \brief Checks for boundedness. */
      tribool bounded() const;
      
#ifdef DOXYGEN
      /*! \brief The vertices of the polytope. */
      PointList<X> vertices() const; 
      /*! \brief An iterator to the vertices of the polytope. */
      vertices_const_iterator vertices_begin() const; 
      /*! \brief An iterator to the vertices of the polytope. */
      vertices_const_iterator vertices_end() const; 
#endif

      /*! \brief Tests if a point is an element of the polytope.
       */
      tribool contains(const Point<X>& point) const;
           
      /*! \brief A rectangle containing the polytope. */
      Rectangle<R> bounding_box() const;
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality. */
      friend tribool Geometry::equal<>(const Polytope<X>& A, 
                                       const Polytope<X>& B);
      /*! \brief Tests disjointness. */
      friend tribool Geometry::disjoint<>(const Polytope<X>& A, 
                                          const Polytope<X>& B);
        
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polytope<X>& A, 
                                        const Rectangle<X>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polytope<X>& A, 
                                        const Polyhedron<X>& B);
      //@}
      

      //@{
      //! \name Geometric binary operations
      /*! \brief The intersection of two polyhedra. */
      friend Polytope<X> closed_intersection<>(const Polytope<X>& A, 
                                               const Polytope<X>& B);
    
      /*! \brief The closure of the intersection of the interiors of two polyhedra. */
      friend Polytope<X> open_intersection<>(const Polytope<X>& A, 
                                             const Polytope<X>& B);
    
      /*! \brief The convex hull of two polyhedra. */
      friend Polytope<X> convex_hull<>(const Polytope<X>& A, 
                                         const Polytope<X>& B);
      
      /*! \brief The Minkowski (pointwise) sum of two polyhedra. */
      friend Polytope<X> minkowski_sum<>(const Polytope<X>& A, 
                                         const Polytope<X>& B);
   
      /*! \brief The Minkowski (pointwise) difference of two polyhedra. */
      friend Polytope<X> minkowski_difference<>(const Polytope<X>& A, 
                                                  const Polytope<X>& B);
      //@}
 #else
     private:
      static void _instantiate_geometry_operators();
     public:
#endif

      //@{ 
      //! \name Input/output operations
      /*! \brief The name of the class. */
      static std::string name();
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
    };
   
    template<class X>
    std::ostream& operator<<(std::ostream& os, const Polytope<X>& p);
    
    template<class X>
    std::istream& operator>>(std::istream& os, Polytope<X>& p);
    
    

    template<class X>  
    Rectangle<typename Polytope<X>::real_type> bounding_box(const Polytope<X>& pltp);

    
    template<class X> 
    tribool equal(const Polytope<X>& pltp1, const Polytope<X>& pltp2);

   
    template<class X>  
    tribool contains(const Polytope<X>& pltp, const Point<X>& pt);
    
    template<class X>  
    tribool disjoint(const Polytope<X>& pltp1, const Polytope<X>& pltp2);
    
    template<class X>
    tribool disjoint(const Polytope<X>& pltp, const Rectangle<X>& r);
    
    template<class X> 
    tribool disjoint(const Rectangle<X>& r, const Polytope<X>& pltp);
    
      
    template<class X> 
    tribool subset(const Polytope<X>& pltp2, const Polytope<X>& pltp2);

    template<class X> 
    tribool subset(const Polytope<X>& pltp2, const Rectangle<X>& r);

    template<class X> 
    tribool subset(const Rectangle<X>& r, const Polytope<X>& pltp);


    template<class X> 
    Polytope<X> convex_hull(const Polytope<X>& pltp1, const Polytope<X>& pltp2);
  
    
    template<class X> 
    Polytope<X> polytope(const Rectangle<X>& r);
  
    template<class X> 
    Polytope<typename Numeric::traits<X>::arithmetic_type> 
    polytope(const Polyhedron<X>& plhd);
  

  }
}

#include "polytope.inline.h"

#endif /* ARIADNE_POLYTOPE_H */
