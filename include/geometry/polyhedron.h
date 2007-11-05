/***************************************************************************
 *            polyhedron.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#ifndef ARIADNE_POLYHEDRON_H
#define ARIADNE_POLYHEDRON_H

#include <iosfwd>
#include <vector>

#include "../base/tribool.h"
#include "../base/iterator.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/zonotope.decl.h"

namespace Ariadne {  
  namespace Geometry {

    class basic_set_tag;
    template<class R> class PolyhedralConstraint;
    template<class R> class Polytope;
    template<class R> class PolyhedronConstraintsIterator;
    
    /*! \ingroup BasicSet
     *  \brief A polyhedron (not necessarily bounded polyhedral set) described by a system of linear inequalities.
     *
     *  The set is described as
     *  \f$ x\in\mathbb{R}^d \mid Ax\leq b , \f$
     *  where \f$A\f$ is a \f$n\times d\f$ matrix and \f$b\f$ is a vector of size \f$n\f$.
     */ 
    template<class X>
    class Polyhedron {
      typedef typename Numeric::traits<X>::number_type R;
      typedef typename Numeric::traits<X>::arithmetic_type F;
      typedef typename Numeric::traits<X>::interval_type I;
     private:
      dimension_type _dimension;
      size_type _number_of_constraints;
      array<X> _data;
     private:
      LinearAlgebra::MatrixSlice<X> _constraints();
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real numbers used to describe the polyhedron. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the polyhedron. */
      typedef Point<X> state_type;
      /*! \brief An iterator over the constraints of the Polyhedron. */
      typedef PolyhedronConstraintsIterator<X> constraints_const_iterator;
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Construct full Euclidean space of dimension \a n.
       */
      explicit Polyhedron(dimension_type n=0);
     
      /*! \brief Construct a polyhedron of dimension \a d with \a nc constraints from the data in the
       *  array beginning at \a data. The jth element of the ith constraint is stored in position i*(d+1)+j, 
       *  and the ith inhomogeneous term is stored in position i*(d+1)+d.
       */
      explicit Polyhedron(dimension_type d, size_type nc, const X* data);
            
      /*! \brief Construct the polyhedron defined by a string literal.
       */
      explicit Polyhedron(const std::string& str);
            
      /*! \brief Construct the polyhedron defined by the matrix equations \f$Ax\leq b\f$.
       */
      explicit Polyhedron(const LinearAlgebra::Matrix<X>& A, const LinearAlgebra::Vector<X>& b);
            
      /*! \brief Construct the polyhedron defined by the matrix of constraints \f$C\hat{x}\geq0\f$.
       */
      explicit Polyhedron(const LinearAlgebra::Matrix<X>& C);
            
      /*! \brief Construct from a list of points. */
      explicit Polyhedron(const PointList<X>& pts);
            
      /*! \brief Convert from a rectangle. */
      template<class XX> explicit Polyhedron(const Rectangle<XX>& rect);
            
      /*! \brief Convert from a polytope. */
      template<class XX> explicit Polyhedron(const Polytope<XX>& plyt);
            
      /*! \brief Copy constructor. */
      template<class XX> Polyhedron(const Polyhedron<XX>& original);
          
      /*! \brief Copy assignment operator. */
      template<class XX> Polyhedron<X>& operator=(const Polyhedron<XX>& original);
      //@}
      
      
      //@{
      //! \name Data access
      /*! \brief The matrix of constraints \f$C\f$ in the inequalities \f$C\left(\begin{array}{c}x\\1\end{array}\right)\geq 0\f$. */
      const LinearAlgebra::MatrixSlice<X> constraints() const;
      /*! \brief The matrix \f$A\f$ in the inequalities \f$Ax\leq b\f$. */
      LinearAlgebra::Matrix<X> A() const;
      /*! \brief The vector \f$b\f$ in the inequalities \f$Ax\leq b\f$. */
      LinearAlgebra::Vector<X> b() const;
      /*! \brief An iterator to the beginning of the constraints. */
      size_type number_of_constraints() const;
      /*! \brief An iterator to the beginning of the constraints. */
      constraints_const_iterator constraints_begin() const;
      /*! \brief An iterator to the end of the constraints. */
      constraints_const_iterator constraints_end() const;
      /*! A reference to the array of real data. */
      array<X>& data();
      /*! A constant reference to the array of real data. */
      const array<X>& data() const;
      /*! A pointer to the beginning of the array of real data. */
      X* begin();
      /*! A constant pointer to the beginning of the array of real data. */
      const X* begin() const;
      //@}


      //@{
      //! \name Modifying operations
      void new_constraint(const PolyhedralConstraint<X>& c);
      //@}


      //@{
      //! \name Geometric operations
      /*! \brief Returns the polyhedron's space dimension.
       */
      dimension_type dimension() const;
      
      /*! \brief Checks for emptyness.
       */
      tribool empty() const;
      
      /*! \brief Checks for boundedness.
       */
      tribool bounded() const;
      
      /*! \brief Tests if a point is an element of the polyhedron.
       */
      template<class XX> tribool contains(const Point<XX>& point) const;

      /*! \brief A rectangle containing the polyhedron. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief An over-approximation of the polyhedron governed by the parameter \a delta.
       *
       * WARNING: The metric error of the approximation may be larger than \a delta.
       */
      Polyhedron<X> over_approximation(const R& delta) const;
      //@}
      
      //@{
      //! \name Conversion operations
      //@}

#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality. */
      friend tribool Geometry::equal<>(const Polyhedron<X>& A, 
                                       const Polyhedron<X>& B);
      /*! \brief Tests disjointness. */
      friend tribool Geometry::disjoint<>(const Polyhedron<X>& A, 
                                          const Polyhedron<X>& B);
        
      /*! \brief Tests disjointness. */
      friend tribool Geometry::disjoint<>(const Polyhedron<X>& A, 
                                          const Rectangle<X>& B);
        
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Rectangle<X>& A, 
                                        const Polyhedron<X>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polyhedron<X>& A, 
                                        const Polyhedron<X>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polyhedron<X>& A, 
                                        const Rectangle<X>& B);
    
      /*! \brief The generators of a polyhedron. (Deprecated) */
      friend Matrix<X> Geometry::generators<>(const Polyhedron<X>& A);
    
      //@}
      

      //@{
      //! \name Geometric binary operations
      /*! \brief The intersection of two polyhedra. */
      friend Polyhedron<X> closed_intersection<>(const Polyhedron<X>& A, 
                                          const Polyhedron<X>& B);
    
      /*! \brief The closure of the intersection of the interiors of two polyhedra. */
      friend Polyhedron<X> open_intersection<>(const Polyhedron<X>& A, 
                                                  const Polyhedron<X>& B);
    
      //@}
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
     private:
      static void _instantiate_geometry_operators();
    };
    
    
    
 
    
    /*! \brief A linear inequality constraint. */
    template<class X>
    class PolyhedralConstraint
    {
      friend class Polyhedron<X>;
      friend class PolyhedronConstraintsIterator<X>;
     public:
      /*! \brief The dimension of the constraint. */
      dimension_type dimension() const;
      /*! \brief Tests if the constraint is satisfied by a point. */
      template<class XV> tribool satisfied_by(const Point<XV>& pt) const;
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      PolyhedralConstraint(const dimension_type d, const X* a);
      dimension_type _d; const X* _a;
    };
    
    
    
    template<class X> 
    tribool empty(const Polyhedron<X>& plhd);

    template<class X> 
    tribool bounded(const Polyhedron<X>& plhd);

    template<class X> 
    Rectangle<typename Polyhedron<X>::real_type> bounding_box(const Polyhedron<X>& plhd);


    template<class X> 
    tribool equal(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2);

   
    template<class X> 
    tribool disjoint(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2);
    
    template<class X> 
    tribool disjoint(const Polyhedron<X>& plhd1, const Rectangle<X>& plhd2);
    
    template<class X> 
    tribool disjoint(const Rectangle<X>& plhd1, const Polyhedron<X>& plhd2);
    
    
    template<class X> 
    tribool subset(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2);

    template<class X> 
    tribool subset(const Polyhedron<X>& plhd1, const Rectangle<X>& plhd2);


    template<class X1,class X2> 
    tribool subset(const Rectangle<X1>& r, const Polyhedron<X2>& plhd);

    template<class X1,class X2> 
    tribool subset(const Zonotope<X1>& z, const Polyhedron<X2>& plhd);

    template<class X1,class X2> 
    tribool subset(const Polytope<X1>& pltp, const Polyhedron<X2>& plhd);



    template<class X> 
    Polyhedron<X> 
    open_intersection(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2);

    template<class X> 
    Polyhedron<X> 
    closed_intersection(const Polyhedron<X>& plhd1, const Polyhedron<X>& plhd2) ;

    template<class X> 
    Polyhedron<X> 
    closed_intersection(const Polyhedron<X>& plhd, const Rectangle<X>& r) ;

    template<class X> 
    Polyhedron<X> 
    closed_intersection(const Rectangle<X>& r, const Polyhedron<X>& plhd) ;

    
    template<class X> 
    Polyhedron<X> 
    polyhedron(const Rectangle<X>& A) ;

    template<class X> 
    Polyhedron<typename Numeric::traits<X>::arithmetic_type> 
    polyhedron(const Polytope<X>& A) ;



    
    template<class X>
    std::ostream& operator<<(std::ostream& os, const PolyhedralConstraint<X>& c);
    
    template<class X>
    std::ostream& operator<<(std::ostream& os, const Polyhedron<X>& p);
    
    template<class X>
    std::istream& operator>>(std::istream& os, Polyhedron<X>& p);



  }
}

#include "polyhedron.inline.h"

#endif /* ARIADNE_POLYHEDRON_H */
