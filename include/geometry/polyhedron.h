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

namespace Ariadne {  
  namespace Geometry {

    template<class R> class Constraint;
    template<class R> class Polytope;
    template<class R> class PolyhedronConstraintsIterator;
    
    /*! \ingroup BasicSet
     *  \brief A polyhedron (not necessarily bounded polyhedral set) described by a system of linear inequalities.
     *
     *  The set is described as
     *  \f$ x\in\mathbb{R}^d \mid Ax\leq b , \f$
     *  where \f$A\f$ is a \f$n\times d\f$ matrix and \f$b\f$ is a vector of size \f$n\f$.
     */ 
    template<class R>
    class Polyhedron {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     private:
      dimension_type _dimension;
      size_type _number_of_constraints;
      array<R> _data;
     private:
      LinearAlgebra::MatrixSlice<R> _constraints();
     public:
      /*! \brief The type of denotable real numbers used to describe the polyhedron. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the polyhedron. */
      typedef Point<R> state_type;
      /*! \brief An iterator over the constraints of the Polyhedron. */
      typedef PolyhedronConstraintsIterator<R> constraints_const_iterator;
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Construct full Euclidean space of dimension \a n.
       */
      explicit Polyhedron<R>(dimension_type n=0);
     
      /*! \brief Construct a polyhedron of dimension \a d with \a nc constraints from the data in the
       *  array beginning at \a data. The jth element of the ith constraint is stored in position i*(d+1)+j, 
       *  and the ith inhomogeneous term is stored in position i*(d+1)+d.
       */
      explicit Polyhedron<R>(dimension_type d, size_type nc, const R* data);
            
      /*! \brief Construct the polyhedron defined by the matrix equations \f$Ax\leq b\f$.
       */
      explicit Polyhedron<R>(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b);
            
      /*! \brief Construct the polyhedron defined by the matrix of constraints \f$C\hat{x}\geq0\f$.
       */
      explicit Polyhedron<R>(const LinearAlgebra::Matrix<R>& C);
            
      /*! \brief Construct from a list of points. */
      explicit Polyhedron<R>(const PointList<R>& pts);
            
      /*! \brief Convert from a rectangle. */
      explicit Polyhedron<R>(const Rectangle<R>& rect);
            
      /*! \brief Convert from a polytope. */
      template<class Rl> explicit Polyhedron<R>(const Polytope<Rl>& plyt);
            
      /*! \brief Copy constructor. */
      template<class Rl> Polyhedron(const Polyhedron<Rl>& original);
          
      /*! \brief Copy assignment operator. */
      Polyhedron<R>& operator=(const Polyhedron<R>& original);
      //@}
      
      
      //@{
      //! \name Data access
      /*! \brief The matrix of constraints \f$C\f$ in the inequalities \f$C\left(\begin{array}{c}x\\1\end{array}\right)\geq 0\f$. */
      const LinearAlgebra::MatrixSlice<R> constraints() const;
      /*! \brief The matrix \f$A\f$ in the inequalities \f$Ax\leq b\f$. */
      LinearAlgebra::Matrix<R> A() const;
      /*! \brief The vector \f$b\f$ in the inequalities \f$Ax\leq b\f$. */
      LinearAlgebra::Vector<R> b() const;
      /*! \brief An iterator to the beginning of the constraints. */
      size_type number_of_constraints() const;
      /*! \brief An iterator to the beginning of the constraints. */
      constraints_const_iterator constraints_begin() const;
      /*! \brief An iterator to the end of the constraints. */
      constraints_const_iterator constraints_end() const;
      /*! A reference to the array of real data. */
      array<R>& data();
      /*! A constant reference to the array of real data. */
      const array<R>& data() const;
      /*! A pointer to the beginning of the array of real data. */
      R* begin();
      /*! A constant pointer to the beginning of the array of real data. */
      const R* begin() const;
      //@}


      //@{
      //! \name Modifying operations
      void new_constraint(const Constraint<R>& c);
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
      tribool contains(const Point<R>& point) const;

      /*! \brief A rectangle containing the polyhedron. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief An over-approximation of the polyhedron governed by the parameter \a delta.
       *
       * WARNING: The metric error of the approximation may be larger than \a delta.
       */
      Polyhedron<R> over_approximation(const R& delta) const;
      //@}
      
      //@{
      //! \name Conversion operations
      //@}

#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality. */
      friend tribool Geometry::equal<>(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
      /*! \brief Tests disjointness. */
      friend tribool Geometry::disjoint<>(const Polyhedron<R>& A, 
                                          const Polyhedron<R>& B);
        
      /*! \brief Tests disjointness. */
      friend tribool Geometry::disjoint<>(const Polyhedron<R>& A, 
                                          const Rectangle<R>& B);
        
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Rectangle<R>& A, 
                                        const Polyhedron<R>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polyhedron<R>& A, 
                                        const Polyhedron<R>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polyhedron<R>& A, 
                                        const Rectangle<R>& B);
    
      /*! \brief The generators of a polyhedron. (Deprecated) */
      friend Matrix<R> Geometry::generators<>(const Polyhedron<R>& A);
    
      //@}
      

      //@{
      //! \name Geometric binary operations
      /*! \brief The intersection of two polyhedra. */
      friend Polyhedron<R> closed_intersection<>(const Polyhedron<R>& A, 
                                          const Polyhedron<R>& B);
    
      /*! \brief The closure of the intersection of the interiors of two polyhedra. */
      friend Polyhedron<R> open_intersection<>(const Polyhedron<R>& A, 
                                                  const Polyhedron<R>& B);
    
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
    template<class R>
    class Constraint
    {
      friend class Polyhedron<R>;
      friend class PolyhedronConstraintsIterator<R>;
     public:
      /*! \brief The dimension of the constraint. */
      dimension_type dimension() const;
      /*! \brief Tests if the constraint is satisfied by a point. */
      template<class RV> tribool satisfied_by(const Point<RV>& pt) const;
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      Constraint(const dimension_type d, const R* a);
      dimension_type _d; const R* _a;
    };
    
    
    
    template<class R> 
    tribool equal(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2);

   
    template<class R> 
    tribool disjoint(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2);
    
    template<class R> 
    tribool disjoint(const Polyhedron<R>& plhd1, const Rectangle<R>& plhd2);
    
    template<class R> 
    tribool disjoint(const Rectangle<R>& plhd1, const Polyhedron<R>& plhd2);
    
    
    template<class R> 
    tribool subset(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2);

    template<class R> 
    tribool subset(const Polyhedron<R>& plhd1, const Rectangle<R>& plhd2);


    template<class R1,class R2> 
    tribool subset(const Rectangle<R1>& r, const Polyhedron<R2>& plhd);

    template<class R1,class R2> 
    tribool subset(const Zonotope<R1>& z, const Polyhedron<R2>& plhd);

    template<class R1,class R2> 
    tribool subset(const Polytope<R1>& pltp, const Polyhedron<R2>& plhd);



    template<class R> 
    Polyhedron<R> 
    open_intersection(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2);

    template<class R> 
    Polyhedron<R> 
    closed_intersection(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2) ;

    template<class R> 
    Polyhedron<R> 
    closed_intersection(const Polyhedron<R>& plhd, const Rectangle<R>& r) ;

    template<class R> 
    Polyhedron<R> 
    closed_intersection(const Rectangle<R>& r, const Polyhedron<R>& plhd) ;

    
    template<class R> 
    Polyhedron<R> 
    polyhedron(const Rectangle<R>& A) ;

    template<class R> 
    Polyhedron<typename Numeric::traits<R>::arithmetic_type> 
    polyhedron(const Polytope<R>& A) ;



    
    template<class R>
    std::ostream& operator<<(std::ostream& os, const Constraint<R>& c);
    
    template<class R>
    std::ostream& operator<<(std::ostream& os, const Polyhedron<R>& p);
    
    template<class R>
    std::istream& operator>>(std::istream& os, Polyhedron<R>& p);



  }
}

#include "polyhedron.inline.h"

#endif /* ARIADNE_POLYHEDRON_H */
