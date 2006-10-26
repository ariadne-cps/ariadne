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

#include "../base/tribool.h"
#include "../base/iterator.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"

namespace Ariadne {  
  namespace Geometry {

    template<class R> class PolyhedronConstraintsIterator;
    
    /*! \ingroup BasicSet
     *  \brief A polyhedron (not necessarily bounded polyhedral set) described by a system of linear inequalities.
     */ 
    template<class R>
    class Polyhedron {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real numbers used to describe the polyhedron. */
      typedef Rational real_type;
      /*! \brief The type of denotable point contained by the polyhedron. */
      typedef Point<R> state_type;
      /*! \brief The type of vector. */
      typedef Ariadne::LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix. */
      typedef Ariadne::LinearAlgebra::Matrix<R> matrix_type;
      /*! \brief An iterator over the constraints of the Polyhedron. */
      typedef PolyhedronConstraintsIterator<R> constraints_const_iterator;
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Construct full Euclidean space of dimension \a n.
       */
      explicit Polyhedron<R>(dimension_type n=0);
     
      /*! \brief Construct the polyhedron defined by the matrix equations \f$Ax\leq b\f$.
       */
      explicit Polyhedron<R>(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b);
            
      /*! \brief Construct from a list of vertices. */
      explicit Polyhedron<R>(const PointList<R>& vertices);
            
      /*! \brief Convert from a rectangle. */
      Polyhedron<R>(const Rectangle<R>& rect);
            
      /*! \brief Convert from a polytope. */
      Polyhedron<R>(const Polytope<R>& plyt);
            
      /*! \brief Copy constructor. 
       */
      template<class Rl> inline Polyhedron(const Polyhedron<Rl>& original)
        : _A(original.A()), _b(original.b()) { }
          
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
      
      
      //@{
      //! \name Data access
      /*! \brief The matrix \f$A\f$ in the inequalities \f$Ax\leq b\f$. */
      const LinearAlgebra::Matrix<R>& A() const { return this->_A; }
      /*! \brief The vector \f$b\f$ in the inequalities \f$Ax\leq b\f$. */
      const LinearAlgebra::Vector<R>& b() const { return this->_b; }
      /*! \brief An iterator to the beginning of the constraints. */
      size_type number_of_constraints() const { return _A.number_of_rows(); }
      /*! \brief An iterator to the beginning of the constraints. */
      constraints_const_iterator constraints_begin() const;
      /*! \brief An iterator to the end of the constraints. */
      constraints_const_iterator constraints_end() const;
      
      
      //@}


      //@{
      //! \name Geometric operations
      /*! \brief Returns the polyhedron's space dimension.
       */
      dimension_type dimension() const;
      
      /*! \brief The vertices of the polyhedron. */
      PointList<Rational>  vertices() const; 
           
      /*! \brief Checks for emptyness.
       */
      tribool empty() const;
      
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
        
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool Geometry::subset<>(const Polyhedron<R>& A, 
                                     const Polyhedron<R>& B);
    
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
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     private:
      static void _instantiate_geometry_operators();
     private:
       LinearAlgebra::Matrix<R> _A;
       LinearAlgebra::Vector<R> _b;
    };
   

    template<class R>
    class Polyhedron< Interval<R> >
    {
      typedef Interval<R> I;
     public:
      template<class Rl1, class Rl2>
      Polyhedron(const LinearAlgebra::Matrix<Rl1> A, const LinearAlgebra::Vector<Rl2> b) 
        : _A(A), _b(b) { }
      template<class Rl>
      Polyhedron(const PointList<Rl> pts);
        
      dimension_type dimension() { return this->_A.number_of_columns(); }
      size_type number_of_constraints() { return this->_A.number_of_rows(); }
      
      const LinearAlgebra::Matrix<I>& A() const { return this->_A; }
      const LinearAlgebra::Vector<I>& b() const { return this->_b; }
      
      tribool contains(const Point<R>& pt) const;
      tribool contains(const Point<I>& pt) const;
      
      static tribool subset(const Polytope<I>& pltp, const Polyhedron<I>& plhd);
     private:
      LinearAlgebra::Matrix<I> _A;
      LinearAlgebra::Vector<I> _b;
    };
  

    
    /*! \brief A linear constraint. */
    template<class R>
    class Constraint
    {
      friend class Polyhedron<R>;
      friend class PolyhedronConstraintsIterator<R>;
     public:
      /*! \brief Write to an output stream. */
      template<class RV> tribool satisfied_by(const Point<RV>& pt) const;
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
     public:
      Constraint(const dimension_type d, const R* a, const R& b) : _d(d), _a(a), _b(&b) { }
      dimension_type _d; const R* _a; const R* _b;
    };

    template<class R1> template<class R2> inline 
    tribool Constraint<R1>::satisfied_by(const Point<R2>& pt) const
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type F;
      F prod=0;
      for(dimension_type i=0; i!=pt.dimension(); ++i) {
        prod+=this->_a[i]*pt[i]; 
      }
      return prod<=*this->_b;
    }
    
    template<class R>
    class PolyhedronConstraintsIterator
      : public boost::iterator_facade<PolyhedronConstraintsIterator<R>,
                                      Constraint<R>,
                                      boost::forward_traversal_tag,
                                      const Constraint<R>&,
                                      const Constraint<R>*
                                     >
    {
     public:
      PolyhedronConstraintsIterator(const Polyhedron<R>& ply, const size_type& n)
        : _c(ply.dimension(),ply.A().begin()+n*ply.dimension(),ply.b().begin()[n]) { }
      bool equal(const PolyhedronConstraintsIterator<R>& other) const { 
        return this->_c._a==other._c._a; }
      const Constraint<R>& dereference() const { return _c; }
      void increment() { _c._a+=_c._d; _c._b+=1; }
     private:
      Constraint<R> _c;
    };
     
      
   
    
   
   


    template<class R> 
    tribool equal(const Polyhedron<R>& A, const Polyhedron<R>& B);
   
    template<class R> 
    tribool disjoint(const Polyhedron<R>& A, const Polyhedron<R>& B);
    
    template<class R> 
    tribool disjoint(const Polyhedron<R>& A, const Rectangle<R>& B);
    
    template<class R> 
    tribool disjoint(const Rectangle<R>& A, const Polyhedron<R>& B);
    
    
    template<class R> 
    tribool subset(const Polyhedron<R>& A, const Polyhedron<R>& B);
    template<class R> 
    tribool subset(const Polyhedron<R>& A, const Rectangle<R>& B);

    template<class R1,class R2> 
    tribool subset(const Rectangle<R1>& A, const Polyhedron<R2>& B);
    template<class R1,class R2> 
    tribool subset(const Zonotope<R1>& A, const Polyhedron<R2>& B);
    template<class R1,class R2> 
    tribool subset(const Polytope<R1>& A, const Polyhedron<R2>& B);

    template<class R> 
    Polyhedron<R> 
    open_intersection(const Polyhedron<R>& A, const Polyhedron<R>& B);

    template<class R> 
    Polyhedron<R> 
    closed_intersection(const Polyhedron<R>& A, const Polyhedron<R>& B) ;

    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Constraint<R>& c) {
      return c.write(os); }
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Polyhedron<R>& p) {
      return p.write(os); }
    template<class R> inline
    std::istream& operator>>(std::istream& os, Polyhedron<R>& p) {
      return p.read(os); }
   



  }
}

#endif /* _ARIADNE_POLYHEDRON_H */
