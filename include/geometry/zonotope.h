/***************************************************************************
 *            zonotope.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file zonotope.h
 *  \brief Zonotopes (affine images of cuboids).
 */

#ifndef ARIADNE_ZONOTOPE_H
#define ARIADNE_ZONOTOPE_H

#include <iosfwd>

#include "../base/iterator.h"
#include "../base/tribool.h"

#include "../numeric/declarations.h"

#include "../linear_algebra/declarations.h"
#include "../linear_algebra/matrix.h"

#include "zonotope.decl.h"

namespace Ariadne {
  namespace Geometry {

    template<class R> class Point;
    template<class R> class PointList;
    template<class E> class RectangleExpression;
    template<class R> class Rectangle;
    template<class R> class Polytope;
    template<class R> class Polyhedron;
    template<class BS> class ListSet;
      
    template<class RC,class RG=RC> class ZonotopeVerticesIterator;
      
    /*!\ingroup BasicSet
     * \brief A zonotope of arbitrary dimension.
     * 
     * A zonotope is a set of the form \f$c+Ge\f$, where \f$||e||_{\infty}\leq1\f$.
     * The columns of the matrix \f$G\f$ are the <em>generators</em> of the 
     * zonotope. 
     *
     * Zonotopes are always bounded.
     * A zonotope always contains its centre point, so can never be empty.
     * However, it may not be regular.
     *
     *
     * The intersection and membership tests may be performed using algorithms from: <br>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  
     * <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i> 
     * (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     *
     * \b Storage: A %Zonotope in dimension d with n generators is described by
     * d(n+1) real numbers. The ith component of the centre is given by the
     * ith element, and the ith component of the kth generator is given by the
     * (kd+i)th element.
     */
    template<class RC, class RG>
    class Zonotope {
      typedef typename Numeric::traits<RC>::number_type R; 
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The type used for to represent the zonotope. */
      typedef RC value_type;
      /*! \brief The real number type. */
      typedef typename Numeric::traits<R>::number_type real_type;
      /*! \brief The type of denotable point contained by the zonotope. */
      typedef Point<R> state_type;
      /*! \brief An iterator through the (possible) vertices of the zonotope. */
      typedef ZonotopeVerticesIterator<RC,RG> vertices_const_iterator;
     private:
      Point<RC> _centre;
      LinearAlgebra::Matrix<RG> _generators;
     private:
      void resize(dimension_type d, size_type m);
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor constructs an empty zonotope. */
      explicit Zonotope();
     
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and no generators. */
      explicit Zonotope(dimension_type d);
     
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
      explicit Zonotope(dimension_type d, size_type m);
     
      /*! \brief Construct from centre and directions. */
      template<class R0, class R1> 
      explicit Zonotope(const Point<R0>& c, const LinearAlgebra::Matrix<R1>& g);

      /*! \brief Construct from centre directions given by two matrices. */
      template<class R0, class R1, class R2> 
      Zonotope(const Point<R0>& c, const LinearAlgebra::Matrix<R1>& g1, const LinearAlgebra::Vector<R2>& g2);
      
      /*! \brief Construct from centre directions given by a matrix and a vector. */
      template<class R0, class R1, class R2> 
      Zonotope(const Point<R0>& c, const LinearAlgebra::Matrix<R1>& g1, const LinearAlgebra::Matrix<R2>& g2);

      /*! \brief Construct from a string literal. */
      explicit Zonotope(const std::string& s);
      
      /*! \brief Copy constructor. */
      Zonotope(const Zonotope<RC,RG>& z);
      
      /*! \brief Copy assignment operator. */
      Zonotope<RC,RG>& operator=(const Zonotope<RC,RG>& z);

      /*! \brief Type conversion copy constructor operator. */
      template<class R0,class R1> Zonotope(const Zonotope<R0,R1>& z);

      /*! \brief Type conversion assignment operator. */
      template<class R0,class R1> Zonotope<RC,RG>& operator=(const Zonotope<R0,R1>& z);

      //@}
      
      
      
      //@{ 
      //! \name Data access
      /*! \brief The ith component of the centre point. */
      const RC& centre(size_type i) const;

      /*! \brief The (i,j)-th component of the matrix of principle directions. */
      const RG& generators(size_type i, size_type j) const;
     
      /*! \brief The centre. */
      const Point<RC>& centre() const;

      /*! \brief The matrix of principle directions. */
      const LinearAlgebra::Matrix<RG>& generators() const;
     
      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const;

      /*! \brief The \a n th of principle direction. */
      LinearAlgebra::Vector<RG> generator(size_type n) const;
      //@}
      
      
      //@{
      //! \name Conversion and approximation operators
      /*! \brief Construct from a rectangle. */
      template<class E> explicit Zonotope(const RectangleExpression<E>& re);
      
      /*! \brief Assign from a rectangle expression. */
      template<class E> Zonotope<RC,RG>& operator=(const RectangleExpression<E>& re);
     
      /*! \brief Convert to a polytope. */
      operator Polytope<typename Numeric::traits<RC,RG>::arithmetic_type> () const;
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<typename Numeric::traits<RC,RG>::arithmetic_type> () const;
      //@}
      

      //@{
      //! \name Geometric operations.
      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      dimension_type dimension() const;
      
      /*! \brief True if the zonotope is empty. */
      tribool empty() const;
      
      /*! \brief Checks for boundedness. */
      tribool bounded() const;
      
      /*! \brief The radius of the zonotope. */
      R radius() const;
      
      /*! \brief Tests if the zonotope contains point. */
      template<class R> tribool contains(const Point<R>& point) const;

      /*! \brief The vertices of the zonotope. */
      PointList<typename Numeric::traits<RC,RG>::arithmetic_type> vertices() const;
      
      /*! \brief An iterator to the beginning of the (possible) vertices. */
      vertices_const_iterator vertices_begin() const;
      
      /*! \brief An iterator to the end of the (possible) vertices. */
      vertices_const_iterator vertices_end() const;
      
      /*! \brief Subdivide into two smaller pieces. */
      ListSet< Zonotope<RC,RG> > divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet< Zonotope<RC,RG> > subdivide() const;
      
      /*! \brief A rectangle containing the given zonotope. */
      Rectangle<R> bounding_box() const;
      
      //@}
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality */
      friend tribool equal(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Rectangle<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Zonotope<RC,RG>& A, const Rectangle<RC,RG>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Rectangle<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Zonotope<RC,RG>& A, const Rectangle<RC,RG>& B);
      //@}
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The Minkoswi sum of two zonotopes */
      friend Zonotope<F> minkowski_sum(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief The Minkoswi sum of a zonotope and a rectangle. */
      friend Zonotope<F> minkowski_sum(const Rectangle<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief The Minkoswi sum of a rectangle and a zonotope. */
      friend Zonotope<F> minkowski_sum(const Zonotope<RC,RG>& A, const Rectangle<RC,RG>& B);
      /*! \brief The Minkoswi difference of two zonotopes */
      friend Zonotope<F> minkowski_difference(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief The Minkoswi difference of a rectangle and a zonotope. */
      friend Zonotope<F> minkowski_difference(const Rectangle<RC,RG>& A, const Zonotope<RC,RG>& B);
      /*! \brief The Minkoswi difference of a zonotope and a rectangle. */
      friend Zonotope<F> minkowski_difference(const Zonotope<RC,RG>& A, const Rectangle<RC,RG>& B);

      /*! \brief Adjoin generators to a zonotope. */
      friend Zonotope<RC,RG> operator+(const Zonotope<RC,RG>& z, const LinearAlgebra::Matrix<RG>& G);
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
     private:
      // Minimize the generator Matrix
      void minimize_generators(void);
      
      // Order the generator Matrix by norm.
      void sort_generators(void);
      
      // A possible vertex is the image of a vertex of the cube in 
      // generator space under the affine transformation
      std::vector< Point<F> > _possible_vertices() const ;

      // A possible vertex is the image of a vertex of the cube in 
      // generator space under the affine transformation
      std::vector< Point<RC> > _approximate_possible_vertices() const ;
    };
  

    template<class R,class RC,class RG> tribool contains(const Zonotope<RC,RG>& z,const Point<R>& pt);
    

    template<class R,class RC,class RG> tribool disjoint(const Rectangle<R>& A, const Zonotope<RC,RG>& B);
    
    template<class R,class RC,class RG> tribool disjoint(const Zonotope<RC,RG>& A, const Rectangle<R>& B);
    
    template<class RC,class RG> tribool disjoint(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
    


    template<class R,class RC,class RG> tribool subset(const Rectangle<R>& A, const Zonotope<RC,RG>& B);
    
    template<class R,class RC,class RG> tribool subset(const Zonotope<RC,RG>& A, const Rectangle<R>& B);
    
    template<class RC,class RG> tribool subset(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);


    template<class R> Zonotope<Numeric::Interval<R>,R> over_approximation(const Zonotope< Numeric::Interval<R>, Numeric::Interval<R> >&);
     
    template<class R> Zonotope<R,R> over_approximation(const Zonotope< Numeric::Interval<R>, R >&);
     
          

    template<class R> Zonotope<R> approximation(const Zonotope< Numeric::Interval<R> >&);

    template<class R> Zonotope<R> approximation(const Zonotope< Numeric::Interval<R>, R >&);

    template<class R> Zonotope<R> approximation(const Zonotope<R>&);

    
    template<class RC,class RG> 
    Zonotope<typename Numeric::traits<RC>::arithmetic_type,RG>
    minkowski_sum(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
    
    template<class RC,class RG> 
    Zonotope<typename Numeric::traits<RC>::arithmetic_type,RG>
    minkowski_difference(const Zonotope<RC,RG>& A, const Zonotope<RC,RG>& B);
    
    
    
    template<class R,class RC,class RG>
    Zonotope<typename Numeric::traits<R,RG>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
    minkowski_sum(const Rectangle<R>& A, const Zonotope<RC,RG>& B);
    
    template<class R,class RC,class RG>
    Zonotope<typename Numeric::traits<R,RC>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
    minkowski_sum(const Zonotope<RC,RG>& A, const Rectangle<R>& B);
    
    template<class R,class RC,class RG>
    Zonotope<typename Numeric::traits<R,RG>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
    minkowski_difference(const Rectangle<R>& A, const Zonotope<RC,RG>& B);
    
    template<class R,class RC,class RG>
    Zonotope<typename Numeric::traits<R,RC>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
    minkowski_difference(const Zonotope<RC,RG>& A, const Rectangle<R>& B);

    
    template<class RC,class RG> 
    ListSet< Zonotope<RC,RG> >
    subdivide(const Zonotope<RC,RG>& z);
    
    template<class RC,class RG> 
    ListSet< Zonotope<RC,RG> >
    divide(const Zonotope<RC,RG>& z);
    
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    zonotope(const Rectangle<R>& r);
    
    template<class RC,class RG> 
    Polytope<typename Numeric::traits<RC,RG>::arithmetic_type> 
    polytope(const Zonotope<RC,RG>& r);
    
    template<class RC,class RG> 
    Polyhedron<typename Numeric::traits<RC,RG>::arithmetic_type> 
    polyhedron(const Zonotope<RC,RG>& r);
    
    
    template<class RC,class RG> 
    Zonotope<RC,RG> 
    operator+(const Zonotope<RC,RG>& z, const LinearAlgebra::Matrix<RG>& A);
    
    template<class RC,class RG>  
    std::ostream& operator<<(std::ostream& os, const Zonotope<RC,RG>& z);
    
    template<class RC,class RG> 
    std::istream& operator>>(std::istream& is, Zonotope<RC,RG>& z);

    
    
    
    template<class RC,class RG> class ZonotopeVerticesIterator;
      
    template<class RC,class RG> 
    std::ostream& operator<<(std::ostream& os, const ZonotopeVerticesIterator<RC,RG>& iter);



    
  

  }
}

#include "zonotope.inline.h"
#include "zonotope.template.h"

#endif /* ARIADNE_ZONOTOPE_H */
