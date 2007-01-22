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

#ifndef _ARIADNE_ZONOTOPE_H
#define _ARIADNE_ZONOTOPE_H

#include <iosfwd>

#include "../base/iterator.h"
#include "../base/tribool.h"

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    template<class R> class ZonotopeVerticesIterator;
      
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
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i> (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     *
     * \b Storage: A %Zonotope in dimension d with n generators is described by
     * d(n+1) real numbers. The ith component of the centre is given by the
     * ith element, and the ith component of the kth generator is given by the
     * (kd+i)th element.
     */
    template<class R>
    class Zonotope {
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the zonotope. */
      typedef Point<R> state_type;
      /*! \brief An iterator through the (possible) vertices of the zonotope. */
      typedef ZonotopeVerticesIterator<R> vertices_const_iterator;
     private:
      dimension_type _dimension;
      size_type _number_of_generators;
      array<R> _data;
     private:
      array<R>& data();
      size_type data_size() const;
      void resize(dimension_type d, size_type m);
      LinearAlgebra::VectorSlice<R> _centre();
      LinearAlgebra::MatrixSlice<R> _generators();
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
      template<class R1, class R2> 
      explicit Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g);

      /*! \brief Construct from centre directions given by two matrices. */
      template<class R1, class R2, class R3> 
      Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Vector<R3>& g2);
      
      /*! \brief Construct from centre directions given by a matrix and a vector. */
      template<class R1, class R2, class R3> 
      Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Matrix<R3>& g2);

      /*! \brief Construct from a string literal. */
      explicit Zonotope(const std::string& s);
      
      /*! \brief Copy constructor. */
      template<class R1> Zonotope(const Zonotope<R1>& z);
      
      /*! \brief Copy assignment operator. */
      template<class Rl> Zonotope<R>& operator=(const Zonotope<Rl>& z);
      //@}
      
      
      
      //@{ 
      //! \name Data access
      /*! \brief The raw data of the zonotope. */
      const array<R>& data() const;

      /*! \brief The centre. */
      Point<R> centre() const;

      /*! \brief The matrix of principle directions. */
      const LinearAlgebra::MatrixSlice<R> generators() const;
     
      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const;

      /*! \brief The \a n th of principle direction. */
      LinearAlgebra::Vector<R> generator(size_type n) const;
      //@}
      
      
      //@{
      //! \name Conversion and approximation operators
      /*! \brief Construct from a rectangle. */
      template<class Rl> explicit Zonotope(const Rectangle<Rl>& r);
      
      /*! \brief Assign from a rectangle. */
      template<class Rl> Zonotope<R>& operator=(const Rectangle<Rl>& z);
     
      /*! \brief Convert to a polytope. */
      operator Polytope<typename Numeric::traits<R>::arithmetic_type> () const;
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<typename Numeric::traits<R>::arithmetic_type> () const;
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
      tribool contains(const Point<R>& point) const;

      /*! \brief The vertices of the zonotope. */
      PointList<typename Numeric::traits<R>::arithmetic_type> vertices() const;
      
      /*! \brief An iterator to the beginning of the (possible) vertices. */
      vertices_const_iterator vertices_begin() const;
      
      /*! \brief An iterator to the end of the (possible) vertices. */
      vertices_const_iterator vertices_end() const;
      
      /*! \brief Subdivide into two smaller pieces. */
      ListSet<R,Geometry::Zonotope> divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet<R,Geometry::Zonotope> subdivide() const;
      
      /*! \brief A rectangle containing the given zonotope. */
      Rectangle<typename Numeric::traits<R>::number_type> bounding_box() const;
      
      /*! \brief Scale the zonotope by a real constant. */
      static Zonotope<F> scale(const Zonotope<R>& z, const R& sf);
           
      //@}
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality */
      friend tribool equal(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Zonotope<R>& A, const Rectangle<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Zonotope<R>& A, const Rectangle<R>& B);
      //@}
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The Minkoswi sum of two zonotopes */
      friend Zonotope<F> minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi sum of a zonotope and a rectangle. */
      friend Zonotope<F> minkowski_sum(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi sum of a rectangle and a zonotope. */
      friend Zonotope<F> minkowski_sum(const Zonotope<R>& A, const Rectangle<R>& B);
      /*! \brief The Minkoswi difference of two zonotopes */
      friend Zonotope<F> minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi difference of a rectangle and a zonotope. */
      friend Zonotope<F> minkowski_difference(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi difference of a zonotope and a rectangle. */
      friend Zonotope<F> minkowski_difference(const Zonotope<R>& A, const Rectangle<R>& B);

      /*! \brief Adjoin generators to a zonotope. */
      friend Zonotope<R> operator+(const Zonotope<R>& z, const LinearAlgebra::Matrix<R>& G);
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
      // Minimize the generator Matrix
      void minimize_generators(void);
      
      // Order the generator Matrix by norm.
      void sort_generators(void);
      
      // A possible vertex is the image of a vertex of the cube in 
      // generator space under the affine transformation
      std::vector< Point<F> > _possible_vertices() const ;

      // A possible vertex is the image of a vertex of the cube in 
      // generator space under the affine transformation
      std::vector< Point<R> > _approximate_possible_vertices() const ;
    };
  


    template<class R> tribool disjoint(const Rectangle<R>& A, const Zonotope<R>& B);
    
    template<class R> tribool disjoint(const Zonotope<R>& A, const Rectangle<R>& B);
    
    template<class R> tribool disjoint(const Zonotope<R>& A, const Zonotope<R>& B);
    
    template<class R> tribool subset(const Rectangle<R>& A, const Zonotope<R>& B);
    
    template<class R> tribool subset(const Zonotope<R>& A, const Rectangle<R>& B);
    
    template<class R> tribool subset(const Zonotope<R>& A, const Zonotope<R>& B);
    
    template<class R> Zonotope<R> over_approximation(const Zonotope<R>&);
 
    template<class R> Zonotope<R> over_approximation(const Zonotope< Interval<R> >&);
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B);
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type>
    minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B);
    
    
    
    template<class R>
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_sum(const Rectangle<R>& A, const Zonotope<R>& B);
    
    template<class R>
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_sum(const Zonotope<R>& A, const Rectangle<R>& B);
    
    template<class R>
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Rectangle<R>& A, const Zonotope<R>& B);

    template<class R>
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Zonotope<R>& A, const Rectangle<R>& B);
    
    
    template<class R> 
    ListSet<R,Zonotope>
    subdivide(const Zonotope<R>& z);
    
    template<class R> 
    ListSet<Interval<R>,Zonotope>
    subdivide(const Zonotope< Interval<R> >& z);
    
    template<class R> 
    ListSet<R,Zonotope>
    divide(const Zonotope<R>& z);
    
    template<class R> 
    ListSet<Interval<R>,Zonotope>
    divide(const Zonotope< Interval<R> >& z);
    
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    zonotope(const Rectangle<R>& r);
    
    template<class R> 
    Polytope<typename Numeric::traits<R>::arithmetic_type> 
    polytope(const Zonotope<R>& r);
    
    template<class R> 
    Polyhedron<typename Numeric::traits<R>::arithmetic_type> 
    polyhedron(const Zonotope<R>& r);
    
    
    template<class R> 
    Zonotope<R> 
    operator+(const Zonotope<R>& z, const LinearAlgebra::Matrix<R>& A);
    
    template<class R> 
    Zonotope<typename Zonotope<R>::F> 
    scale(const Zonotope<R>& z, const R& sf);
    
    template<class R>  
    std::ostream& operator<<(std::ostream& os, const Zonotope<R>& z);
    
    template<class R> 
    std::istream& operator>>(std::istream& is, Zonotope<R>& z);

    
    
    
    template<class R> class ZonotopeVerticesIterator;
      
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const ZonotopeVerticesIterator<R>& iter);



    
  

  }
}

#include "zonotope.inline.h"
#include "zonotope.template.h"

#endif /* _ARIADNE_ZONOTOPE_H */
