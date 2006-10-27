/***************************************************************************
 *            zonotope.h
 *
 *  6 February 2006
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

    template<> 
    inline bool is_a<Zonotope,Zonotope>() { return true; }
    template<> 
    inline bool is_a<Zonotope,Polyhedron>() { return true; }

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
      /*! \brief The type of vectors. */
      typedef Ariadne::LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix giving principal directions. */
      typedef Ariadne::LinearAlgebra::Matrix<R> matrix_type;
      /*! \brief An iterator through the (possible) vertices of the zonotope. */
      typedef ZonotopeVerticesIterator<R> vertices_iterator;
      typedef ZonotopeVerticesIterator<R> vertices_const_iterator;

     private:
      /* Zonotope's centre. */
      Point<R> _centre;
      
      /* Zonotope's principal directions. */
      matrix_type _generators;
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor constructs an empty zonotope. */
      explicit Zonotope()
        : _centre(0),  _generators(0,0) { }
     
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and no generators. */
      explicit Zonotope(dimension_type n)
        : _centre(n),  _generators(n,0) { }
     
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
       explicit Zonotope(dimension_type n, size_type m)
        : _centre(n),  _generators(n,m) { }
     
     /*! \brief Construct from centre position vectorvand directions. */
      explicit Zonotope(const LinearAlgebra::Vector<R>& c, const LinearAlgebra::Matrix<R>& g)
        : _centre(c), _generators(g)
      {
        if (c.size()!=g.number_of_rows()) {
          throw std::domain_error(
              "The the Matrix of principal directions does not have the same number of rows as the point dimension.");
        }
        this->minimize_generators();
      }
       
      /*! \brief Construct from centre and directions. */
      explicit Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& g)
        : _centre(c), _generators(g)
      {
        if (c.dimension()!=g.number_of_rows()) {
          throw std::domain_error(
              "The the Matrix of principal directions does not have the same number of rows as the point dimension.");
        }
        this->minimize_generators();
      }

      /*! \brief Construct from a string literal. */
      explicit Zonotope(const std::string& s);
      
      /*! \brief Copy constructor. */
      template<class Rl>
      Zonotope(const Zonotope<Rl>& original)
        : _centre(original.centre()),
          _generators(original.generators())
      { }
      
      /*! \brief Copy assignment operator. */
      Zonotope<R>& operator=(const Zonotope<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      //@}
      
      
      
      //@{ 
      //! \name Data access
      /*! \brief The centre. */
      Point<R> centre() const { return this->_centre; }

      /*! \brief The matrix of principle directions. */
      const LinearAlgebra::Matrix<R>& generators() const {
        return this->_generators;
      }
     
      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const {
        return this->_generators.number_of_columns();
      }

      /*! \brief The \a n th of principle direction. */
      LinearAlgebra::Vector<R> generator(size_type n) const {
        return column(this->_generators,n);
      }
      
      //@}
      
      
      //@{
      //! \name Conversion and approximation operators
       
      /*! \brief Convert from a rectangle. */
      Zonotope(const Rectangle<real_type>& r);
      
      /*! \brief Convert to a polytope. */
      //operator Polytope<typename Numeric::traits<R>::arithmetic_type> () const;
      /*! \brief Convert to a polytope. */
      operator Polytope<Rational> () const;
      
      /*! \brief Convert to a polyhedron. */
      //operator Polyhedron<typename Numeric::traits<R>::arithmetic_type> () const;
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<Rational> () const;
      
      /*! \brief Construct an over-approximating zonotope from an interval zonotope. */
      static Zonotope<R> over_approximation(const Zonotope<I>& iz);
      //@}
      

      //@{
      //! \name Geometric operations.
      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      dimension_type dimension() const {
        return this->_centre.dimension();
      }
      
      /*! \brief True if the zonotope is empty. */
      tribool empty() const;
      
      /*! \brief The radius of the zonotope. */
      real_type radius() const {
        return this->bounding_box().radius();
      }
      
      /*! \brief Tests if the zonotope contains point. */
      tribool contains(const Point<R>& point) const;

      /*! \brief The vertices of the zonotope. */
      PointList<Rational> vertices() const;
      
      /*! \brief An iterator to the beginning of the (possible) vertices. */
      vertices_const_iterator vertices_begin() const;
      
      /*! \brief An iterator to the end of the (possible) vertices. */
      vertices_const_iterator vertices_end() const;
      
      /*! \brief Subdivide into two smaller pieces. */
      ListSet<R,Geometry::Zonotope> divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet<R,Geometry::Zonotope> subdivide() const;
      
      /*! \brief A rectangle containing the given zonotope. */
      Rectangle<R> bounding_box() const;
      
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

      /*! \brief Scale the zonotope by a real constant. */
      static Zonotope<F> scale(const Zonotope<R>& z, const R& sf);
           
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
  

    template<class R>
    class Zonotope< Interval<R> >
    {
      typedef Interval<R> I;
     public:
      template<class R1, class R2> Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g)
        : _centre(c), _generators(g) { }
      dimension_type dimension() const { return _centre.dimension(); }
      const Point<I>& centre() const { return _centre; }
      const LinearAlgebra::Matrix<I>& generators() const { return _generators; }
      
      Zonotope<R> over_approximation() const;
      
      std::ostream& write(std::ostream& os) const;
     private:
      Point<I> _centre;
      LinearAlgebra::Matrix<I> _generators;
    };



    //template<class R>
    //inline
    //tribool 
    //equal(const Zonotope<R>& A, const Zonotope<R>& B) 
    //{
    //  return Zonotope<R>::equal(A,B);
    //}
    
    template<class R> tribool disjoint(const Rectangle<R>& A, const Zonotope<R>& B);
    
    template<class R> tribool disjoint(const Zonotope<R>& A, const Rectangle<R>& B);
    
    template<class R> tribool disjoint(const Zonotope<R>& A, const Zonotope<R>& B);
    
    template<class R> tribool subset(const Rectangle<R>& A, const Zonotope<R>& B);
    
    template<class R> tribool subset(const Zonotope<R>& A, const Rectangle<R>& B);
    
    template<class R> tribool subset(const Zonotope<R>& A, const Zonotope<R>& B);
    
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B);
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type>
    minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B);
    
    
    
    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_sum(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::minkowski_sum(Zonotope<R>(A),B);
    }

    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_sum(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::minkowski_sum(A,Zonotope<R>(B));
    }

    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::minkowski_difference(Zonotope<R>(A),B);
    }

    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::minkowski_difference(A,Zonotope<R>(B));
    }

    
    template<class R> 
    inline
    Zonotope<R> 
    operator+(const Zonotope<R>& z, const LinearAlgebra::Matrix<R>& A) 
    {
      return Zonotope<R>(z.centre(),concatenate_columns(z.generators(),A));
    }

    template<class R> 
    inline
    Zonotope<typename Zonotope<R>::F> 
    scale(const Zonotope<R>& z, const R& sf) {
      return Zonotope<R>::scale(z,sf);
    }

    
    template<class R>
    class ZonotopeVerticesIterator 
      : public boost::iterator_facade<ZonotopeVerticesIterator<R>,
                                      Point<typename Numeric::traits<R>::arithmetic_type>,
                                      boost::forward_traversal_tag,
                                      Point<typename Numeric::traits<R>::arithmetic_type> const&,
                                      Point<typename Numeric::traits<R>::arithmetic_type> const*
                                     >
    {
      friend class Zonotope<R>;
      typedef typename Numeric::traits<R>::arithmetic_type F;
      const Zonotope<R>* _z; long unsigned int _i; bool _parity; Point<F> _pt;
     public:
      ZonotopeVerticesIterator(const Zonotope<R>& z, bool end) 
        : _z(&z), _i(end ? (1u<<(z.number_of_generators()-1))*3 : 0), _parity(0), _pt(z.centre())
      {
        if(end) { return; }
        for(dimension_type i=0; i!=z.dimension(); ++i) {
          for(dimension_type j=0; j!=z.number_of_generators(); ++j) {
            this->_pt[i]-=z.generators()(i,j); } }
      }
      bool equal(const ZonotopeVerticesIterator<R>& other) const {
        //std::cerr << "ZonotopeVerticesIterator<R>::equal" << std::endl;
        return this->_i==other._i && this->_z==other._z; }
      const Point<F>& dereference() const { 
        //std::cerr << "ZonotopeVerticesIterator<R>::dereference" << std::endl;
        return this->_pt; }
      void increment() { 
        //std::cerr << "ZonotopeVerticesIterator<R>::increment" << std::endl;
        uint j=0; uint m=1; if(this->_parity) { while(!(m&(this->_i))) { ++j; m*=2u; } ++j; m*=2u; }
        this->_parity=!this->_parity;
        if(j==this->_z->number_of_generators()) { this->_i+=m; return; }
        if(m&(this->_i)) { this->_pt=this->_pt-R(2)*this->_z->generator(j); this->_i-=m; }
        else { this->_pt=this->_pt+R(2)*this->_z->generator(j); this->_i+=m; }
      }

      std::ostream& write(std::ostream& os) const {
        return os << _z << " " << _i << " " << _parity << std::endl; }
    }; 

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ZonotopeVerticesIterator<R>& iter) {
      return iter.write(os);
    }

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const Zonotope<R>& z) {
      return z.write(os);
    }
    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Zonotope<R>& z) {
      return z.read(is);
    }



    
  

  }
}

#endif /* _ARIADNE_ZONOTOPE_H */
