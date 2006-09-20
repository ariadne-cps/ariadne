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

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../numeric/interval.h"

#include "../geometry/ppl_polyhedron.h"
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    template<> 
    inline bool is_a<Zonotope,Zonotope>() { return true; }
    template<> 
    inline bool is_a<Zonotope,Polyhedron>() { return true; }

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Zonotope<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Zonotope<R>&);

    /*! \ingroup BasicSet
     *  \brief A zonotope of arbitrary dimension.
     * 
     * A zonotope is a set of the form \f$c+Ae\f$, where \f$||e||_{infty}\leq1\f$.
     * The intersection and membership tests are performed using algorithms from: <br>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i> (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     */
    template <typename R>
    class Zonotope {
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the zonotope. */
      typedef Point<R> state_type;
      /*! \brief The type of vectors. */
      typedef Ariadne::LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix giving principal directions. */
      typedef Ariadne::LinearAlgebra::Matrix<R> matrix_type;

     private:
      /* Zonotope's centre. */
      Rectangle<R> _central_block;
      
      /* Zonotope's principal directions. */
      matrix_type _generators;
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor constructs an empty zonotope. */
      explicit Zonotope()
        : _central_block(0),  _generators(0,0) { }
     
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and no generators. */
      explicit Zonotope(dimension_type n)
        : _central_block(n),  _generators(n,0) { }
     
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
       explicit Zonotope(dimension_type n, size_type m)
        : _central_block(n),  _generators(n,m) { }
     
      /*! \brief Construct a zonotope centred at zero from generators \a g. */
      explicit Zonotope(const matrix_type& g)
        : _central_block(g.number_of_rows()), _generators(g)
      {
        this->minimize_generators();
      }
       
     /*! \brief Construct from centre position vectorvand directions. */
      explicit Zonotope(const vector_type& c, const matrix_type& g)
        : _central_block(c), _generators(g)
      {
        if (c.size()!=g.number_of_rows()) {
          throw std::domain_error(
              "The the Matrix of principal directions does not have the same number of rows as the point dimension.");
        }
        this->minimize_generators();
      }
       
     /*! \brief Construct from centre and directions. */
      explicit Zonotope(const state_type& c, const matrix_type& g)
        : _central_block(c), _generators(g)
      {
        if (c.dimension()!=g.number_of_rows()) {
          throw std::domain_error(
              "The the Matrix of principal directions does not have the same number of rows as the point dimension.");
        }
        this->minimize_generators();
      }
       
      /*! \brief Construct from central block and directions. */
      explicit Zonotope(const Rectangle<R>& c, const LinearAlgebra::Matrix<R>& m)
        : _central_block(c), _generators(m)
      {
        if (c.dimension()!=m.number_of_rows()) {
          throw std::domain_error(
              "The the Matrix of principal directions does not have the same number of rows as the point dimension.");
        }
        this->minimize_generators();
      }
       
      /*! \brief Construct from a rectangle. */
      explicit Zonotope(const Rectangle<real_type>& r);

      /*! \brief Construct from a string literal. */
      explicit Zonotope(const std::string& s);
      
      /*! \brief Copy constructor. */
      Zonotope(const Zonotope<R>& original)
        : _central_block(original._central_block),
          _generators(original._generators)
      { }
      
      /*! \brief Copy assignment operator. */
      Zonotope<R>& operator=(const Zonotope<R>& original) {
        if(this != &original) {
          this->_central_block = original._central_block;
          this->_generators = original._generators;
        }
        return *this;
      }
      //@}
      
      
      
      //@{ 
      //! \name Data access
      /*! \brief The central block, usually used to store small roundoff errors. */
      Rectangle<R> central_block() const { return this->_central_block; }

      /*! \brief The matrix of principle directions. */
      const matrix_type& generators() const {
        return this->_generators;
      }
     
      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const {
        return this->_generators.number_of_columns();
      }

      /*! \brief The \a n th of principle direction. */
      vector_type generator(size_type n) const {
        return column(this->_generators,n);
      }
      
      //@}
      
      
      //@{
      //! \name Conversion and approximation operators
      /*! \brief Convert to a Parma Polyhedra Library closed polyhedron. */
      operator Parma_Polyhedra_Library::C_Polyhedron () const;
      
      /*! \brief Construct a parallelopic over-approximation. */
      Parallelotope<R> over_approximating_parallelotope() const;
      //@}
      

      //@{
      //! \name Geometric operations.
      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      dimension_type dimension() const {
        return this->_central_block.dimension();
      }
      
      /*! \brief True if the zonotope is empty. */
      bool empty() const { return false; }
      
      /*! \brief True if the zonotope has empty interior. */
      bool empty_interior() const;
      
      /*! \brief The centre of the zonotope. */
      state_type centre() const {
        return this->_central_block.centre();
      }
      
      /*! \brief The radius of the zonotope. */
      real_type radius() const {
        return this->bounding_box().radius();
      }
      
      /*! \brief Tests if the zonotope contains point. */
      bool contains(const state_type& point) const;
     
      /*! \brief Tests if the interior of the zonotope contains point. */
      bool interior_contains(const state_type& point) const;

      /*! \brief The vertices of the zonotope. */
      PointList<Rational> vertices() const;
      
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
      friend bool equal(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests disjointness */
      friend bool disjoint(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests disjointness */
      friend bool disjoint(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief Tests disjointness */
      friend bool disjoint(const Zonotope<R>& A, const Rectangle<R>& B);
      /*! \brief Tests intersection of interiors */
      friend bool interiors_intersect(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests intersection of interiors */
      friend bool interiors_intersect(const Rectangle<R>& A, const Zonotope<R>& B);
       /*! \brief Tests intersection of interiors */
      friend bool interiors_intersect(const Zonotope<R>& A, const Rectangle<R>& B);
      /*! \brief Tests inclusion of \a A in the interior of \a B. */
      friend bool inner_subset(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests inclusion of \a A in the interior of \a B. */
      friend bool inner_subset(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief Tests inclusion of \a A in the interior of \a B. */
      friend bool inner_subset(const Zonotope<R>& A, const Rectangle<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend bool subset(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend bool subset(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend bool subset(const Zonotope<R>& A, const Rectangle<R>& B);
      //@}
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The Minkoswi sum of two zonotopes */
      friend Zonotope<R> minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi sum of a zonotope and a rectangle. */
      friend Zonotope<R> minkowski_sum(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi sum of a rectangle and a zonotope. */
      friend Zonotope<R> minkowski_sum(const Zonotope<R>& A, const Rectangle<R>& B);
      /*! \brief The Minkoswi difference of two zonotopes */
      friend Zonotope<R> minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi difference of a rectangle and a zonotope. */
      friend Zonotope<R> minkowski_difference(const Rectangle<R>& A, const Zonotope<R>& B);
      /*! \brief The Minkoswi difference of a zonotope and a rectangle. */
      friend Zonotope<R> minkowski_difference(const Zonotope<R>& A, const Rectangle<R>& B);

      /*! \brief Adjoin generators to a zonotope. */
      friend Zonotope<R> operator+(const Zonotope<R>& z, const LinearAlgebra::Matrix<R>& G);
      //@}
#endif

      /*! \brief Computes an over approximation from an "interval zonotope". */
      static Zonotope<R> over_approximation(const Rectangle<R>& c, const LinearAlgebra::IntervalMatrix<R>& A);
      /*! \brief Scale the zonotope by a real constant. */
      static Zonotope<R> scale(const Zonotope<R>& z, const R& sf);
           
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     public:
      static bool equal(const Zonotope<R>&, const Zonotope<R>&);
      static bool disjoint(const Zonotope<R>&, const Zonotope<R>&);
      static bool disjoint(const Zonotope<R>&, const Rectangle<R>&);
      static bool interiors_intersect(const Zonotope<R>&, const Zonotope<R>&);
      static bool interiors_intersect(const Zonotope<R>&, const Rectangle<R>&);
      static bool subset(const Zonotope<R>&, const Zonotope<R>&);
      static bool subset(const Zonotope<R>&, const Rectangle<R>&);
      static bool subset(const Rectangle<R>&, const Zonotope<R>&);
      static bool inner_subset(const Zonotope<R>&, const Zonotope<R>&);
      static bool inner_subset(const Zonotope<R>&, const Rectangle<R>&);
      static bool inner_subset(const Rectangle<R>&, const Zonotope<R>&);
      static Zonotope<R> minkowski_sum(const Zonotope<R>&, const Zonotope<R>&);
      static Zonotope<R> minkowski_difference(const Zonotope<R>&, const Zonotope<R>&);
    private:
      // Minimize the generator Matrix
      void minimize_generators(void);
      
      // Order the generator Matrix by norm.
      void sort_generators(void);
      
      // The linear inequalities defining the zonotope.
      void compute_linear_inequalities(matrix_type&, vector_type&) const;

      // An extended generator is obtained by treating the errors in the central 
      // block as generators aligned to some coordinate axis.
      LinearAlgebra::Matrix<R> _extended_generators() const;
      
      // A possible vertex is the image of a vertex of the cube in 
      // generator space under the affine transformation
      std::vector< Point<Rational> > _possible_vertices() const ;

      // A possible vertex is the image of a vertex of the cube in 
      // generator space under the affine transformation
      std::vector< Point<R> > _approximate_possible_vertices() const ;
    };
  
    





    template <typename R>
    inline
    bool 
    equal(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::equal(A,B);
    }
    
    template <typename R>
    inline
    bool 
    disjoint(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::disjoint(A,B);
    }
    
    template <typename R>
    inline
    bool 
    disjoint(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::disjoint(B,A);
    }

    template <typename R>
    inline
    bool 
    disjoint(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Zonotope<R>::disjoint(A,B);
    }
    
    
    
    template <typename R>
    inline
    bool 
    interiors_intersect(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::interiors_intersect(A,B);
      
    }
   
    template <typename R>
    inline
    bool 
    interiors_intersect(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::interiors_intersect(B,A);
    }

    template <typename R>
    inline
    bool 
    interiors_intersect(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Zonotope<R>::interiors_intersect(A,B);
    }
    
    
    
    template <typename R>
    inline
    bool 
    inner_subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::inner_subset(A,B);
    }

    template <typename R>
    inline
    bool 
    inner_subset(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::inner_subset(A,B);
    }

    template <typename R>
    inline
    bool 
    inner_subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Zonotope<R>::inner_subset(A,B);
    }

    

    template <typename R>
    inline
    bool 
    subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::subset(A,B);
    }

    template <typename R>
    inline
    bool 
    subset(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::subset(A,B);
    }

    template <typename R>
    inline
    bool 
    subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Zonotope<R>::subset(A,B);
    }

    
    
    template<typename R> 
    inline
    Zonotope<R> 
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::minkowski_sum(A,B);
    }

    template<typename R> 
    inline
    Zonotope<R> 
    minkowski_sum(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::minkowski_sum(Zonotope<R>(A),B);
    }


    template<typename R> 
    inline
    Zonotope<R> 
    minkowski_sum(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Zonotope<R>::minkowski_sum(A,Zonotope<R>(B));
    }

    
    
    template<typename R> 
    inline
    Zonotope<R> 
    minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::minkowski_difference(A,B);
    } 

    template<typename R> 
    inline
    Zonotope<R> 
    minkowski_difference(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Zonotope<R>::minkowski_difference(Zonotope<R>(A),B);
    }

    template<typename R> 
    inline
    Zonotope<R> 
    minkowski_difference(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Zonotope<R>::minkowski_difference(A,Zonotope<R>(B));
    }

    template<typename R> 
    inline
    Zonotope<R> 
    operator+(const Zonotope<R>& z, const LinearAlgebra::Matrix<R>& A) 
    {
      return Zonotope<R>::minkowski_sum(z,Zonotope<R>(A));
    }

    template<typename R> 
    inline
    Zonotope<R> 
    scale(const Zonotope<R>& z, const R& sf) {
      return Zonotope<R>::scale(z,sf);
    }

    
    template<typename R> inline 
    std::ostream& operator<<(std::ostream& os, const Zonotope<R>& z) {
      return z.write(os);
    }
    
    template<typename R> inline
    std::istream& operator>>(std::istream& is, Zonotope<R>& z) {
      return z.read(is);
    }



    
  

  }
}

#endif /* _ARIADNE_ZONOTOPE_H */
