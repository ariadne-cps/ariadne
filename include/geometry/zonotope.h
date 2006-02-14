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
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include <ariadne.h>

#include <linear_algebra/linear_algebra.h>

#include <base/utility.h>
#include <base/interval.h>

#include <geometry/point.h>
#include <geometry/rectangle.h>
#include <geometry/list_set.h>
#include <geometry/polyhedron.h>

namespace Ariadne {
  namespace Geometry {
    template < typename R > class Zonotope;

    template < typename R > class Rectangle;
    template < typename R > class Polyhedron;
    template < typename R, template <typename> class BS > class ListSet;

    template<typename R> Zonotope<R> minkowski_sum(const Zonotope<R>& A, const Zonotope& b);
    template<typename R> Zonotope<R> minkowski_difference(const Zonotope<R>& A, const Zonotope& b);

    template<typename R> bool interiors_intersect(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool disjoint(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool inner_subset(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool subset(const Zonotope<R>& A, const Zonotope<R>& B);

    template<typename R> bool subset_of_open_cover(const Zonotope<R>& A, const ListSet<R, Zonotope >& U);
    template<typename R> bool inner_subset(const Zonotope<R>& A, const ListSet<R,Zonotope>& U);
    template<typename R> bool subset(const Zonotope<R>& A, const ListSet<R,Zonotope>& B);
    
    template<typename R> std::ostream& operator<<(std::ostream&, const Zonotope<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Zonotope<R>&);

    /*! \brief A zonotope of arbitrary dimension.
     * 
     * A zonotope is a set of the form \f$c+Ae\f$, where \f$||e||_{infty}\leq1\f$.
     * The intersection and membership tests are performed using algorithms from: <br/>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  <it>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</it> (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     */
    template <typename R>
    class Zonotope {
       /*! \brief Tests intersection of interiors. */
      friend bool interiors_intersect <> (const Parallelopiped<R>& A,
                                          const Zonotope<R>& B);

       /*! \brief Tests disjointness */
      friend bool disjoint <> (const Zonotope<R>& A,
                               const Zonotope<R>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a B. */
      friend bool inner_subset <> (const Zonotope<R>& A,
                                   const Zonotope<R>& B);

      /*! \brief Tests inclusion. */
      friend bool subset <> (const Zonotope<R>& A,
                             const Zonotope<R>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a B. */
      friend bool inner_subset <> (const Zonotope<R>& A,
                                   const ListSet<R,::Ariadne::Geometry::Zonotope>& B);

      /*! \brief Tests if \a A is a subset of \a B. */
      friend bool subset <> (const Zonotope<R>& A,
                             const ListSet<R,::Ariadne::Geometry::Zonotope>& B);


      /*! \brief Tests inclusion in an open cover, represented as a ListSet.
       */
      friend bool subset_of_open_cover <> (const Zonotope<R>& A,
                                           const ListSet<R,::Ariadne::Geometry::Zonotope>& B);

     public:
      /*! \brief The unsigned integer type used to denote the array positions. */
      typedef size_t size_type;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R Real;
      /*! \brief The type of denotable state contained by the rectangle. */
      typedef Point<R> State;
      /*! \brief The type of matrix giving principal directions. */
      typedef ::Ariadne::LinearAlgebra::vector<R> Vector;
      /*! \brief The type of matrix giving principal directions. */
      typedef ::Ariadne::LinearAlgebra::matrix<R> Matrix;
     private:
      /* Zonotope's centre. */
      State _centre;
      
      /* Zonotope's principal directions. */
      Matrix _generators;
      
     public:
      /*! \brief Default constructor constructs an empty parallelopiped of dimension \a n. */
      explicit Zonotope(size_type n = 0)
        : _centre(n),  _generators(n,0) 
      {
      }
      
      /*! \brief Construct from centre and directions. */
      explicit Zonotope(const State& c, const Matrix& m)
        : _centre(c), _generators(m)
      {
        if (c.dimension()!=LinearAlgebra::number_of_rows(m)) {
          throw std::domain_error(
              "The the matrix of principal directions does not have the same number of rows as the state dimension.");
        }
        
      }
      
      /*! \brief Construct from a Rectangle. */
      explicit Zonotope(const Rectangle<Real>& r);
      
      /*! \brief Construct from a string literal. */
      explicit Zonotope(const std::string& s)
        : _centre(), _generators()
      {
        std::stringstream ss(s);
        ss >> *this;
      }
      
      /*! \brief Copy constructor. */
      Zonotope(const Zonotope<R>& original)
        : _centre(original._centre),
          _generators(original._generators)
      { }
      
      /*! \brief Copy assignment operator. */
      Zonotope<R>& operator=(const Zonotope<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      
      /*! \brief A rectangle containing the given parallelopiped. */
      inline Rectangle<R> bounding_box() const {
        Vector offset(this->dimension());
        for(size_type i=0; i!=this->dimension(); ++i) {
          for(size_type j=0; i!=this->dimension(); ++j) {
            offset[i] += abs(this->_generators(i,j));
          }
        }
        return Rectangle<R>(this->centre()+offset, this->centre()-offset);
      }
      
      /*! \brief Convert to a polyhedron. */
      inline operator Polyhedron<R> () const;
      
      /*! \brief The dimension of the Euclidean space the parallelopiped lies in. */
      inline size_type dimension() const {
        return (this->_centre).dimension();
      }
      
      /*! \brief True if the parallelopiped is empty. */
      inline bool empty() const {
        throw std::domain_error("Zonotope::empty() not implemented.");
      }
      
      /*! \brief True if the parallelopiped has empty interior. */
      inline bool empty_interior() const {
        throw std::domain_error("Zonotope::empty_interior() not implemented.");
      }
      
      /*! \brief The centre of the parallelopiped. */
      inline State centre() const {
        return this->_centre;
      }
      
      /*! \brief The \a n th of principle direction. */
      inline Vector principle_direction(size_type n) const {
        return matrix_column(this->_generators);
      }
      
      /*! \brief The matrix of principle directions. */
      inline Matrix principle_directions() const {
        return this->_generators;
      }
      
      /*! \brief Tests if the parallelopiped contains \a state. */
      inline bool contains(const State& state) const {
        if (state.dimension()!=this->dimension()) {
          throw std::domain_error("This object and parameter have different space dimensions");
        }  
        
        if (this->empty()) { return false; }
          
        Vector v=coordinates(state);
        for (size_type i=0; i<v.size(); ++i) {
          if(v[i]<-1 || v[i]>+1) {
            return false;
          }
        }
      
        return true;
      }
      
      /*! \brief Tests if the interior of the parallelopiped contains \a state. */
      inline bool interior_contains(const State& state) const {
        if (state.dimension()!=this->dimension()) {
         throw std::domain_error("This object and parameter have different space dimensions");
        } 
      
        if (this->empty()) { return false; }
        
        Vector v=coordinates(state);
        for (size_type i=0; i<v.size(); ++i) {
          if(v[i]<=-1 || v[i]>=+1) {
            return false;
          }
        }
      
        return true;
      }
      
      /*! \brief The equality operator (not implemented).
       *
       * Not currently implemented, since it requires matching the columns of 
       * the matrix of principal directions. 
       */
      inline bool operator==(const Zonotope<Real>& A) const
      {
        throw std::domain_error("Zonotope::operator==(...)  not implemented");
      }
      
      /*! \brief The inequality operator */
      inline bool operator!=(const Zonotope<Real>& A) const {
        throw std::domain_error("Zonotope::operator!=(...)  not implemented");
        return !(*this == A);
      }

      friend std::ostream&
      operator<< <> (std::ostream& os, 
                     const Zonotope<R>& r);
      
      friend std::istream&
      operator>> <> (std::istream& is, 
                     Zonotope<R>& r);
      
     private:
      // The linear inequalities defining the zonotope.
      inline void compute_linear_inequalities(Matrix&, Vector&, Vector&) const;
    };
  
    using namespace ::Ariadne::LinearAlgebra;
    
    template<typename R>
    inline 
    void Zonotope<R>::compute_linear_inequalities(Matrix& A, Vector& o, Vector& b) const
    {
      using namespace ::Ariadne::LinearAlgebra;
      size_type n=this->dimension();
      
      Vector c=this->centre() - State(n,0);
      //FIXME: can't use matrix inverse here
      A=inverse(this->principle_directions());
      o=A*c;
      b=vector<R>(n);
      for(size_type i=0; i!=n; ++i) {
        b[i]=1;
      }
    }
      


      
    template <typename R>
    Zonotope<R>::Zonotope(const Rectangle<Real>& r)
      : _centre(r.dimension()), _generators(r.dimension(),r.dimension())
    {
      for(size_type i=0; i!=dimension(); ++i) {
        _centre[i] = (r.lower_bound(i)+r.upper_bound(i))/2;
        _generators(i,i) = (r.upper_bound(i)-r.lower_bound(i))/2;
      }
    }
    
    template <typename R>
    Zonotope<R>::operator Polyhedron<R>() const 
    {
      using namespace ::Ariadne::LinearAlgebra;
      
      typedef typename Zonotope<R>::Real Real;
      typedef typename Zonotope<R>::State State;
      
      size_type n = this->dimension();
      
      /* Express in form invs * x - offst in [-bnds,+bnds] */
      matrix<R> invs;
      vector<R> offst;
      vector<R> bnds;
      this->compute_linear_inequalities(invs,offst,bnds);
      
     
      matrix<R> A(2*n,n);
      vector<R> b(2*n);
      
      for(uint i=0; i!=n; ++i) {
        for(uint j=0; j!=n; ++j) {
          A(i,j) = -invs(i,j);
          A(i+n,j) = invs(i,j);
        }
        b(i) = bnds(i)-offst(i);
        b(i+n) = bnds(i)+offst(i);
      }
      return Polyhedron<R>(A,b);
    }



    /*! \brief Tests disjointness */
    template <typename R>
    bool disjoint(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return !minkowski_difference(A,B).contains(State(A.dimension(),0))
      }
        
    /*! \brief Tests intersection of interiors */
    template <typename R>
    bool interiors_intersect(const Zonotope<R>& A,
                             const Zonotope<R>& B) 
    {
      return !minkowski_sum(A,B).interior_contains(State(A.dimension(),0));
    }
    
    
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    bool inner_subset(const Zonotope<R>& A,
                      const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }


    /*! \brief Tests inclusion in an open cover.  */
    template <typename R>
    bool subset_of_open_cover(const Zonotope<R>& A,
                              const ListSet<R, Zonotope >& cover) 
    {
      throw std::domain_error("subset_of_open_cover(Zonotope, std::vector<Zonotope>) not implemented");
    }

    
    /*! \brief Tests inclusion of \a A om the interior of \a B. */
    template <typename R>
    bool inner_subset(const Zonotope<R>& A,
                      const ListSet<R,Zonotope>& B) 
    {
      throw std::domain_error("subset_of_closed_cover(Zonotope, std::vector<Zonotope>) not implemented");
    }



    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Zonotope<R>& z) 
    {
      if(z.dimension() > 0) {
        os << "Zonotope(\n  centre=" << z.centre();
        os << "\n  directions=" << z.principle_directions();
        os << "\n) ";
      }

      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Zonotope<R>& z)
    {
      throw std::domain_error("Not implemented");
    }
      
    

  }
}

#endif /* _ARIADNE_ZONOTOPE_H */
