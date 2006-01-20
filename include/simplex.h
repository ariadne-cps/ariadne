/***************************************************************************
 *            simplex.h
 *
 *  6 January 2006
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
 
/*! \file simplex.h
 *  \brief Simplices.
 */

#ifndef _ARIADNE_SIMPLEX_H
#define _ARIADNE_SIMPLEX_H

#include <iosfwd>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "ariadne.h"
#include "array.h"
#include "utility.h"
#include "linear_algebra.h"
#include "interval.h"
#include "point.h"
#include "rectangle.h"
#include "polyhedron.h"
#include "list_set.h"

namespace Ariadne {
  namespace Geometry {
    template < typename R > class Simplex;
    template < typename R, template <typename> class BS > class ListSet;

    template <typename R> Simplex<R> intersection(const Simplex<R> &A, const Simplex<R> &B);
    template <typename R> Simplex<R> regular_intersection(const Simplex<R> &A, const Simplex<R> &B);

    template<typename R> bool interiors_intersect(const Simplex<R> &A, const Simplex<R> &B);
    template<typename R> bool disjoint(const Simplex<R> &A, const Simplex<R> &B);
    template<typename R> bool inner_subset(const Simplex<R> &A, const Simplex<R> &B);
    template<typename R> bool subset(const Simplex<R> &A, const Simplex<R> &B);
    template<typename R> bool subset_of_open_cover(const Simplex<R> &A, const std::vector< Simplex<R> > &list);
    template<typename R> bool subset_of_closed_cover(const Simplex<R> &A, const std::vector< Simplex<R> > &list);

    template<typename R> bool inner_subset(const Simplex<R> &rect, const ListSet<R,Simplex> &A);
    template<typename R> bool subset(const Simplex<R> &rect, const ListSet<R,Simplex> &A);
    
    template<typename R> std::ostream& operator<<(std::ostream&, const Simplex<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Simplex<R>&);

    /*! \brief A simplex of arbitrary dimension.
     */
    template <typename R>
    class Simplex {
      /*! \brief Makes intersection */
      friend Simplex<R> intersection <> (const Simplex<R>& A,
                                         const Simplex<R>& B);

      /*! \brief Makes intersection of interiors */
      friend Simplex<R> regular_intersection <> (const Simplex<R>& A,
                                                   const Simplex<R>& B);

       /*! \brief Tests intersection of interiors. */
      friend bool interiors_intersect <> (const Simplex<R>& A,
                                          const Simplex<R>& B);

       /*! \brief Tests disjointness */
      friend bool disjoint <> (const Simplex<R>& A,
                               const Simplex<R>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a B. */
      friend bool inner_subset <> (const Simplex<R>& A,
                                   const Simplex<R>& B);

      /*! \brief Tests inclusion. */
      friend bool subset <> (const Simplex<R>& A,
                             const Simplex<R>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a DS. */
      friend bool inner_subset <> (const Simplex<R>& A,
                                   const ListSet<R,Ariadne::Geometry::Simplex>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a DS. */
      friend bool subset <> (const Simplex<R>& A,
                             const ListSet<R,::Ariadne::Geometry::Simplex>& B);


      /*! \brief Tests inclusion in an open cover.
       *  \internal We shouldn't restrict to a std::list.
       */
      friend bool subset_of_open_cover <> (const Simplex<R>& A,
                                           const std::vector< Simplex<R> >& list);

      /*! \brief Tests inclusion in an closed cover.
       *  \internal We shouldn't restrict to a std::list.
       */
      friend bool subset_of_closed_cover <> (const Simplex<R>& A,
                                             const std::vector< Simplex<R> >& list);

     public:
      /*! \brief The unsigned integer type used to denote the array positions. */
      typedef size_t size_type;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R Real;
      /*! \brief The type of denotable state contained by the simplex. */
      typedef Point<R> State;
     private:
      /* Simplex's vertices. */
      array<State> _vertices;
   
     public:
      /*! \brief Default constructor constructs standard simplex of dimension \a n. */
      Simplex(size_type n = 0)
        : _vertices(n+1,State(n)) 
      {
        for(size_type i=0; i!=n; ++i) {
          _vertices[i][i]=1;
        }
      }
    
      /*! \brief Construct from list of vertices. */
      Simplex(const std::vector<State>& v)
        : _vertices(v.begin(),v.end())
      {
        size_type d=_vertices.size()-1;
        for(size_type i=0; i!=d+1; ++i) {
          if(_vertices[i].dimension()!=d) {
            throw std::domain_error("The the list of vertices is invalid");
          }
        }
      }
      
      /*! \brief Construct from list of vertices. */
      Simplex(const array<State>& v)
        : _vertices(v.begin(),v.end())
      {
        size_type d=_vertices.size()-1;
        for(size_type i=0; i!=d+1; ++i) {
          if(_vertices[i].dimension()!=d) {
            throw std::domain_error("The the list of vertices is invalid");
          }
        }
      }
      
      /*! \brief Construct from a string literal. */
      explicit Simplex(const std::string& s)
        : _vertices()
      {
        std::stringstream ss(s);
        ss >> *this;
      }
      
      /*! \brief Copy constructor. */
      Simplex(const Simplex<R>& original)
        : _vertices(original._vertices)
      { }
      
      /*! \brief Copy assignment operator. */
      Simplex<R>& operator=(const Simplex<R>& original) {
        if(this != &original) {
          this->_vertices = original._vertices;
        }
        return *this;
      }
      
      operator Polyhedron<R> () const {
        std::vector<State> vert_vec(_vertices.begin(),_vertices.end());
        return Polyhedron<R>(vert_vec);
      }
      
      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      inline size_type dimension() const {
        return (this->_vertices)[0].dimension();
      }
      
      /*! \brief True if the rectangle is empty. */
      inline bool empty() const {
        throw std::domain_error("Simplex::empty() not implemented.");
      }
      
      /*! \brief True if the rectangle has empty interior. */
      inline bool empty_interior() const {
        throw std::domain_error("Simplex::empty_interior() not implemented.");
      }
      
      /*! \brief The array of vertices. */
      inline const array<State>& vertices() const {
        return this->_vertices;
      }
      
      /*! \brief The @a n th vertex. */
      inline const State& vertex(size_type n) const {
        return this->_vertices[n];
      }
      
      /*! \brief Tests if \a state is included into a simplex. */
      inline bool contains(const State& state) const {
        return Polyhedron<R>(*this).contains(state);
      }
      
      /*! \brief Tests if \a state is included into the interior a simplex. */
      inline bool interior_contains(const State& state) const {
        return Polyhedron<R>(*this).interior_contains(state);
      }
    
      /*! \brief The equality operator (not implemented).
       *
       * Not currently implemented, since it requires matching the columns of 
       * the matrix of principal directions. 
       */
      inline bool operator==(const Simplex<Real>& A) const
      {
        throw std::domain_error("Simplex::operator==(...)  not implemented");
      }
      
      /*! \brief The inequality operator */
      inline bool operator!=(const Simplex<Real>& A) const {
        throw std::domain_error("Simplex::operator!=(...)  not implemented");
        return !(*this == A);
      }
      
      friend std::ostream&
      operator<< <> (std::ostream& os, 
                     const Simplex<R>& r);
      
      friend std::istream&
      operator>> <> (std::istream& is, 
                     Simplex<R>& r);
      
    };
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline bool disjoint(const Simplex<R>& A, const Simplex<R>& B) 
    {
      return disjoint(Polyhedron<R>(A),Polyhedron<R>(B));
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline bool disjoint(const Simplex<R>& A, const Rectangle<R>& B) 
    {
      return disjoint(Polyhedron<R>(A),Polyhedron<R>(B));
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline bool disjoint(const Rectangle<R>& A, const Simplex<R>& B) 
    {
      return disjoint(B,A);
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline bool disjoint(const Simplex<R>& A, const Polyhedron<R>& B) 
    {
      return disjoint(Polyhedron<R>(A),B);
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline bool disjoint(const Polyhedron<R>& A, const Simplex<R>& B) 
    {
      return disjoint(B,A);
    }
    
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline bool interiors_intersect(const Simplex<R>& A,
                                    const Simplex<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A),Polyhedron<R>(B));
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline bool interiors_intersect(const Simplex<R>& A,
                                    const Rectangle<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A),Polyhedron<R>(B));
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline bool interiors_intersect(const Rectangle<R>& A,
                                    const Simplex<R>& B) 
    {
      return interiors_intersect(B,A);
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline bool interiors_intersect(const Simplex<R>& A,
                                    const Polyhedron<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A),B);
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline bool interiors_intersect(const Polyhedron<R>& A,
                                    const Simplex<R>& B) 
    {
      return interiors_intersect(B,A);
    }
    
    
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline bool inner_subset(const Simplex<R>& A,
                             const Simplex<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A),Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline bool inner_subset(const Simplex<R>& A,
                             const Rectangle<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A),Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline bool inner_subset(const Rectangle<R>& A,
                             const Simplex<R>& B) 
    {
      return inner_subset(B,A);
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline bool inner_subset(const Simplex<R>& A,
                             const Polyhedron<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A),B);
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline bool inner_subset(const Polyhedron<R>& A,
                             const Simplex<R>& B) 
    {
      return inner_subset(B,A);
    }

    
    /*! \brief Tests inclusion of \a A in \a B. */
    template <typename R>
    bool subset(const Simplex<R>& A, 
                const Simplex<R>& B) 
    {
      return subset(Polyhedron<R>(A),Polyhedron<R>(B));
    }
    

    //*! \brief Tests inclusion in an open cover.  */
    template <typename R>
    bool subset_of_open_cover(const Simplex<R>& A,
                              const std::vector< Simplex<R> >& cover) 
    {
      throw std::domain_error("subset_of_open_cover(Simplex, std::vector<Simplex>) not implemented");
    }

    
    //*! \brief Tests inclusion in a closed cover.  */
    template <typename R>
    bool subset_of_closed_cover(const Simplex<R>& A,
                                const std::vector< Simplex<R> >& cover) 
    {
      throw std::domain_error("subset_of_closed_cover(Simplex, std::vector<Simplex>) not implemented");
    }


    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Simplex<R>& s) 
    {
//      if(p.empty()) {
//        os << "Empty";
//     }
//      else 
      if(s.dimension() > 0) {
        os << "Simplex(vertices=" << s.vertices() << ") ";
      }

      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Simplex<R>& s)
    {
      throw std::domain_error("Not implemented");
    }
      
  }
}

#endif /* _ARIADNE_SIMPLEX_H */
