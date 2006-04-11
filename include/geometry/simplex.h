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

#include "../declarations.h"

#include "../base/array.h"
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Simplex<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Simplex<R>&);

    /*! \brief A simplex of arbitrary dimension.
     */
    template <typename R>
    class Simplex {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the simplex. */
      typedef Point<R> state_type;
     
     private:
      /* Simplex's vertices. */
      array<state_type> _vertices;
   
     public:
      /*! \brief Default constructor constructs standard simplex of dimension \a n. */
      Simplex(size_type n = 0);
    
      /*! \brief Construct from list of vertices. */
      explicit Simplex(const std::vector<state_type>& v);
     
      /*! \brief Construct from list of vertices. */
      explicit Simplex(const array<state_type>& v);
      
      /*! \brief Construct from a string literal. */
      explicit Simplex(const std::string& s);
      
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
      
      /*! \brief The equality operator (not implemented).
       *
       * Not currently implemented, since it requires matching the columns of 
       * the matrix of principal directions. 
       */
      bool operator==(const Simplex<real_type>& A) const
      {
        throw std::domain_error("Simplex::operator==(...)  not implemented");
      }
      
      /*! \brief The inequality operator */
      bool operator!=(const Simplex<real_type>& A) const {
        throw std::domain_error("Simplex::operator!=(...)  not implemented");
        return !(*this == A);
      }
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<R> () const;
      
      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const {
        return (this->_vertices)[0].dimension();
      }
      
      /*! \brief True if the rectangle is empty. */
      bool empty() const {
        throw std::domain_error("Simplex::empty() not implemented.");
      }
      
      /*! \brief True if the rectangle has empty interior. */
      bool empty_interior() const {
        throw std::domain_error("Simplex::empty_interior() not implemented.");
      }
      
      /*! \brief The array of vertices. */
      const array<state_type>& vertices() const {
        return this->_vertices;
      }
      
      /*! \brief The \a n th vertex. */
      const state_type& vertex(size_type n) const {
        return this->_vertices[n];
      }
      
      /*! \brief Tests if \a point is included into a simplex. */
      bool contains(const state_type& point) const {
        return Polyhedron<R>(*this).contains(point);
      }
      
      /*! \brief Tests if \a point is included into the interior a simplex. */
      bool interior_contains(const state_type& point) const {
        return Polyhedron<R>(*this).interior_contains(point);
      }
    
      friend std::ostream& operator<< <> (std::ostream& os, const Simplex<R>& r);
      friend std::istream& operator>> <> (std::istream& is, Simplex<R>& r);
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
    inline bool subset(const Simplex<R>& A, 
                const Simplex<R>& B) 
    {
      return subset(Polyhedron<R>(A),Polyhedron<R>(B));
    }
    
    

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Simplex<R>& s); 

    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Simplex<R>& s);

      
  }
}

#endif /* _ARIADNE_SIMPLEX_H */
