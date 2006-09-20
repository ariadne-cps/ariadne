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
#include "../geometry/point_list.h"
#include "../geometry/polytope.h"

namespace Ariadne {
  namespace Geometry {

    template<> 
    inline bool is_a<Simplex,Simplex>() { return true; }
    template<> 
    inline bool is_a<Simplex,Polytope>() { return true; }

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Simplex<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Simplex<R>&);

    /*! \ingroup BasicSet
     *  \brief A simplex of arbitrary dimension.
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
      LinearAlgebra::Matrix<R> _vertices;
   
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor constructs standard simplex of dimension \a n. */
      Simplex(size_type n = 0);
    
      /*! \brief Construct from matrix giving the vertices in column form. */
      explicit Simplex(const LinearAlgebra::Matrix<R>& A);
      
      /*! \brief Construct from list of vertices. */
      explicit Simplex(const PointList<R>& v);
      
      /*! \brief Construct from a string literal. */
      explicit Simplex(const std::string& s);
      
      /*! \brief Copy constructor. */
      Simplex(const Simplex<R>& original)
        : _vertices(original._vertices)
      { }

      /*! \brief Construct from a rectangle. 
       *
       * It returns an error since, it can not be doen without errors.
       */
      Simplex(const Rectangle<R>& r){ 
        throw std::runtime_error("Simplex(const Rectangle<R>& r): errorless conversion can not be accomplished.");
      }
      
      /*! \brief Copy assignment operator. */
      Simplex<R>& operator=(const Simplex<R>& original) {
        if(this != &original) {
          this->_vertices = original._vertices;
        }
        return *this;
      }
      //@}
      
    
      //@{
      //! \name Conversion operators
      /*! \brief Convert to a Parma Polyhedra Library polyhedron. */
      operator Parma_Polyhedra_Library::C_Polyhedron() const;

      /*! \brief Convert to a polytope. */
      operator Polytope<R>() const;
      //@}
      
      
      //@{
      //! \name Geometric operations
      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const {
        return this->_vertices.size1();
      }
      
      /*! \brief True if the simplex is empty. */
      bool empty() const {
        throw std::domain_error("Simplex::empty() not implemented.");
      }
      
      /*! \brief True if the simplex has empty interior. */
      bool empty_interior() const {
        throw std::domain_error("Simplex::empty_interior() not implemented.");
      }
      
      /*! \brief The array of vertices. */
      PointList<R> vertices() const {
        return PointList<R>(this->_vertices);
      }
      
      /*! \brief The \a n th vertex. */
      state_type vertex(size_type n) const {
        return Point<R>(LinearAlgebra::Vector<R>(column(this->_vertices,n)));
      }
      
      /*! \brief Tests if \a point is included into a simplex. */
      bool contains(const state_type& point) const;
      
      /*! \brief Tests if \a point is included into the interior a simplex. */
      bool interior_contains(const state_type& point) const;
      //@}

      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
    };


    template<typename R> inline 
    std::ostream& operator<<(std::ostream& os, const Simplex<R>& s) {
      return s.write(os);
    }
    
    template<typename R> inline
    std::istream& operator>>(std::istream& is, Simplex<R>& s) {
      return s.read(is);
    }

     
  }
}

#endif /* _ARIADNE_SIMPLEX_H */
