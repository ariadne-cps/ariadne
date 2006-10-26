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
    template<class R> std::ostream& operator<<(std::ostream&, const Simplex<R>&);
    template<class R> std::istream& operator>>(std::istream&, Simplex<R>&);

    /*! \ingroup BasicSet
     *  \brief A simplex of arbitrary dimension.
     */
    template<class R>
    class Simplex : public Polytope<R>
    {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the simplex. */
      typedef Point<R> state_type;
     
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor constructs simplex of dimension \a 0. */
      Simplex();
    
      /*! \brief Construct from matrix giving the vertices in column form. */
      explicit Simplex(const LinearAlgebra::Matrix<R>& A);
      
      /*! \brief Construct from list of vertices. */
      explicit Simplex(const PointList<R>& v);
      
      /*! \brief Copy constructor. */
      Simplex(const Simplex<R>& s) : Polytope<R>(s) { }

      /*! \brief Copy assignment operator. */
      Simplex<R>& operator=(const Simplex<R>& s) {
        if(this != &s) {
          this->Polytope<R>::operator=(s); 
        }
        return *this;
      }
      //@}
    
      //! \name Geometric predicates
      //! \brief Specialized containment predicate
      tribool contains(const Point<R>& pt) const;
      //@}
      
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     public:
      LinearAlgebra::Vector<typename Numeric::traits<R>::arithmetic_type> coordinates(const Point<R>& s) const;

    };

  
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const Simplex<R>& s) {
      return s.write(os);
    }
    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Simplex<R>& s) {
      return s.read(is);
    }

     
  }
}

#endif /* _ARIADNE_SIMPLEX_H */
