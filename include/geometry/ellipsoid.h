/***************************************************************************
 *            ellipsoid.h
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
 
/*! \file ellipsoid.h
 *  \brief Solid ellipsoids in Euclidean space.
 */

#ifndef _ARIADNE_ELLIPSOID_H
#define _ARIADNE_ELLIPSOID_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/array.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Ellipsoid<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Ellipsoid<R>&);

    template<typename R> R euclidean_norm_square(const LinearAlgebra::Vector<R>&);
    
    /*! \ingroup BasicSet
     *  \brief An ellipsoid \f$(x-c)^T A (x-c)\leq 1\f$ of arbitrary dimension.
     */
    template <typename R>
    class Ellipsoid {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the simplex. */
      typedef Point<R> state_type;
      /*! \brief The type of vector used to describe differences between points. */
      typedef LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix used to describe lengths and directions of the axes. */
      typedef LinearAlgebra::Matrix<R> matrix_type;
     
     private:
      state_type _centre;
      matrix_type _bilinear_form;
     
     public:
      //! \name Constructors 
      //@{ 
      /*! \brief Default constructor constructs unit sphere centred at the origin. \a n. */
      Ellipsoid(size_type n = 0);
    
      /*! \brief Construct from centre and bilinear form. */
      explicit Ellipsoid(const state_type& c, const matrix_type& r);
     
      /*! \brief Construct from a string literal. */
      explicit Ellipsoid(const std::string& s);
      
      /*! \brief Conversion constructor from a sphere. */
      /* FIXME: Conversion only works if R is a field; otherwise division need not be exact */
      Ellipsoid(const Sphere<R>& s);
      
      /*! \brief Copy constructor. */
      Ellipsoid(const Ellipsoid<R>& original)
        : _centre(original._centre), _bilinear_form(original._bilinear_form)
      { }
      
      /*! \brief Copy assignment operator. */
      Ellipsoid<R>& operator=(const Ellipsoid<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_bilinear_form = original._bilinear_form;
        }
        return *this;
      }
      //@}
      
      //@{
      //! \brief Comparison operators
      /*! \brief The equality operator (not implemented).
       *
       * Not currently implemented, since it requires matching the columns of 
       * the Matrix of principal directions. 
       */
      bool operator==(const Ellipsoid<R>& other) const
      {
        return this->_centre==other._centre && this->_bilinear_form==other._bilinear_form;
      }
      
      /*! \brief The inequality operator. */
      bool operator!=(const Ellipsoid<R>& other) const {
        return !(*this == other);
      }
      //@}
      
      //! \name Data elements
      //@{
      /*! \brief The centre of the ellipsoid. */
      const state_type& centre() const {
        return this->_centre;
      }
      
      /*! \brief The bilinear form describing axes of the ellipsoid. */
      const matrix_type& bilinear_form() const {
        return this->_bilinear_form;
      }
      //@}
      
      //! \name Geometric Operators
      //@{
      /*! \brief The dimension of the Euclidean space the ellipsoid lies in. */
      size_type dimension() const {
        return this->_centre.dimension();
      }
      
      /*! \brief True if the ellipsoid is empty. */
      bool empty() const {
        return false;
      }
      
      /*! \brief True if the ellipsoid has empty interior. */
      bool empty_interior() const {
        return _bilinear_form.singular();
      }
      
      /*! \brief A rectangle containing the ellipsoid. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief Tests if \a point is contained in the ellipsoid. */
      bool contains(const state_type& point) const {
        vector_type vec=point - this->_centre;
        return inner_product(vec,vector_type(this->_bilinear_form*vec))<=1;
      }
      
      /*! \brief Tests if \a point is contained in the interior of the ellipsoid. */
      bool interior_contains(const state_type& point) const {
        vector_type vec=point-this->_centre;
        return inner_product(vec,this->_bilinear_form*vec)<1;
      }
      //@}
      
#ifdef DOXYGEN
    //@{
    //! \name Geometric binary predicates
    /*! \brief Tests disjointness */
    friend bool disjoint(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
    friend bool disjoint(const Rectangle<R>& A, const Ellipsoid<R>& B);
    friend bool disjoint(const Ellipsoid<R>& A, const Rectangle<R>& B);
    /*! \brief Tests intersection of interiors */
    friend bool interiors_intersect(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
    friend bool interiors_intersect(const Rectangle<R>& A, const Ellipsoid<R>& B);
    friend bool interiors_intersect(const Ellipsoid<R>& A, const Rectangle<R>& B);
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    friend bool inner_subset(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
    friend bool inner_subset(const Rectangle<R>& A, const Ellipsoid<R>& B);
    friend bool inner_subset(const Ellipsoid<R>& A, const Rectangle<R>& B);
    /*! \brief Tests inclusion of \a A in \a B. */
    friend bool subset(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
    friend bool subset(const Rectangle<R>& A, const Ellipsoid<R>& B);
    friend bool subset(const Ellipsoid<R>& A, const Rectangle<R>& B);
    //@}
    
    //@{
    //! \name Geometric binary operations
    /*! \brief The Minkoswi sum of two ellipsoids */
    friend Ellipsoid<R> minkowski_sum(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
    /*! \brief The Minkoswi difference of two ellipsoids */
    friend Ellipsoid<R> minkowski_difference(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
    //@}
#endif
      friend std::ostream& operator<< <> (std::ostream& os, const Ellipsoid<R>& r);
      friend std::istream& operator>> <> (std::istream& is, Ellipsoid<R>& r);
    };
    
     
    template <typename R>
    inline bool disjoint(const Ellipsoid<R>& A, const Ellipsoid<R>& B) 
    {
      throw std::runtime_error("bool disjoint(const Ellipsoid<R>&, const Ellipsoid<R>&) not implemented");
    }
    
    template <typename R>
    inline bool disjoint(const Ellipsoid<R>& A, const Rectangle<R>& B) 
    {
      throw std::runtime_error("bool disjoint(const Ellipsoid<R>&, const Rectangle<R>&) not implemented");
    }
    
    template <typename R>
    inline bool disjoint(const Rectangle<R>& A, const Ellipsoid<R>& B) 
    {
      return disjoint(B,A);
    }
    
    
    template <typename R>
    inline bool interiors_intersect(const Ellipsoid<R>& A,
                                    const Ellipsoid<R>& B) 
    {
      throw std::runtime_error("bool interiors_intersect(const Ellipsoid<R>&, const Ellipsoid<R>&) not implemented");
    }
    
    template <typename R>
    inline bool interiors_intersect(const Ellipsoid<R>& A,
                                    const Rectangle<R>& B) 
    {
      throw std::runtime_error("bool disjoint(const Ellipsoid<R>&, const Rectangle<R>&) not implemented");
    }
    
    template <typename R>
    inline bool interiors_intersect(const Rectangle<R>& A,
                                    const Ellipsoid<R>& B) 
    {
      return interiors_intersect(B,A);
    }
    
    
    template <typename R>
    inline bool inner_subset(const Ellipsoid<R>& A,
                             const Ellipsoid<R>& B) 
    {
      throw std::runtime_error("bool inner_subset(const Ellipsoid<R>&, const Ellipsoid<R>&) not implemented");
    }

    template <typename R>
    inline bool inner_subset(const Ellipsoid<R>& A,
                             const Rectangle<R>& B) 
    {
      for(dimension_type i=0; i!=A.dimension(); ++i) {
        if(! inner_subset(A.centre()[i]-A.radius(),A.centre()[i]+A.radius(),B[i]) ) {
          return false;
        }
      }
      return true;
    }

    template <typename R>
    inline bool inner_subset(const Rectangle<R>& A,
                             const Ellipsoid<R>& B) 
    {
      array< Point<R> > vertices=A.vertices();
      for(typename Rectangle<R>::vertex_iterator vertex_iter=vertices.begin(); vertex_iter!=vertices.end(); ++vertex_iter) {
        if(! B.interior_contains(*vertex_iter) ) {
          return false;
        }
      }
      return true;
    }

    
    template <typename R>
    inline bool subset(const Ellipsoid<R>& A, 
                       const Ellipsoid<R>& B) 
    {
      throw std::runtime_error("bool subset(const Ellipsoid<R>&, const Ellipsoid<R>&) not implemented");
    }
    
    template <typename R>
    inline bool subset(const Ellipsoid<R>& A, 
                       const Rectangle<R>& B) 
    {
      return subset(A.bounding_box(),B);
    }
    
    template <typename R>
    inline bool subset(const Rectangle<R>& A, 
                       const Ellipsoid<R>& B) 
    {
      array< Point<R> > vertices=A.vertices();
      for(typename Rectangle<R>::vertex_iterator vertex_iter=vertices.begin(); vertex_iter!=vertices.end(); ++vertex_iter) {
        if(! B.contains(*vertex_iter) ) {
          return false;
        }
      }
      return true;
    }

    template<typename R>
    inline
    Geometry::Ellipsoid<R> 
    scale(const Geometry::Ellipsoid<R>& s, const R& scale_factor) {

      const Geometry::Point<R>& centre=s.centre();
      const LinearAlgebra::Matrix<R>& bilinear_form=s.bilinear_form();
      
      Geometry::Point<R> new_centre(s.dimension());

      for(size_type i=0; i!=s.dimension(); ++i) {
        new_centre[i]=scale_factor*centre[i];
      }

      return Geometry::Ellipsoid<R>(new_centre, scale_factor*bilinear_form);
    }

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Ellipsoid<R>& s); 

    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Ellipsoid<R>& s);

      
  }
}

#endif /* _ARIADNE_SPHERE_H */
