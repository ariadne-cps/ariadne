/***************************************************************************
 *            sphere.h
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
 
/*! \file sphere.h
 *  \brief Solid spheres in Euclidean space.
 */

#ifndef _ARIADNE_SPHERE_H
#define _ARIADNE_SPHERE_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/array.h"
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    /* Forward declaration of friends. */
    template<class R> std::ostream& operator<<(std::ostream&, const Sphere<R>&);
    template<class R> std::istream& operator>>(std::istream&, Sphere<R>&);

    template<class R> R euclidean_distance_square_down(const Point<R>&,const Point<R>&);
    template<class R> R euclidean_distance_square_up(const Point<R>&,const Point<R>&);
    
    /*! \ingroup BasicSet
     *  \brief A ball \f$||x-c||\leq r\f$ of arbitrary dimension.
     */
    template <class R>
    class Sphere {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the simplex. */
      typedef Point<R> state_type;
     
     private:
      /* Simplex's vertices. */
      state_type _centre;
      real_type _radius;
     
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Default constructor constructs standard simplex of dimension \a n. */
      Sphere(size_type n = 0);
    
      /*! \brief Construct from centre and radius. */
      explicit Sphere(const state_type& c, const real_type& r);
     
      /*! \brief Construct from centre and radius. */
      explicit Sphere(const state_type& c, real_type r);
     
      /*! \brief Construct from a string literal. */
      explicit Sphere(const std::string& s);
      
      /*! \brief Copy constructor. */
      Sphere(const Sphere<R>& original)
        : _centre(original._centre), _radius(original._radius)
      { }
      
      /*! \brief Copy assignment operator. */
      Sphere<R>& operator=(const Sphere<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_radius = original._radius;
        }
        return *this;
      }
      //@}
      
      
      //@{ 
      //! \name Comparison operators
      /*! \brief The equality operator.
       */
      bool operator==(const Sphere<R>& other) const
      {
        return this->_centre==other._centre && this->_radius==other._radius;
      }
      
      /*! \brief The inequality operator */
      bool operator!=(const Sphere<R>& other) const {
        return !(*this == other);
      }
      //@}
      
      
      //@{ 
      //! \name Data elements
      /*! \brief The centre of the ball. */
      const state_type& centre() const {
        return this->_centre;
      }
      
      /*! \brief The radius of the ball. */
      const real_type& radius() const {
        return this->_radius;
      }
      //@}
      
      
      //@{ 
      //! \name Geometric operations
      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const {
        return this->_centre.dimension();
      }
      
      /*! \brief True if the ball is empty. */
      bool empty() const {
        return this->_radius<R(0);
      }
      
      /*! \brief True if the simplex has empty interior. */
      bool empty_interior() const {
        return this->_radius <= R(0);
      }
      
      /*! \brief Tests if \a point is included into a simplex. */
      bool contains(const state_type& point) const {
        return euclidean_distance_square_down(point,this->_centre) <= mul_approx(this->_radius , this->_radius);
      }
      
      /*! \brief Tests if \a point is included into the interior a simplex. */
      bool interior_contains(const state_type& point) const {
        return euclidean_distance_square_up(point,this->_centre) < mul_approx(this->_radius , this->_radius);
      }
      //@}
      
#ifdef DOXYGEN
    //@{
    //! \name Geometric binary predicates
    /*! \brief Tests disjointness */
    friend bool disjoint(const Sphere<R>& A, const Sphere<R>& B);
    /*! \brief Tests disjointness */
    friend bool disjoint(const Rectangle<R>& A, const Sphere<R>& B);
    /*! \brief Tests disjointness */
    friend bool disjoint(const Sphere<R>& A, const Rectangle<R>& B);
    /*! \brief Tests intersection of interiors */
    friend bool interiors_intersect(const Sphere<R>& A, const Sphere<R>& B);
    /*! \brief Tests intersection of interiors */
    friend bool interiors_intersect(const Rectangle<R>& A, const Sphere<R>& B);
    /*! \brief Tests intersection of interiors */
    friend bool interiors_intersect(const Sphere<R>& A, const Rectangle<R>& B);
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    friend bool inner_subset(const Sphere<R>& A, const Sphere<R>& B);
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    friend bool inner_subset(const Rectangle<R>& A, const Sphere<R>& B);
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    friend bool inner_subset(const Sphere<R>& A, const Rectangle<R>& B);
    /*! \brief Tests inclusion of \a A in \a B. */
    friend bool subset(const Sphere<R>& A, const Sphere<R>& B);
    /*! \brief Tests inclusion of \a A in \a B. */
    friend bool subset(const Rectangle<R>& A, const Sphere<R>& B);
    /*! \brief Tests inclusion of \a A in \a B. */
    friend bool subset(const Sphere<R>& A, const Rectangle<R>& B);
    //@}
    
    //@{
    //! \name Geometric binary operations
    /*! \brief The Minkoswi sum of two spheres */
    friend Sphere<R> minkowski_sum(const Sphere<R>& A, const Sphere<R>& B);
   /*! \brief The Minkoswi difference of two spheres */
    friend Sphere<R> minkowski_difference(const Sphere<R>& A, const Sphere<R>& B);
    //@}
#endif      
      friend std::ostream& operator<< <> (std::ostream& os, const Sphere<R>& r);
      friend std::istream& operator>> <> (std::istream& is, Sphere<R>& r);
    };
          
    template<class R>
    inline
    R
    euclidean_distance_square_down(const Point<R>& pt1, const Point<R>& pt2) 
    {
      R result=R(0);
      R tmp;
      for(dimension_type i=0; i!=pt1.dimension(); ++i) {
        tmp = (pt1[i]>=pt2[i]) ? sub_down(pt1[i],pt2[i]) : sub_down(pt2[i],pt1[i]);
        result=add_down(result,mul_down(tmp,tmp));
      }
      return result;
    }
      
    template<class R>
    inline
    R
    euclidean_distance_square_up(const Point<R>& pt1, const Point<R>& pt2) 
    {
      R result=R(0);
      R tmp;
      for(dimension_type i=0; i!=pt1.dimension(); ++i) {
        tmp = (pt1[i]>=pt2[i]) ? sub_up(pt1[i],pt2[i]) : sub_up(pt2[i],pt1[i]);
        result=add_up(result,mul_up(tmp,tmp));
      }
      return result;
    }
      
    template <class R>
    inline bool disjoint(const Sphere<R>& A, const Sphere<R>& B) 
    {
      return euclidean_distance_down(A.centre(),B.centre()) > 
        pow_up(add_up(A.radius(),B.radius()),2);
    }
    
    template <class R>
    inline bool disjoint(const Sphere<R>& A, const Rectangle<R>& B) 
    {
      throw std::runtime_error("bool disjoint(const Sphere<R>&, const Rectangle<R>&) not implemented");
    }
    
    template <class R>
    inline bool disjoint(const Rectangle<R>& A, const Sphere<R>& B) 
    {
      return disjoint(B,A);
    }
    
    
    template <class R>
    inline bool interiors_intersect(const Sphere<R>& A,
                                    const Sphere<R>& B) 
    {
      return euclidean_distance_up(A.centre(),B.centre()) < 
        pow_down(add_down(A.radius(),B.radius()),2);
    }
    
    template <class R>
    inline bool interiors_intersect(const Sphere<R>& A,
                                    const Rectangle<R>& B) 
    {
      throw std::runtime_error("bool disjoint(const Sphere<R>&, const Rectangle<R>&) not implemented");
    }
    
    template <class R>
    inline bool interiors_intersect(const Rectangle<R>& A,
                                    const Sphere<R>& B) 
    {
      return interiors_intersect(B,A);
    }
    
    
    template <class R>
    inline bool inner_subset(const Sphere<R>& A,
                             const Sphere<R>& B) 
    {
      throw std::runtime_error("bool inner_subset(const Sphere<R>&, const Sphere<R>&) not implemented");
      //return A.radius()<B.radius && euclidean_distance_square(A.centre()-B.centre()) < square(B.centre()-A.centre());
    }

    template <class R>
    inline bool inner_subset(const Sphere<R>& A,
                             const Rectangle<R>& B) 
    {
      throw std::runtime_error("bool inner_subset(const Sphere<R>&, const Rectangle<R>&) not implemented");
      //for(dimension_type i=0; i!=A.dimension(); ++i) {
      //  if(! inner_subset(A.centre()[i]-A.radius(),A.centre()[i]+A.radius(),B[i]) ) {
      //    return false;
      //  }
      //}
      return true;
    }

    template <class R>
    inline bool inner_subset(const Rectangle<R>& A,
                             const Sphere<R>& B) 
    {
      throw std::runtime_error("bool inner_subset(const Rectangle<R>&, const Sphere<R>&) not implemented");
      //array< Point<R> > vertices=A.vertices();
      //for(class Rectangle<R>::vertex_iterator vertex_iter=vertices.begin(); vertex_iter!=vertices.end(); ++vertex_iter) {
      //  if(! B.interior_contains(*vertex_iter) ) {
      //    return false;
      //  }
      //}
      return true;
    }

    
    template <class R>
    inline bool subset(const Sphere<R>& A, 
                       const Sphere<R>& B) 
    {
      throw std::runtime_error("bool subset(const Sphere<R>&, const Sphere<R>&) not implemented");
      //return A.radius()<=B.radius && euclidean_distance_square(A.centre()-B.centre()) <= square(B.centre()-A.centre());
    }
    
    template <class R>
    inline bool subset(const Sphere<R>& A, 
                       const Rectangle<R>& B) 
    {
      throw std::runtime_error("bool subset(const Sphere<R>&, const Rectangle<R>&) not implemented");
      //return subset(A.bounding_box(),B);
    }
    
    template <class R>
    inline bool subset(const Rectangle<R>& A, 
                       const Sphere<R>& B) 
    {
      throw std::runtime_error("bool subset(const Rectangle<R>&, const Sphere<R>&) not implemented");
      //array< Point<R> > vertices=A.vertices();
      //for(class Rectangle<R>::vertex_iterator vertex_iter=vertices.begin(); vertex_iter!=vertices.end(); ++vertex_iter) {
      //  if(! B.contains(*vertex_iter) ) {
      //    return false;
      //  }
      //}
      return true;
    }

    template<class R>
    inline
    Geometry::Sphere<R> 
    scale(const Geometry::Sphere<R>& s, const R& scale_factor) {

      const Geometry::Point<R>& centre=s.centre();
      Geometry::Point<R> new_centre(s.dimension());

      for(size_type i=0; i!=s.dimension(); ++i) {
        new_centre[i]=scale_factor*centre[i];
      }

      return Geometry::Sphere<R>(new_centre, scale_factor*s.radius());
    }

    template <class R>
    std::ostream&
    operator<<(std::ostream& os, const Sphere<R>& s); 

    
    template <class R>
    std::istream& 
    operator>>(std::istream& is, Sphere<R>& s);

      
  }
}

#endif /* _ARIADNE_SPHERE_H */
