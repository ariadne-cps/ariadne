/***************************************************************************
 *            transformation_system.h
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
 
/*! \file transformation_system.h
 *  \brief Zonotopic Vectors (affine images of cuboids).
 */

#ifndef _ARIADNE_ZONOTOPIC_VECTOR_H
#define _ARIADNE_ZONOTOPIC_VECTOR_H

#include <iosfwd>

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/interval_matrix.h"


namespace Ariadne {
  namespace LinearAlgebra {
    template < typename R > class TransformationSystem;

    template<typename R> std::ostream& operator<<(std::ostream&, const TransformationSystem<R>&);

    /*! \brief An affine transformation of the unit ball in the supremum norm, representing a zonotopic vector.
     *
     */
    template <typename R>
    class TransformationSystem {
     public:
      /*! \brief Default constructor constructs a zonotopic Vector representation of the origin of dimension \a n. */
      explicit TransformationSystem(size_type n = 0)
        : _centre(n),  _generators(n,0) { }
      
      /*! \brief Construct from centre and directions. */
      explicit TransformationSystem(const Vector<R>& c, const Matrix<R>& m)
        : _centre(c), _generators(m)
      {
        if (c.size()!=m.size1()) {
          throw std::domain_error(
              "The the Matrix of principal directions does not have the same number of rows as the point dimension.");
        }

        this->minimize_generators();
      }
       
      /*! \brief Construct from a Vector. */
      TransformationSystem(const Vector<R>& v)
        : _centre(v), _generators(v.size(),0)
      {
      }

      /*! \brief Construct from an interval Vector. */
      TransformationSystem(const IntervalVector<R>& iv)
        : _centre(iv.size()), _generators(iv.size(),iv.size())
      {
        for(size_type i=0; i!=this->size(); ++i) {
          this->_centre(i) = (iv(i).lower()+iv(i).upper())/2;
          this->_generators(i,i) = (iv(i).upper()-iv(i).lower())/2;
        }

        this->minimize_generators();
      }

      /*! \brief Construct from directions. */
      explicit TransformationSystem(const Matrix<R>& m)
        : _centre(m.size1()), _generators(m)
      {
        this->minimize_generators();
      }
       
      /*! \brief Copy constructor. */
      TransformationSystem(const TransformationSystem<R>& v)
        : _centre(v._centre), _generators(v._generators)
      {
      }

      /*! \brief Copy assignment operator. */
      TransformationSystem<R>& operator=(const TransformationSystem<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      
      /*! \brief A rectangle containing the given zonotope. */
      inline IntervalVector<R> bounding_box() const {
        IntervalVector<R> unit(this->size());
        for(size_type i=0; i!=this->size(); ++i) {
          unit[i]=Interval<R>(-1,1);
        }
        return this->_centre+this->_generators*unit;
      }
      
      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      inline size_type size() const {
        return this->_centre.size();
      }
      
      /*! \brief The number of generators. */
      inline size_type number_of_generators() const {
        return this->_generators.size2();
      }
      
      /*! \brief The centre of the zonotope. */
      inline Vector<R> centre() const {
        return this->_centre;
      }
      
      /*! \brief The Matrix of principle directions. */
      inline Matrix<R> generators() const {
        return this->_generators;
      }
     private:
      // Minimize the generator Matrix
      void minimize_generators(void);
     private:
      /* Zonotopic Vector's centre. */
      Vector<R> _centre;
      /* Zonotopic Vector's principal directions. */
      Matrix<R> _generators;
    
    };
  

  
    /*! \brief The scalar multiple of a zonotopic Vector. */
    template<typename R> 
    inline
    TransformationSystem<R> operator*(const R& s, const TransformationSystem<R>& v)
    {
      return TransformationSystem<R>(s*v.centre(),s*v.generators());
    }
    
    /*! \brief The scalar multiple of a zonotopic Vector. */
    template<typename R> 
    inline
    TransformationSystem<R> operator*(const TransformationSystem<R>& v, const R& s)
    {
      return s*v;
    }
    
    /*! \brief The scalar multiple of a zonotopic Vector. */
    template<typename R> 
    inline
    TransformationSystem<R> operator/(const TransformationSystem<R>& v, const R& s)
    {
      return (1/s)*v;
    }
    
    /*! \brief The sum of two zonotopic Vectors. */
    template<typename R> 
    inline
    TransformationSystem<R> operator+(const TransformationSystem<R>& u, 
                                  const TransformationSystem<R>& v)
    {
      if (u.size()!=v.size()) {
        throw std::domain_error(
          "operator+(TransformationSystem<R>,const TransformationSystem<R>& v): the two zonotopes have different dimension.");
      }
      size_type m=u.generators().size2();
      size_type n=v.generators().size2();
      
      Matrix<R> gen(u.size(),m+n);
      
      for (size_type i=0; i!=u.size(); ++i) {
        for (size_type j=0; j!=m; ++j) {
          gen(i,j)=u.generators()(i,j);
        }
        for (size_type j=0; j!=n; ++j) {
          gen(i,j+m)=v.generators()(i,j);
        }
      }
      return TransformationSystem<R>(u.centre()+v.centre(),gen);
    }
   
    /*! \brief The sum of a zonotopic Vector and a Vector. */
    template<typename R> 
    inline
    TransformationSystem<R> operator+(const TransformationSystem<R>& u, 
                                  const Vector<R>& v)
    {
      return TransformationSystem<R>(u.centre()+v,u.generators());
    }
    
    /*! \brief The sum of a zonotopic Vector and an interval Vector. */
    template<typename R> 
    inline
    TransformationSystem<R> operator+(const TransformationSystem<R>& u, 
                                  const IntervalVector<R>& v)
    {
      return u+TransformationSystem<R>(v);
    }
    
    /*! \brief The sum of an interval Vector and a zonotopic Vector. */
    template<typename R> 
    inline
    TransformationSystem<R> operator+(const IntervalVector<R>& u, 
                                  const TransformationSystem<R>& v)
    {
      return TransformationSystem<R>(u)+v;
    }
    
    /*! \brief The sum of a Vector and a TransformationSystem. */
    template<typename R> 
    inline
    TransformationSystem<R> operator+(const Vector<R>& u, 
                                  const TransformationSystem<R>& v)
    {
      return TransformationSystem<R>(u+v.centre(),v.generators());
    }
    
    /*! \brief The product of a Matrix and a TransformationSystem. */
    template<typename R> 
    inline
    TransformationSystem<R> operator*(const Matrix<R>& A, 
                                  const TransformationSystem<R>& v)
    {
      return TransformationSystem<R>(A*v.centre(),A*v.generators());
    }
    

    /*! The convex hull of \f$[-1,1]*v\f$. */
    template<typename R> 
    LinearAlgebra::TransformationSystem<R>
    symmetrise(const LinearAlgebra::IntervalVector<R>& iv);
      
    
  }
}

#endif /* _ARIADNE_ZONOTOPIC_VECTOR_H */
