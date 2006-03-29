/***************************************************************************
 *            zonotopic_vector.h
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
 
/*! \file zonotopic_vector.h
 *  \brief Zonotopic vectors (affine images of cuboids).
 */

#ifndef _ARIADNE_ZONOTOPIC_VECTOR_H
#define _ARIADNE_ZONOTOPIC_VECTOR_H

#include <ostream>
#include <string>

#include <vector>

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../base/utility.h"
#include "../base/interval.h"

namespace Ariadne {
  namespace LinearAlgebra {
    template < typename R > class zonotopic_vector;

    template<typename R> std::ostream& operator<<(std::ostream&, const zonotopic_vector<R>&);
    template<typename R> std::istream& operator>>(std::istream&, zonotopic_vector<R>&);

    /*! \brief A zonotopic vector set. 
     * 
     * A zonotope is a set of the form \f$c+Ae\f$, where \f$||e||_{infty}\leq1\f$.
     * The intersection and membership tests are performed using algorithms from: <br/>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  <it>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</it> (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     */
    template <typename R>
    class zonotopic_vector {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R Real;
     public:
      /*! \brief Default constructor constructs a zonotopic vector representation of the origin of dimension \a n. */
      explicit zonotopic_vector(size_type n = 0)
        : _centre(n),  _generators(n,0) { }
      
      /*! \brief Construct from centre and directions. */
      explicit zonotopic_vector(const vector<R>& c, const matrix<R>& m)
        : _centre(c), _generators(m)
      {
        if (c.size()!=m.size1()) {
          throw std::domain_error(
              "The the matrix of principal directions does not have the same number of rows as the point dimension.");
        }

        this->minimize_generators();
      }
       
      /*! \brief Construct from a vector. */
      zonotopic_vector(const vector<R>& v)
        : _centre(v), _generators(v.size(),0)
      {
      }

      /*! \brief Construct from an interval vector. */
      zonotopic_vector(const vector< Interval<R> >& iv)
        : _centre(iv.size()), _generators(iv.size(),iv.size())
      {
        for(size_type i=0; i!=this->size(); ++i) {
          this->_centre(i) = (iv(i).lower()+iv(i).upper())/2;
          this->_generators(i,i) = (iv(i).upper()-iv(i).lower())/2;
        }

        this->minimize_generators();
      }

      /*! \brief Construct from a string literal. */
      explicit zonotopic_vector(const std::string& s)
        : _centre(), _generators()
      {
        std::stringstream ss(s);
        ss >> *this;
      }
      
      /*! \brief Copy constructor. */
      zonotopic_vector(const zonotopic_vector<R>& v)
        : _centre(v._centre), _generators(v._generators)
      {
      }

      /*! \brief Copy assignment operator. */
      zonotopic_vector<R>& operator=(const zonotopic_vector<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      
      /*! \brief A rectangle containing the given zonotope. */
      inline vector< Interval<R> > bounding_box() const {
        vector< Interval<R> > unit(this->size());
        for(size_type i=0; i!=this->size(); ++i) {
          unit[i]=Interval<R>(-1,1);
        }
        return this->_centre+this->_generators*unit;
      }
      
      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      inline size_type size() const {
        return this->_centre.size();
      }
      
      /*! \brief The centre of the zonotope. */
      inline vector<R> centre() const {
        return this->_centre;
      }
      
      /*! \brief The matrix of principle directions. */
      inline matrix<R> generators() const {
        return this->_generators;
      }
      
      friend std::istream& operator>> <> (std::istream& is, zonotopic_vector<R>& r);
     private:
      // Minimize the generator matrix
      inline void minimize_generators(void);
     private:
      /* Zonotopic vector's centre. */
      vector<R> _centre;
      /* Zonotopic vector's principal directions. */
      matrix<R> _generators;
    
    };
  
    template<typename R>
    inline 
    void zonotopic_vector<R>::minimize_generators() 
    {
      /* Remove all zero columns from generators. */
      /* TODO: Remove all linear combinations. */
      std::vector<size_type> nzcols;
      for(size_type j=0; j!=_generators.size2(); ++j) {
        for(size_type i=0; i!=_generators.size1(); ++i) {
          if(_generators(i,j)!=0) {
            nzcols.push_back(i);
            break;
          }
        }
      }
      
      matrix<R> new_gens(this->_generators.size1(),nzcols.size());
      for(size_type j=0; j!=new_gens.size2(); ++j) {
        size_type k=nzcols[j];
        for(size_type i=0; i!=new_gens.size1(); ++i) {
          new_gens(i,j)=this->_generators(i,k);
        }
      }
      
      this->_generators=new_gens;
    }
  
    /*! \brief The scalar multiple of a zonotopic vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator*(const R& s, const zonotopic_vector<R>& v)
    {
      return zonotopic_vector<R>(s*v.centre(),s*v.generators());
    }
    
    /*! \brief The scalar multiple of a zonotopic vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator*(const zonotopic_vector<R>& v, const R& s)
    {
      return s*v;
    }
    
    /*! \brief The scalar multiple of a zonotopic vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator/(const zonotopic_vector<R>& v, const R& s)
    {
      return (1/s)*v;
    }
    
    /*! \brief The sum of two zonotopic vectors. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator+(const zonotopic_vector<R>& u, 
                                  const zonotopic_vector<R>& v)
    {
      if (u.size()!=v.size()) {
        throw std::domain_error(
          "operator+(zonotopic_vector<R>,const zonotopic_vector<R>& v): the two zonotopes have different dimension.");
      }
      size_type m=u.generators().size2();
      size_type n=v.generators().size2();
      
      matrix<R> gen(u.size(),m+n);
      
      for (size_type i=0; i!=u.size(); ++i) {
        for (size_type j=0; j!=m; ++j) {
          gen(i,j)=u.generators()(i,j);
        }
        for (size_type j=0; j!=n; ++j) {
          gen(i,j+m)=v.generators()(i,j);
        }
      }
      return zonotopic_vector<R>(u.centre()+v.centre(),gen);
    }
   
    /*! \brief The sum of a zonotopic vector and a vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator+(const zonotopic_vector<R>& u, 
                                  const vector<R>& v)
    {
      return zonotopic_vector<R>(u.centre()+v,u.generators());
    }
    
    /*! \brief The sum of a zonotopic vector and an interval vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator+(const zonotopic_vector<R>& u, 
                                  const vector< Interval<R> >& v)
    {
      return u+zonotopic_vector<R>(v);
    }
    
    /*! \brief The sum of an interval vector and a zonotopic vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator+(const vector< Interval<R> >& u, 
                                  const zonotopic_vector<R>& v)
    {
      return zonotopic_vector<R>(u)+v;
    }
    
    /*! \brief The sum of a vector and a zonotopic_vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator+(const vector<R>& u, 
                                  const zonotopic_vector<R>& v)
    {
      return zonotopic_vector<R>(u+v.centre(),v.generators());
    }
    
    /*! \brief The product of a matrix and a zonotopic_vector. */
    template<typename R> 
    inline
    zonotopic_vector<R> operator*(const matrix<R>& A, 
                                  const zonotopic_vector<R>& v)
    {
      return zonotopic_vector<R>(A*v.centre(),A*v.generators());
    }
    
    template<typename R> 
    std::ostream&
    operator<<(std::ostream& os, const zonotopic_vector<R>& zv) 
    {
      return os << zv.centre() << "+" << zv.generators() << "[-1,1]^" << zv.size();
    }
    
  }
}

#endif /* _ARIADNE_ZONOTOPIC_VECTOR_H */
