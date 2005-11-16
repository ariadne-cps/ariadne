/***************************************************************************
 *            state.h
 *
 *  Sun Jan 23 18:00:21 2005
 *  Copyright  2005  Alberto Casagrande
 *  Email casagrande@dimi.uniud.it
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

/*! \file state.h
 *  \brief A state in Euclidean space.
 */

#ifndef _ARIADNE_STATE_H
#define _ARIADNE_STATE_H

#include <vector>
#include <iostream>
#include <stdexcept>

#include "utility.h"
#include "linear_algebra.h"

namespace Ariadne {
  namespace Geometry {

    template <typename R = Rational> class State;

    template <typename R>
    class State {
     public:
      typedef R Real;
      typedef Real value_type;
      typedef size_t size_type;
     private:
      /*! \brief The vector defining the state */
      boost::numeric::ublas::vector<Real> _vector;

     public:
      State() : _vector(0) { }

      State(size_type dim) : _vector(dim) {
        for(size_type i=0; i!=dimension(); ++i) {
          _vector[i]=Real(0);
        }
      }

      State(size_type dim, const Real default_value) : _vector(dim) {
        for(size_type i=0; i!=dimension(); ++i) {
          _vector[i]=default_value;
        }
      }

      template<class ForwardIterator>
      State(ForwardIterator b, ForwardIterator e) : _vector(std::distance(b,e))
      {
        for(size_type i=0; i!=dimension(); ++i) {
          _vector[i]=*b;
          ++b;
        }

      }

      State(const State& original) : _vector(original._vector) { }

      inline Real& operator[] (size_type index) {
        if (((this->_vector).size() <= index)||(index<0)) {
          throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }

      inline const Real& operator[](size_t index) const {
        if (((this->_vector).size() <= index)||(index<0)) {
            throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }


      /*! \brief Checks equivalence between two states. */
      inline bool operator==(const State<Real> &A) const {
        /* Return false if states have different dimensions */
        if (this->dimension()!=A.dimension()) { return false; }

        /* for each dimension i */
        for (size_t i=0; i<this->dimension(); i++) {
          if (this->_vector[i]!=A._vector[i]) { return false; }
        }

        return true;
      }

      /*! \brief Checks equivalence between two states. */
      inline bool operator!=(const State<Real> &A) const {
        return !( *this == A );
      }

      /*! \brief The dimension of the Euclidean space the state lies in. */
      inline size_t dimension() const {
        return (this->_vector).size();
      }

      inline Real get(size_t index) const {
        if (((this->_vector).size() <= index)||(index<0)) {
            throw std::out_of_range("Out of the vector's range.");
        }
        return  (this->_vector[index]);
      }

      inline void set(size_t index, const Real& r) {
        if (((this->_vector).size() <= index) || (index<0)) {
            throw std::out_of_range("Out of the vector's range.");
        }
        this->_vector[index]=r;
      }

      inline State<Real> &operator=(const State<Real> &A) {
        this->_vector=A._vector;
        return *this;
      }

      template <typename RType>
      friend std::ostream& operator<<(std::ostream &os, const State<RType> &state);

      template <typename RType>
      friend std::istream& operator>>(std::istream &is, State<RType> &state);
    };


    template <typename R>
    std::ostream& operator<< (std::ostream &os, const State<R> &state)
    {
      os << "[";
      if(state.dimension() > 0) {
        os << state[0] ;
        for (size_t i=1; i<state.dimension(); i++) {
          os << ", " << state[i];
        }
      }
      os << "]" ;

      return os;
    }

    template <typename R>
    std::istream& operator>> (std::istream &is, State<R> &state)
    {
      static size_t last_size;

      std::vector<R> v;
      v.reserve(last_size);
      is >> v;
      last_size = v.size();

      state._vector.resize(v.size());
      for(size_t i=0; i!=v.size(); ++i) {
        state._vector[i]=v[i];
      }
      return is;
    }

  }
}

#endif /* _ARIADNE_STATE_H */
