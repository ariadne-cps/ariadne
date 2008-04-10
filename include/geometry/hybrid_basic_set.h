/***************************************************************************
 *            hybrid_basic_set.h
 *
 *  Copyright  2006-7  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#ifndef ARIADNE_HYBRID_BASIC_SET_H
#define ARIADNE_HYBRID_BASIC_SET_H

#include <string>
#include <iostream>

#include "geometry/declarations.h"
#include "geometry/discrete_state.h"
#include "geometry/hybrid_space.h"

namespace Ariadne {
  namespace Geometry {

    class basic_set_tag;

    class HybridSystemError : public std::runtime_error {
     public:
      HybridSystemError(const std::string& what) : std::runtime_error(what) { }
    };

    template<class R> class HybridPoint;
    


  /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a single basic set for in a discrete mode.
     */
    template<class BS> 
    class HybridBasicSet 
    {
      typedef typename BS::real_type R;
     public:
      /*! \brief */
      typedef typename BS::real_type real_type;
      /*! \brief */
      typedef HybridSpace space_type;
      /*! \brief */
      typedef HybridPoint<R> state_type;
      /*! \brief */
      typedef DiscreteState discrete_state_type;
      /*! \brief */
      typedef BS continuous_set_type;
      /*! \brief */
      typedef basic_set_tag set_category;

      // No default constructor as original set need not have default constructor
      // HybridBasicSet() : BS(), _discrete_state() { }
      /*! \brief */
      template<class M, class S> HybridBasicSet(const M& q, const S& s) : _state(q), _set(s) { }
      /*! \brief */
      const DiscreteState& state() const { return this->_state; }
      /*! \brief */
      const BS& set() const { return this->_set; }
    
      /*! \brief */
      dimension_type dimension() const { return this->_set.dimension(); }
    
      /*! \brief */
      bool operator==(const HybridBasicSet<BS>& other) const { 
        return this->_state()==other.state() && this->set()==other.set(); }
      /*! \brief */
      bool operator!=(const HybridBasicSet<BS>& other) const { return !(*this==other); }
      /*! \brief */
      bool operator<(const HybridBasicSet<BS>& other) const { return this->_state < other->_state; }
     private:
      discrete_state_type _state;
      continuous_set_type  _set;
    };
  
    template<class BS> inline 
    std::ostream& operator<<(std::ostream& os, const HybridBasicSet<BS>& hs) {
      return os << "{ q=" << hs.state() << ", s=" << hs.set() << " }";
    }

    template<class R>
    class HybridBox
      : public HybridBasicSet< Box<R> >
    {
     public:
      template<class HS> HybridBox(const HS& hs) : HybridBasicSet< Box<R> >(hs) { }
      template<class M, class S> HybridBox(const M& q, const S& s) : HybridBasicSet< Box<R> >(q,s) { }
    };

    template<class R>
    class HybridZonotope
      : public HybridBasicSet< Zonotope<R> >
    {
     public:
      template<class HS> HybridZonotope(const HS& hs) : HybridBasicSet< Box<R> >(hs) { }
      template<class M, class S> HybridZonotope(const M& q, const S& s) : HybridBasicSet< Zonotope<R> >(q,s) { }
    };

    template<class R>
    class HybridPolytope
      : public HybridBasicSet< Polytope<R> >
    {
     public:
      template<class HS> HybridPolytope(const HS& hs) : HybridBasicSet< Box<R> >(hs) { }
      template<class M, class S> HybridPolytope(const M& q, const S& s) : HybridBasicSet< Polytope<R> >(q,s) { }
    };


    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a single basic set for in a discrete mode.
     */
    template<class BS> 
    class HybridTimedBasicSet 
    {
      typedef typename BS::real_type R;
     public:
      typedef typename BS::real_type real_type;
      typedef DiscreteState discrete_state_type;
      typedef basic_set_tag set_category;

      /*! \brief */
      template<class A> HybridTimedBasicSet(const time_type& t, const discrete_time_type& n, const discrete_state_type& q, const A& a)
        : _time(t), _steps(n), _discrete_state(q), _continuous_state_set(a) { }
      /*! \brief */
      const time_type& time() const { return this->_time; }
      /*! \brief */
      const discrete_time_type& steps() const { return this->_steps; }
      /*! \brief */
      const discrete_state_type& discrete_state() const { return this->_discrete_state; }
      /*! \brief */
      const BS& continuous_state_set() const { return this->_continuous_state_set; } 

      /*! \brief */
      Box<R> bounding_box() const { return this->_continuous_state_set.bounding_box(); }

      /*! \brief */
      bool operator==(const HybridTimedBasicSet<BS>& other) const { 
        return this->_time==other._time
          && this->_steps==other._steps
          && this->_discrete_state==other._discrete_state 
          && this->_continuous_state_set==other._continuous_state_set; }
      /*! \brief */
      bool operator!=(const HybridTimedBasicSet<BS>& other) const { return !(*this==other); }
      /*! \brief */
      bool operator<(const HybridTimedBasicSet<BS>& other) const { return this->_time<other._time; }
     private:
      time_type _time;
      discrete_time_type _steps;
      discrete_state_type _discrete_state;
      BS _continuous_state_set;
    };
  
    template<class BS> inline 
    std::ostream& operator<<(std::ostream& os, const HybridTimedBasicSet<BS>& hs) {
      return os << "{ t=" << hs.time() << ", n=" << hs.steps()
                << ", q=" << hs.discrete_state() << ", s=" << hs.continuous_state_set() 
                << "}";
    }

   
  }
}

#include "hybrid_basic_set.inline.h"

#endif /* ARIADNE_HYBRID_BASIC_SET_H */
