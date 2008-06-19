/***************************************************************************
 *            hybrid_timed_set.h
 *
 *  Copyright  2007  Alberto Casagrande,  Pieter Collins
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
 
#ifndef ARIADNE_HYBRID_TIMED_SET_H
#define ARIADNE_HYBRID_TIMED_SET_H

namespace Ariadne {  
  

    /*! \brief A (basic) set with an associated time value. */
    template<class BS> 
    class HybridTimedSet
      : public tuple<Rational, Integer, DiscreteState, BS>
    {
      typedef Rational Q; typedef Integer Z; typedef DiscreteState DS;
      typedef tuple<Q,Z,DS,BS> _base;
     public:
      typedef typename BS::real_type real_type;
      /*! \brief Constructor. */
      template<class BST> HybridTimedSet(const DiscreteState& ds, const BST& bs) : _base(Q(0),Z(0),ds,BS(bs)) { }
      /*! \brief Constructor. */
      template<class BST> HybridTimedSet(const HybridBasicSet<BST>& hbs) : _base(Q(0),Z(0),hbs.state(),BS(hbs.set())) { }
      /*! \brief Constructor. */
      template<class BST> HybridTimedSet(const Q& t, const Z& n, const DS& ds, const BST& bs)
        : _base(t,n,ds,BS(bs)) { }
      /*! \brief The real time associated to the set. */
      const Q& time() const { return this->first; }
      /*! \brief The discrete time associated to the set. */
      const Z& steps() const { return this->second; }
      /*! \brief The discrete state. */
      const DS& state() const { return this->third; } 
      /*! \brief The untimed basic set. */
      const BS& set() const { return this->fourth; } 
      /*! \brief The dimension set. */
      dimension_type dimension() const { return this->set().dimension(); } 
      
      /*! \brief Comparison operator by time value. */
      bool operator<(const HybridTimedSet<BS>& other) const { 
        return this->time() < other.time()
          || ( this->time()==other.time() && ( this->steps() < other->steps()
            || ( this->steps()==other.steps() && (this->state()<other.state()) ) ) ); }
    };
  
    template<class BS> 
    std::ostream& operator<<(std::ostream& os, const HybridTimedSet<BS>& ts) {
      return os << "{ (" << ts.time() << ","<<ts.set()<<") : (" << ts.state() << ","<< ts.set() << ") }";
    }
   
  
} // namespace Ariadne

#endif /* ARIADNE_TIMED_SET_H */
