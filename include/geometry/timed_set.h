/***************************************************************************
 *            timed_set.h
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
 
#ifndef ARIADNE_TIMED_SET_H
#define ARIADNE_TIMED_SET_H

namespace Ariadne {  
  namespace Geometry {

    /*! \brief A (basic) set with an associated time value. */
    template<class T, class S> 
    class TimedSet 
    {
     public:
      /*! \brief The type used to record the time. */
      typedef T time_type;
      /*! \brief The type used for the set. */
      typedef S set_category;
      /*! \brief The type used for real numbers. */
      typedef typename S::real_type real_type;

      /*! \brief Constructor. */
      TimedSet(const S& s)
        : _time(0), _set(s) { }

      /*! \brief Constructor. */
      template<class A>
      TimedSet(const time_type& t, const A& a)
        : _time(t), _set(a) { }

      /*! \brief Constructor. */
      template<class A0, class A1>
      TimedSet(const time_type& t, const A0& a0, const A1& a1)
        : _time(t), _set(a0,a1) { }

      /*! \brief The time associated to the set. */
      const T& time() const { return this->_time; }
      /*! \brief The untimed set. */
      const S& set() const { return this->_set; } 
      /*! \brief The untimed set. */
      dimension_type dimension() const { return this->_set.dimension(); } 
      
      /*! \brief Equality operator. */
      bool operator==(const TimedSet<T,S>& other) const { 
        return this->_time == other._time && this->_set==other._set; }
      /*! \brief Inequality operator. */
      bool operator!=(const TimedSet<T,S>& other) const { return !(*this==other); }
      /*! \brief Comparison operator by time value. */
      bool operator<(const TimedSet<T,S>& other) const { return this->_time < other._time; }
     private:
      T _time;
      S _set;
    };
  
    template<class T, class S> 
    std::ostream& operator<<(std::ostream& os, const TimedSet<T,S>& ts) {
      return os << "{ " << ts.time() << " : " << ts.set() << " }";
    }
   
  }
}

#endif /* ARIADNE_TIMED_SET_H */
