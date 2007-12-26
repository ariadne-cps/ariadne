/***************************************************************************
 *            orbit.h
 *
 *  Copyright 2007  Pieter Collins
 * 
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
 
/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */

#ifndef ARIADNE_ORBIT_H
#define ARIADNE_ORBIT_H

#include <iosfwd>
#include <vector>

#include "base/exceptions.h"
#include "geometry/declarations.h"

namespace Ariadne {

  namespace Output { class epsstream; }


  namespace Geometry {

    template<class T, class S> class TimedSet;

    template<class T, class ES, class RS> class Orbit;


    template<class T, class S> 
    class TimedSet 
    {
     public:
      template<class TT, class SS> TimedSet(const TT& t, const SS& s) : _t(t), _s(s) { }
      const T& time() const { return this->_t; }
      const S& set() const { return this->_s; }
     private:
      T _t;
      S _s;
    };

    template<class T, class ES, class RS=ES> 
    class dSet 
      : public TimedSet<T,ES>
    {
     public:
      template<class TT, class EST, class RST> dSet(const TT& t, const RST& rs, const EST& es) : TimedSet<T,ES>(t,es), _rs(rs) { }
      const ES& evolved_set() const { return this->TimedSet<T,ES>::set(); }
      const RS& reached_set() const { return this->_rs; }
     private:
      RS _rs;
    };


    /*! \ingroup Evolution
     *  \brief An orbit of a dynamic system.
     */
    template<class T>
    class OrbitInterface
    {
     public:
      virtual ~OrbitInterface();
      virtual Output::epsstream& write(Output::epsstream& eps) const = 0;
      virtual std::ostream& write(std::ostream& eps) const = 0;
    };

    template<class T>
    OrbitInterface<T>::~OrbitInterface() {
    }

    template<class T>
    Output::epsstream& operator<<(Output::epsstream& eps, const OrbitInterface<T>& orbit) {
      return orbit.write(eps);
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const OrbitInterface<T>& orbit) {
      return orbit.write(os);
    }


    /*! \ingroup Evolution
     *  \brief An orbit of a discrete-time dynamic system.
     * 
     */
    template<class BS>
    class Orbit<Numeric::Integer,BS>
      : public OrbitInterface<Numeric::Integer>
    {
      typedef Numeric::Integer T;
     private:
      std::vector< TimedSet<T,BS> > _data;
     public:
      //@{
      //! \name Constructors
      /*! \brief Construct an orbit with an initial set. */
      explicit Orbit() :  _data() { }

      /*! \brief Construct an orbit with an initial set. */
      explicit Orbit(const BS& initial_set) :  _data(1u,TimedSet<T,BS>(0,initial_set)) { }
      
      /*! \brief Copy constructor. */
      Orbit(const Orbit<T,BS>& orbit) : _data(orbit._data) { }
      
      /*! \brief Copy assignment. */
      Orbit<T,BS>& operator=(const Orbit<T,BS>& orbit) { this->_data=orbit._data; return *this; }
      
      /*! \brief Conversion constructor. */
      template<class S> explicit Orbit(const Orbit<T,S>& orbit) : _data() { 
        for(uint i=0; i<=orbit.size(); ++i) { this->_data.push_back(TimedSet<BS,T>(orbit[i])); } }

      //@}
     
      //@{
      //! \name Data access
      /*! \brief Append a time and set to the orbit. */
      template<class S> void push_back(const T& t, const S& s) { this->_data.push_back(TimedSet<T,BS>(t,s)); } 

      /*! \brief Clear the orbit. */
      void clear() { this->_data.clear(); }

      /*! The number of steps used to compute the orbit, excluding the inital set. */
      size_type steps() const { return this->_data.size()-1; }

      /*! Return the initial time. */
      const T& initial_time() const { return this->_data.front().time(); }

      /*! Return the initial set. */
      const BS& initial_set() const { return this->_data.front().set(); }

      /*! Return the final time. */
      const T& final_time() const { return this->_data.back().time(); }

      /*! Return the final set. */
      const BS& final_set() const { return this->_data.back().set(); }

      /*! Return the \a i<sup>th</sup> timed set. */
      const TimedSet<T,BS>& operator[](size_type i) const { return this->_data[i]; }

      /*! Return the reach set between times \a i-1 and \a i. */
      const TimedSet<T,BS>& reach(size_type i) const { return this->_data[i]; }

      /*! Return the union of the timed sets. */
      const ListSet<BS> reach() const { 
        ListSet<BS> result; for(size_type i=0; i!=this->_data.size(); ++i) { result.adjoin(this->_data[i].set()); } return result; }

      //@}

      //@{
      //! \name Input-output operators
      std::ostream& write(std::ostream& os) const;
      //! \name Input-output operators
      virtual Output::epsstream& write(Output::epsstream& eps) const;
      //@}
    };

    
    template<class T, class BS>
    std::ostream& 
    operator<<(std::ostream& os, const TimedSet<T,BS> ts) 
    {
      return os<<ts.time()<<":"<<ts.set();
    }
    
    template<class BS>
    std::ostream& 
    Orbit<Numeric::Integer,BS>::write(std::ostream& os) const
    {
      return os<<"Orbit"<<this->_data;
    }
    

    /*! \ingroup Evolution
     *  \brief An orbit of a discrete-time dynamic system.
     * 
     */
    template<class ES, class RS>
    class Orbit<Numeric::Rational, ES, RS>
      : public OrbitInterface<Numeric::Rational>
    {
      typedef Numeric::Rational T;
     private:
      std::vector< dSet<T,ES,RS> > _data;
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor. */
      explicit Orbit()
        :  _data() { }

      /*! \brief Construct an orbit with an initial set. */
      explicit Orbit(const ES& initial_set)
        :  _data(1u,dSet<T,ES,RS>(T(0),initial_set,initial_set)) { }
      
      /*! \brief Copy constructor. */
      Orbit(const Orbit<T,ES,RS>& orbit)
        : _data(orbit._data) { }
      
      /*! \brief Copy assignment. */
      Orbit<T,ES,RS>& operator=(const Orbit<T,ES,RS>& orbit) { 
        if(this!=&orbit) { this->_data=orbit._data; } return *this; }
      
      //@}
     
      //@{
      //! \name Data access
      /*! \brief Append the new time, the intermediate reach set and the final set to the orbit. */
      void push_back(const T& t, const RS& rs, const ES& es) { 
        this->_data.push_back(dSet<T,ES,RS>(t,rs,es)); }

      /*! \brief Clear the orbit. */
      void clear() { this->_data.clear(); }

      /*! The number of steps used to compute the orbit, excluding the inital set. */
      size_type steps() const { return this->_data.size()-1; }

      /*! Return the initial timed set. */
      const TimedSet<T,ES>& initial() const { return this->_data.front(); }

      /*! Return the initial time. */
      const T& initial_time() const { return this->_data.front().time(); }

      /*! Return the initial set. */
      const ES& initial_set() const { return this->_data.front().evolved_set(); }

      /*! Return the final timed set. */
      const TimedSet<T,ES>& final() const { return this->_data.back(); }

      /*! Return the final time. */
      const T& final_time() const { return this->_data.back().time(); }

      /*! Return the final set. */
      const ES& final_set() const { return this->_data.back().evolved_set(); }

      /*! Return the \a i<sup>th</sup> timed set. */
      const TimedSet<T,ES>& operator[](size_type i) const { return this->_data[i]; }

      /*! Return the union of the reached sets. */
      const ListSet<RS> reach() const { 
        ListSet<RS> result; for(size_type i=1; i<this->_data.size(); ++i) { result.adjoin(this->_data[i].reached_set()); } return result; }
      //@}

      //@{
      //! \name Input-output operators
      std::ostream& write(std::ostream& os) const;
      //! \name Input-output operators
      virtual Output::epsstream& write(Output::epsstream& eps) const;
      //@}

    };

  }
}

#include "orbit.template.h"

#endif /* ARIADNE_ORBIT_H */
