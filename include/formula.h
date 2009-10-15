/***************************************************************************
 *            formula.h
 *
 *  Copyright 2008-9  Pieter Collins
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


/*! \file formala.h
 *  \brief Formulae over variables
 */

#ifndef ARIADNE_FORMULA_H
#define ARIADNE_FORMULA_H

#include <cstdarg>
#include <iostream>
#include <string>


#include "macros.h"
#include "pointer.h"
#include "container.h"
#include "stlio.h"

#include "operators.h"
#include "variables.h"
#include "expression.h"
#include "assignment.h"
#include "space.h"

#include "numeric.h"

namespace Ariadne {





class Enumeration;
class Integer;
class Real;

//! \brief An ASCII string.
typedef std::string String;

class StateSpace;

class EnumeratedType;
class EnumeratedValue;



template<class T> class Variable;
template<class R> class Expression;

typedef Constant<Real> RealConstant;

typedef Variable<EnumeratedValue> EnumeratedVariable;
typedef Variable<String> StringVariable;
typedef Variable<Integer> IntegerVariable;
typedef Variable<Real> RealVariable;

typedef Space<Real> RealSpace;

typedef Expression<EnumeratedValue> EnumeratedExpression;
typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;

typedef Expression<bool> DiscretePredicate;
typedef Expression<tribool> ContinuousPredicate;

typedef String Identifier;
typedef Set<Identifier> VariableSet;
typedef Set< Variable<Real> > RealVariableSet;



//! \brief A discrete event.
struct Event {
  public:
    Event() {  } // Constructs an invalid event.
    Event(const std::string& name) { _names.push_back(name); this->_id=_names.size()-1; }
    Event(int n) { std::stringstream ss; ss<<"e"<<n; _names.push_back(ss.str()); this->_id=_names.size()-1; ; }
    int id() const { return this->_id; }
    const std::string& name() const { return _names[this->_id]; }
    bool operator==(const Event& other) const { return this->id()==other.id(); }
    bool operator!=(const Event& other) const { return this->id()!=other.id(); }
    bool operator<(const Event& other) const { return this->id()<other.id(); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Event& self) { return os << self.name(); }
  private:
    uint _id;
  private:
    static std::vector<std::string> _names;
};


//! \brief A finite or cofinite set of discrete events.
struct EventSet {
  public:
    EventSet() : _events(), _is_complement(false) { }
    EventSet(const std::set<Event>& s) : _events(s), _is_complement(false) { }
    EventSet(const std::vector<Event>& v) : _events(v.begin(),v.end()), _is_complement(false) { }
    EventSet(const Event& e) : _events(), _is_complement(false) { this-> _events.insert(e); }

    static EventSet none() { EventSet r; return r; }
    static EventSet all() { EventSet r; r._is_complement=true; return r; }

    bool operator==(const EventSet& other) const {
        return this->_events==other._events && this->_is_complement==other._is_complement; }
    bool empty() const { return !_is_complement && _events.empty(); }
    bool finite() const { return !_is_complement; }
    uint size() const { ARIADNE_ASSERT_MSG(this->finite(),"EventSet "<<*this<<" is infinite"); return _events.size(); }
    const Event& front() const { ARIADNE_ASSERT_MSG(this->finite(),"EventSet "<<*this<<" is infinite"); return *this->_events.begin(); }

    bool contains(const Event& e) const {
        return Ariadne::contains(this->_events,e) xor this->_is_complement; }
    EventSet& insert(const Event& e) {
        if(_is_complement) { this->_events.erase(e); } else { this->_events.insert(e); } return *this; }
    EventSet& erase(const Event& e) {
        if(_is_complement) { this->_events.insert(e); } else { this->_events.erase(e); } return *this; }
    EventSet& adjoin(const EventSet& s) { assert(!_is_complement && !s._is_complement);
        Ariadne::adjoin(this->_events,s._events);
        return *this; }

    EventSet& restrict(const EventSet& s) {
        if(this->_is_complement) {
            if(s._is_complement) { Ariadne::adjoin(this->_events,s._events); }
            else { this->_events=difference(s._events,this->_events); this->_is_complement=false; }
        } else {
            if(s._is_complement) { Ariadne::restrict(this->_events,s._events); }
            else { Ariadne::remove(this->_events,s._events); }
        }
        return *this; }

    EventSet& remove(const EventSet& s) {
        if(this->_is_complement) {
            if(s._is_complement) { this->_events=difference(s._events,this->_events); this->_is_complement=false; }
            else { Ariadne::adjoin(this->_events,s._events); }
        } else {
            if(s._is_complement) { Ariadne::remove(this->_events,s._events); }
            else { Ariadne::restrict(this->_events,s._events); }
        } return *this; }
    friend EventSet operator!(const Event& e);
    friend EventSet operator!(const EventSet& s);
    friend EventSet operator!(const std::set<Event>& s);

    friend std::ostream& operator<<(std::ostream& os, const EventSet& self);
  private:
    std::set<Event> _events;
    bool _is_complement;
};

inline EventSet operator!(const Event& e) { return !(EventSet(e)); }
inline EventSet operator!(const EventSet& s) { EventSet r; r._events=s._events; r._is_complement=!s._is_complement; return r; };
inline EventSet operator!(const std::set<Event>& s) { return !EventSet(s); }

inline EventSet join(const EventSet& s1, const EventSet& s2) { EventSet r(s1); r.adjoin(s2); return r; }
inline EventSet intersection(const EventSet& s1, const EventSet& s2) { EventSet r(s1); r.restrict(s2); return r; }
inline EventSet difference(const EventSet& s1, const EventSet& s2) { EventSet r(s1); r.remove(s2); return r; }

inline std::ostream& operator<<(std::ostream& os, const EventSet& self) {
    if(self._is_complement) { os<<"!"; } return os << self._events; }





/*! \brief An internal enumerated type, suitable for use when a finite type is needed.
 *  \sa EnumeratedType, EnumeratedVariable, DiscretePredicate
 */
struct EnumeratedType {
  public:
    //! \brief Create a new enumerated type with name \a name and values \a values.
    EnumeratedType(const std::string& name, const std::vector<std::string>& values) : _name(name), _values(values) {
        for(uint i=0; i!=this->_values.size(); ++i) {
            assert(this->_values[i]!="");
            for(uint j=i+1; j!=this->_values.size(); ++j) {
                assert(this->_values[i]!=this->_values[j]); } } }
    //! \brief The name of the type.
    const std::string& name() const { return this->_name; }
    //! \brief The number of values the type can take.
    uint size() const { return _values.size(); }
    //! \brief Returns \c true  if \a str is a valid value for the type.
    bool has_value(const std::string& str) const {
        for(uint i=0; i!=_values.size(); ++i) { if(_values[i]==str) { return true; } } return false; }
    //! \brief The \a i<sup>th</sup> value the type can take.
    const std::string& value(uint i) const { return _values.at(i); }
    int index(const std::string& str) const {
        for(uint i=0; i!=_values.size(); ++i) { if(_values[i]==str) { return i; } } assert(false); }
    //! \brief Test if two type variables describe the same type.
    bool operator==(const EnumeratedType& other) const;
    //! \brief Test if two types are different.
    bool operator!=(const EnumeratedType& other) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const EnumeratedType& self) {
        os << self.name() << "=";
        for(uint i=0; i!=self._values.size(); ++i) { os<<(i==0?'{':',')<<self._values[i]; } return os << '}'; }
  private:
    std::string _name;
    std::vector<std::string> _values;
};

//! \brief A predicate over some discrete variables.
//! \details \sa EnumeratedValue, EnumeratedVariable, DiscretePredicate
struct EnumeratedValue {
  public:
    //! \brief Construct a discrete value based on the string \a str.
    EnumeratedValue(const String& str) : _value(str) { }
    //! \brief The string representation of the value.
    operator const String& () const { return this->_value; }
    //! \brief Equality operator.
    bool operator==(const EnumeratedValue& s) const { return this->_value==s._value; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const EnumeratedValue& self) {
        return os << self._value; }
  private:
    std::string _value;
};





template<class R> class ExpressionInterface;
template<class R> class VariableExpression;
template<class R> class ConstantExpression;
template<class R, class Op, class A> class UnaryExpression;
template<class R, class Op, class A1, class A2> class BinaryExpression;
template<class R> class Expression;
template<class R> std::ostream& operator<<(std::ostream&, const Expression<R>& f);
template<class R> std::ostream& operator<<(std::ostream&, const ExpressionInterface<R>& f);









typedef ExtendedVariable<Real> ExtendedRealVariable;
typedef DottedVariable<Real> DottedRealVariable;
typedef PrimedVariable<Real> PrimedRealVariable;
typedef PrimedVariable<String> PrimedStringVariable;
typedef PrimedVariable<Integer> PrimedIntegerVariable;
typedef PrimedVariable<EnumeratedValue> PrimedEnumeratedVariable;

typedef Assignment<EnumeratedVariable,EnumeratedExpression> EnumeratedAssignment;
typedef Assignment<PrimedEnumeratedVariable,EnumeratedExpression> EnumeratedUpdate;
typedef Assignment<IntegerVariable,IntegerExpression> IntegerAssignment;
typedef Assignment<PrimedIntegerVariable,IntegerExpression> IntegerUpdate;
typedef Assignment<StringVariable,StringExpression> StringAssignment;
typedef Assignment<PrimedStringVariable,StringExpression> StringUpdate;
typedef Assignment<RealVariable,RealExpression> RealAssignment;
typedef Assignment<PrimedRealVariable,RealExpression> RealUpdate;
typedef Assignment<DottedRealVariable,RealExpression> RealDynamic;

typedef Assignment<RealVariable,RealExpression> RealAlgebraicAssignment;
typedef Assignment<PrimedRealVariable,RealExpression> RealUpdateAssignment;
typedef Assignment<DottedRealVariable,RealExpression> RealDifferentialAssignment;
typedef Assignment<ExtendedRealVariable,RealExpression> ExtendedRealAssignment;
typedef Assignment<DottedRealVariable,RealExpression> DottedRealAssignment;
typedef Assignment<PrimedRealVariable,RealExpression> PrimedRealAssignment;








} // namespace Ariadne

#endif // ARIADNE_FORMULA_H
