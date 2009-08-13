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

#include "function.h"
#include "operators.h"

#include "macros.h"
#include "pointer.h"
#include "container.h"

#include "expression.h"

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

typedef Variable<EnumeratedValue> EnumeratedVariable;
typedef Variable<String> StringVariable;
typedef Variable<Integer> IntegerVariable;
typedef Variable<Real> RealVariable;

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




//! \brief An assignment statement.
template<class LHS, class RHS>
struct Assignment
{
    Assignment(const LHS& l, const RHS& r) : lhs(l), rhs(r) { }
    LHS lhs; RHS rhs;
};

template<class LHS, class RHS> bool operator<(const Assignment<LHS,RHS>& a1, const Assignment<LHS,RHS>& a2) {
    return a1.lhs < a2.lhs;
}

template<class LHS, class RHS> inline std::ostream& operator<<(std::ostream& os, const Assignment<LHS,RHS>& a) {
    return os<<a.lhs<<"="<<a.rhs;
}



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



class Space;
std::ostream& operator<<(std::ostream& os, const Space& spc);

/*! \brief A space defined as a list of named variables.
 *  \details \sa Variable
 */
struct Space
{
  public:
    //! \brief The trivial space \f$\R^0\f$.
    Space() : _variables() { }
    Space(const std::set<Identifier>& vs) : _variables(vs.begin(),vs.end()) { }
    unsigned int size() const { return _variables.size(); }
    //! \brief The dimension of the space.
    unsigned int dimension() const { return _variables.size(); }
    //! \brief The \a i<sup>th</sup> named variable.
    const Identifier& operator[](unsigned int i) const { return _variables.at(i); }
    const Identifier& variable(unsigned int i) const { return _variables.at(i); }

    //! \brief The index of the named variable \a v.
    unsigned int index(const Identifier& v) const {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(v==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<v<<" is not in the Space "<<*this);
        return _variables.size(); }
    //! \brief Append the named variable \a v to the variables defining the space; ignores if the variable is already in the space.
    Space& insert(const Identifier& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(_variables[i]==v) { return *this; } }
        _variables.push_back(v); return *this; }
    //! \brief Adjoins the variables in \a spc.
    Space& adjoin(const Space& spc) {
        for(uint i=0; i!=spc._variables.size(); ++i) { this->insert(spc._variables[i]); } return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space& append(const Identifier& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(_variables[i]!=v,"Variable "<<v<<" is already a variable of the StateSpace "<<*this);
        }
        _variables.push_back(v); return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space& operator,(const Identifier& v) { return this->append(v); }
  public:
    //! \brief The space containing the union of the variables.
    friend Space join(const Space&, const Space&);
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const Space&);
  private:
    std::vector<Identifier> _variables;
};

inline std::ostream& operator<<(std::ostream& os, const Space& spc) { return os << spc._variables; }

inline Space operator,(const Identifier& v1, const Identifier& v2) {
    Space r; r,v1,v2; return r; }

inline Space join(const Space& spc1, const Space& spc2) {
    Space r(spc1); r.adjoin(spc2); return r; }













template<class R> class ExpressionInterface;
template<class R> class VariableExpression;
template<class R> class ConstantExpression;
template<class R, class O, class A> class UnaryExpression;
template<class R, class O, class A1, class A2> class BinaryExpression;
template<class R> class Expression;
template<class R> std::ostream& operator<<(std::ostream&, const Expression<R>& f);
template<class R> std::ostream& operator<<(std::ostream&, const ExpressionInterface<R>& f);










template<class VAR> struct Next;
template<class T> Next< Variable<T> > next(const Variable<T>& v);



template<class T>
struct Next< Variable<T> >
{
  private:
    typedef Assignment< Next< Variable<T> >, Expression<T> > AssignmentType;
    Variable<T> _base;
  public:
    AssignmentType operator=(const T& val) { return AssignmentType(*this,Expression<T>(val)); }
    AssignmentType operator=(const Variable<T>& var) { return AssignmentType(*this,Expression<T>(var)); }
    AssignmentType operator=(const Expression<T>& expr) { return AssignmentType(*this,expr); }
    friend Next< Variable<T> > next<>(const Variable<T>&);
  private:
    Next(const Variable<T>& v) : _base(v) { }
};

template<>
struct Next<EnumeratedVariable>
{
    typedef EnumeratedValue T;
    EnumeratedVariable base;
    typedef Assignment<Next<EnumeratedVariable>,EnumeratedExpression> EnumeratedAssignment;
//     friend Next<EnumeratedVariable> next(const EnumeratedVariable& v);
    EnumeratedAssignment operator=(const String& str) { return EnumeratedAssignment(*this,EnumeratedExpression(EnumeratedValue(str))); }
    EnumeratedAssignment operator=(const EnumeratedValue& val) { return EnumeratedAssignment(*this,EnumeratedExpression(val)); }
    EnumeratedAssignment operator=(const EnumeratedVariable& var) { return EnumeratedAssignment(*this,EnumeratedExpression(var)); }
    EnumeratedAssignment operator=(const EnumeratedExpression& expr) { return EnumeratedAssignment(*this,expr); }
    friend std::ostream& operator<<(std::ostream& os, const Next<EnumeratedVariable>& self) {
        return os << "next(" << self.base.name() << ")"; }
    friend Next< Variable<T> > next<>(const Variable<T>&);
  private:
    Next(const EnumeratedVariable& v) : base(v) { }
};

template<>
struct Next<RealVariable>
{
    typedef Real T;
    RealVariable base;
    typedef Assignment<Next<RealVariable>,RealExpression> RealAssignment;
    friend Next<RealVariable> next(const RealVariable& v);
    RealAssignment operator=(const double& x) { return RealAssignment(*this,RealExpression(ConstantExpression<Real>(x))); }
    RealAssignment operator=(const RealVariable& var) { return RealAssignment(*this,RealExpression(var)); }
    RealAssignment operator=(const RealExpression& fn) { return RealAssignment(*this,fn); }
    friend std::ostream& operator<<(std::ostream& os, const Next<RealVariable>& self) {
        return os << "next(" << self.base.name() << ")"; }
    friend Next< Variable<T> > next<>(const Variable<T>&);
  private:
    Next(const RealVariable& v) : base(v) { }
};

template<class T> inline Next< Variable<T> > next(const Variable<T>& v) { return Next< Variable<T> >(v); }

template<class T> inline std::ostream& operator<<(std::ostream& os, const Next< Variable<T> >& nv) { return nv.write(os); }



template<class T> class DottedVariable;

template<class T>
struct DottedVariable
{
    Variable<T> base;
    Assignment<DottedVariable<T>,Expression<T> > operator=(const Expression<T>& f) {
        return Assignment<DottedVariable<T>,Expression<T> >(*this,f); }
    Assignment<DottedVariable<T>,Expression<T> > operator=(const Variable<T>& v) {
        return Assignment<DottedVariable<T>,Expression<T> >(*this,Expression<T>(v)); }
    friend std::ostream& operator<<(std::ostream& os, const DottedVariable<T>& self) {
        return os << "dot(" << self.base.name() << ")"; }
    explicit DottedVariable(const Variable<T>& v) : base(v) { }
};

inline DottedVariable<Real> dot(const Variable<Real>& v) { return DottedVariable<Real>(v); }

typedef Assignment<EnumeratedVariable,EnumeratedExpression> EnumeratedAssignment;
typedef Assignment<Next<EnumeratedVariable>,EnumeratedExpression> EnumeratedUpdate;
typedef Assignment<RealVariable,RealExpression> RealAssignment;
typedef Assignment<Next<RealVariable>,RealExpression> RealUpdate;
typedef Assignment<DottedVariable<Real>,Expression<Real> > RealDynamic;







} // namespace Ariadne

#endif // ARIADNE_FORMULA_H
