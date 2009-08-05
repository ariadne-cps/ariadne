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

#include "polynomial.h"

namespace Ariadne {





class Enumeration;
class Integer;
class Real;

//! \brief An ASCII string.
typedef std::string String;

class StateSpace;

class EnumeratedType;
class EnumeratedValue;



template<class T> class Quantity;
template<class T> class Variable;
template<class R> class Formula;

typedef Variable<EnumeratedValue> EnumeratedVariable;
typedef Variable<Integer> IntegerVariable;
typedef Variable<Real> RealVariable;

typedef Formula<EnumeratedValue> EnumeratedFormula;
typedef Formula<Integer> IntegerFormula;
typedef Formula<Real> RealFormula;

typedef Formula<bool> DiscretePredicate;
typedef Formula<tribool> ContinuousPredicate;





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



//! \brief A named variable.
struct NamedVariable {
  public:
    //! \brief The internal name of the variable.
    const String& name() const { return *this->_name_ptr; }
    //! \brief Equality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator==(const NamedVariable& v) const { return this->_name_ptr==v._name_ptr; }
    //! \brief Inequality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator!=(const NamedVariable& v) const { return !(*this==v); }
    //! \brief Less-than comparison orders %Variable objects lexicographically by name.
    bool operator<(const NamedVariable& v) const { return *this->_name_ptr<*v._name_ptr; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const NamedVariable& v) {
        return os << v.name(); }
  protected:
    //! \brief Create an internal enumerated variable with the given \a name.
    NamedVariable(const String& name)
        : _name_ptr(new std::string(name)) { }
    NamedVariable(const String& name, const EnumeratedType& type)
        : _name_ptr(new std::string(name)), _type_ptr(new EnumeratedType(type)) { }
    //! \brief The internal type of the variable.
    const EnumeratedType& type() const { return *this->_type_ptr; }
  private:
    shared_ptr<std::string> _name_ptr;
    shared_ptr<EnumeratedType> _type_ptr;

};


//! \brief A named discrete variable.
struct DiscreteVariable : public NamedVariable {
  protected:
    DiscreteVariable(const std::string& name) : NamedVariable(name) { }
    DiscreteVariable(const std::string& name, const EnumeratedType& type) : NamedVariable(name,type) { }
};


//! \brief A named integer variable, suitable for use in formulae and in defining the discrete state.
//! Equality is performed using references by name, so different variables may have the same "name".
//! \tparam T The type of mathematical object the variable represents; e.g. String, Integer, Real
//! \sa Formula.
template<class T> struct Variable : public NamedVariable
{
  public:
    //! \brief Construct a new named variable with name \a name. If another variable has the same \c name,
    //! the two are \e not considered to be the same quantity.
    Variable(String name);
    //! \brief Copy constructer. Returns a C++ variable representing the same named variable.
    Variable(const Variable<T>& v);
    //! \brief Makes an assignment object.
    Assignment< Variable<T>, Formula<T> > operator=(const T& val) const;
};


//! \brief A named variable which can take values of a given enumerated type.
//! \details \sa EnumeratedType, EnumeratedVariable, DiscretePredicate
template<> struct Variable<EnumeratedValue> : public DiscreteVariable {
  public:
    //! \brief Create an internal enumerated variable with the given \a name and \a type.
    Variable(const String& name, const EnumeratedType& type) : DiscreteVariable(name,type) { }
    //! \brief The internal type of the variable.
    const EnumeratedType& type() const { return this->NamedVariable::type(); }
    //! \brief Makes a discrete assignment object.
    Assignment<EnumeratedVariable,EnumeratedFormula> operator=(const String& val) const;
    //! \brief Makes a discrete assignment object.
    Assignment<EnumeratedVariable,EnumeratedFormula> operator=(const EnumeratedValue& val) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const EnumeratedVariable& v) {
        return os << v.name(); }
        //return os << v.name() << ":" << v.type().name(); }
};


//! \brief A named integer variable, suitable for use in formulae and in defining the discrete state.
//! Equality is performed using references by name, so different variables may have the same "name".
//! \details \sa StateSpace, Formula.
template<> struct Variable<Integer> : public DiscreteVariable
{
  public:
    //! \brief Construct a new named integer variable with name \a name. If another variable has the same \c name,
    //! the two are \e not considered to be the same quantity.
    Variable(std::string name) : DiscreteVariable(name) { }
    //! \brief Copy constructer. Returns a C++ variable representing the same named variable.
    Variable(const IntegerVariable& v) : DiscreteVariable(v) { }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const IntegerVariable& v) {
        return os << v.name(); }
        //return os << v.name() << ":Integer"; }
};

//! \brief A named continuous variable.
struct ContinuousVariable : public NamedVariable {
    ContinuousVariable(const std::string& name) : NamedVariable(name) { }
};


template<> struct Quantity<Real> : public ContinuousVariable {
    Quantity(const std::string& name) : ContinuousVariable(name) { }
};

//! \brief A named real variable, suitable for use in formulae and in defining state spaces.
//! Equality is performed using references by name, so different variables may have the same "name".
//! \details \sa StateSpace, Formula.
template<> struct Variable<Real> : public Quantity<Real>
{
  public:
    //! \brief Construct a new named variable with name \a name. If another variable has the same \c name,
    //! the two are \e not considered to be the same quantity.
    Variable(std::string name) : Quantity<Real>(name) { }
    //! \brief Copy constructer. Returns a C++ variable representing the same named variable.
    Variable(const RealVariable& v) : Quantity<Real>(v) { }
    //! \brief A dynamically-allocated ProjectionExpression representing the variable.
    ExpressionInterface* expression(const StateSpace& spc) const;

    //! \brief Create an assignment object representing the variable's defining equation.
    Assignment<RealVariable,RealFormula> operator=(const RealFormula& f) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const RealVariable& v) {
        return os << v.name(); }
        //return os << v.name() << ":Real"; }
};


//! \brief A finite set of variables.
typedef Set<NamedVariable> VariableSet;

/*! \brief A space defined as a list of named variables.
 *  \details \sa Variable
 */
struct StateSpace
{
  public:
    //! \brief The trivial space \f$\R^0\f$.
    StateSpace() : _variables() { }
    StateSpace(const std::set<NamedVariable>& vs) : _variables(vs.begin(),vs.end()) { }
    unsigned int size() const { return _variables.size(); }
    //! \brief The dimension of the space.
    unsigned int dimension() const { return _variables.size(); }
    //! \brief The \a i<sup>th</sup> named variable.
    const NamedVariable& operator[](unsigned int i) const { return _variables.at(i); }
    const NamedVariable& variable(unsigned int i) const { return _variables.at(i); }

    //! \brief The index of the named variable \a v.
    unsigned int index(const NamedVariable& v) const {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(v==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<v<<" is not in the StateSpace "<<*this);
        return _variables.size(); }
    //! \brief Append the named variable \a v to the variables defining the space; ignores if the variable is already in the space.
    StateSpace& insert(const NamedVariable& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(_variables[i]==v) { return *this; } }
        _variables.push_back(v); return *this; }
    //! \brief Adjoins the variables in \a spc.
    StateSpace& adjoin(const StateSpace& spc) {
        for(uint i=0; i!=spc._variables.size(); ++i) { this->insert(spc._variables[i]); } return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    StateSpace& append(const NamedVariable& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(_variables[i]!=v,"Variable "<<v<<" is already a variable of the StateSpace "<<*this);
        }
        _variables.push_back(v); return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    StateSpace& operator,(const NamedVariable& v) { return this->append(v); }
  public:
    //! \brief The space containing the union of the variables.
    friend StateSpace join(const StateSpace&, const StateSpace&);
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const StateSpace&);
  private:
    std::vector<NamedVariable> _variables;
};

inline std::ostream& operator<<(std::ostream& os, const StateSpace& spc) { return os << spc._variables; }

inline StateSpace operator,(const NamedVariable& v1, const NamedVariable& v2) {
    StateSpace r; r,v1,v2; return r; }

inline StateSpace join(const StateSpace& spc1, const StateSpace& spc2) {
    StateSpace r(spc1); r.adjoin(spc2); return r; }



/*! \brief A space defined as a list of named variables.
 *  \details \sa Variable
 */
template<class Var>
struct Space
{
    typedef Var VariableType;
  public:
    //! \brief The trivial space \f$\R^0\f$.
    Space() : _variables() { }
    unsigned int size() const { return _variables.size(); }
    //! \brief The dimension of the space.
    unsigned int dimension() const { return _variables.size(); }
    //! \brief The \a i<sup>th</sup> named variable.
    const VariableType& operator[](unsigned int i) const { return _variables.at(i); }
    const VariableType& variable(unsigned int i) const { return _variables.at(i); }

    //! \brief The index of the named variable \a v.
    unsigned int index(const VariableType& v) const {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(v==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<v<<" is not in the StateSpace "<<*this);
        return _variables.size(); }
    //! \brief Append the named variable \a v to the variables defining the space; ignores if the variable is already in the space.
    Space<VariableType>& insert(const EnumeratedVariable& v) {
        for(uint i=0; i!=_variables.size(); ++i) { if(_variables[i]==v) { return *this; } } _variables.push_back(v); return *this; }
    //! \brief Adjoins the variables in \a spc.
    Space<VariableType>& adjoin(const Space<VariableType>& spc) {
        for(uint i=0; i!=spc._variables.size(); ++i) { this->insert(spc._variables[i]); } return *this; }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space<VariableType>& append(const VariableType& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(_variables[i]!=v,"Variable "<<v<<" is already a variable of the StateSpace "<<*this);
        }
        _variables.push_back(v); return *this; }
  public:
    //! \brief The space containing the union of the variables.
    friend Space<VariableType> join(const Space<VariableType>& spc1, const Space<VariableType>& spc2) {
        Space<VariableType> r(spc1); r.adjoin(spc2); return r; };
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Space<VariableType>& spc) { return os<<spc._variables; }
  private:
    std::vector<VariableType> _variables;
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

//! \brief
struct Valuation {
    void set(const Variable<EnumeratedValue>& v, const String& s) {
        this->enumerated_values.insert(std::make_pair(v,EnumeratedValue(s))); }
    //void set(const Variable<Integer>& v, const Integer& n) { this->integer_values[v]=n; }
    void set(const Variable<Real>& v, const Interval& x) { this->real_values[v]=x; }
    EnumeratedValue get(const Variable<EnumeratedValue>& v) const { return enumerated_values.find(v)->second; }
    //Integer get(const Variable<Integer>& v) const { return integer_values.find(v)->second; }
    Interval get(const Variable<Real>& v) const { return real_values.find(v)->second; }

    template<class T> T get(const NamedVariable& v) const;
  private:
  public:
    std::map<NamedVariable,EnumeratedValue> enumerated_values;
    //std::map<Variable<Integer>,Integer> integer_values;
    std::map<NamedVariable,Interval> real_values;
};

inline std::ostream& operator<<(std::ostream& os, const Valuation& v) {
    return os<<v.enumerated_values<<v.real_values; }


template<class T> inline T _get(const Valuation& p, const NamedVariable& v);
template<> inline EnumeratedValue _get<EnumeratedValue>(const Valuation& p, const NamedVariable& v) {
    return p.enumerated_values.find(v)->second; }
template<> inline Interval _get<Interval>(const Valuation& p, const NamedVariable& v) {
    return p.real_values.find(v)->second; }

template<class T> inline T Valuation::get(const NamedVariable& v) const { return Ariadne::_get<T>(*this,v); }



template<class R> class FormulaInterface;
template<class R> class VariableFormula;
template<class R, class C=R> class ConstantFormula;
template<class R, class O, class A=R> class UnaryFormula;
template<class R, class O, class A1=R, class A2=R> class BinaryFormula;
template<class R> class Formula;
template<class R> std::ostream& operator<<(std::ostream&, const Formula<R>& f);

template<class R>
class FormulaInterface {
  public:
    virtual ~FormulaInterface() { }
    virtual FormulaInterface<R>* clone() const = 0;
    virtual VariableSet arguments() const = 0;
    virtual R evaluate(const Valuation& x) const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
};

template<class R>
inline std::ostream& operator<<(std::ostream& os, const FormulaInterface<R>& f) { return f.write(os); }

template<class R> struct VariableFormula : public FormulaInterface<R> {
    VariableFormula(const Variable<R>& v) : var(v) { }
    virtual VariableFormula<R>* clone() const { return new VariableFormula<R>(*this); }
    virtual VariableSet arguments() const { VariableSet spc; spc.insert(var); return spc; }
    virtual R evaluate(const Valuation& x) const { return x.get(var); }
    virtual std::ostream& write(std::ostream& os) const { return os<<var; }
    Variable<R> var;
};

template<class R,class C> struct ConstantFormula : public FormulaInterface<R> {
    ConstantFormula(const C& c) : val(c) { }
    virtual ConstantFormula<R,C>* clone() const { return new ConstantFormula<R,C>(*this); }
    virtual VariableSet arguments() const { return VariableSet(); }
    virtual R evaluate(const Valuation& x) const { return R(val); }
    virtual std::ostream& write(std::ostream& os) const { return os<<val; }
    C val;
};


template<class R,class O,class A> struct UnaryFormula : public FormulaInterface<R> {
    UnaryFormula(O o, shared_ptr< FormulaInterface<A> > a) : op(o), arg_ptr(a) { }
    UnaryFormula(O o, const FormulaInterface<A>& a) : op(o), arg_ptr(a.clone()) { }
    virtual UnaryFormula<R,O,A>* clone() const { return new UnaryFormula<R,O,A>(*this); }
    virtual VariableSet arguments() const { return arg_ptr->arguments(); }
    virtual R evaluate(const Valuation& x) const { return op(arg_ptr->evaluate(x)); }
    virtual std::ostream& write(std::ostream& os) const { return os<<op<<"("<<*arg_ptr<<")"; }
    O op; shared_ptr< FormulaInterface<A> > arg_ptr;
};

template<class R,class O,class A1,class A2> struct BinaryFormula : public FormulaInterface<R> {
    BinaryFormula(O o, shared_ptr< FormulaInterface<A1> > a1, shared_ptr< FormulaInterface<A2> > a2)
        : op(o), arg1_ptr(a1), arg2_ptr(a2) { }
    BinaryFormula(O o, FormulaInterface<A1> const& a1, FormulaInterface<A2> const& a2)
        : op(o), arg1_ptr(a1.clone()), arg2_ptr(a2.clone()) { }
    BinaryFormula(O o, Formula<A1> a1, Formula<A2> a2)
        : op(o), arg1_ptr(a1.ptr), arg2_ptr(a2.ptr) { }
    virtual BinaryFormula<R,O,A1,A2>* clone() const { return new BinaryFormula<R,O,A1,A2>(*this); }
    virtual VariableSet arguments() const { return join(arg1_ptr->arguments(),arg2_ptr->arguments()); }
    virtual R evaluate(const Valuation& x) const { return op(arg1_ptr->evaluate(x),arg2_ptr->evaluate(x)); }
    virtual std::ostream& write(std::ostream& os) const { return os<<*arg1_ptr<<op<<*arg2_ptr; }
    O op; shared_ptr< FormulaInterface<A1> > arg1_ptr; shared_ptr< FormulaInterface<A2> > arg2_ptr;
};

//! \brief A formula taking values of type \a R.
//! \tparam R The type of mathematical object the formula represents; e.g. String, Integer, Real
//! \sa Variable
template<class R> class Formula {
  public:
    //! \brief The type returned by the formula.
    typedef R result_type;
    //! \brief The constant \a c.
    Formula(const R& c) : ptr(new ConstantFormula<R>(c)) { }
    //! \brief The formula "v" in named variable \a v.
    Formula(const Variable<R>& v) : ptr(new VariableFormula<R>(v)) { }
    Formula(const FormulaInterface<R>& e) : ptr(e.clone()) { }
    explicit Formula(FormulaInterface<R>* p) : ptr(p) { }
    Formula(shared_ptr<FormulaInterface<R> > p) : ptr(p) { }
    //! \brief The variables used in the formula.
    VariableSet arguments() const { return ptr->arguments(); }
    //! \brief Evaluate the formula on the valuation of variables \a x.
    R evaluate(const Valuation& x) const { return ptr->evaluate(x); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<< <>(std::ostream&, const Formula<R>&);
  public:
    shared_ptr< FormulaInterface<R> > ptr;
};

template<class R> inline std::ostream& operator<<(std::ostream& os, const Formula<R>& f) { return f.ptr->write(os); }






template<>
class FormulaInterface<Real> {
  public:
    virtual ~FormulaInterface() { }
    virtual FormulaInterface<Real>* clone() const = 0;
    virtual VariableSet arguments() const = 0;
    virtual ExpressionInterface* expression(const StateSpace& spc) const = 0;
    virtual Interval evaluate(const Valuation& x) const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
};

template<> struct VariableFormula<Real> : public FormulaInterface<Real> {
    VariableFormula(const RealVariable& v) : var(v) { }
    virtual VariableFormula* clone() const { return new VariableFormula(*this); }
    virtual VariableSet arguments() const { VariableSet spc; spc.insert(var); return spc; }
    virtual ExpressionInterface* expression(const StateSpace& spc) const;
    virtual Interval evaluate(const Valuation& x) const { return x.get(var); }
    virtual std::ostream& write(std::ostream& os) const { return os<<var; }
    RealVariable var;
};

template<class C> struct ConstantFormula<Real,C> : public FormulaInterface<Real> {
    ConstantFormula(const C& c) : val(c) { }
    virtual ConstantFormula<Real,C>* clone() const { return new ConstantFormula<Real,C>(*this); }
    virtual VariableSet arguments() const { return VariableSet(); }
    virtual Interval evaluate(const Valuation& x) const { return Interval(val); }
    virtual ExpressionInterface* expression(const StateSpace& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os<<val; }
    C val;
};


template<class Op> struct UnaryFormula<Real,Op> : public FormulaInterface<Real> {
    typedef shared_ptr<FormulaInterface<Real> > FormulaPointer;
    UnaryFormula(Op o, FormulaPointer a) : op(o), arg_ptr(a) { }
    UnaryFormula(Op o, const FormulaInterface<Real>& a) : op(o), arg_ptr(a.clone()) { }
    virtual UnaryFormula<Real,Op>* clone() const { return new UnaryFormula<Real,Op>(*this); }
    virtual VariableSet arguments() const { return arg_ptr->arguments(); }
    virtual Interval evaluate(const Valuation& x) const { return op(arg_ptr->evaluate(x)); }
    virtual ExpressionInterface* expression(const StateSpace& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os<<op<<"("<<*arg_ptr<<")"; }
    //template<class Res> Res evaluate() const { return op(arg.evaluate<Res>()); }
    Op op; FormulaPointer arg_ptr;
};

template<class Op> struct BinaryFormula<Real,Op> : public FormulaInterface<Real> {
    typedef shared_ptr<FormulaInterface<Real> > FormulaPointer;
    BinaryFormula(Op o, FormulaPointer a1, FormulaPointer a2) : op(o), arg1_ptr(a1), arg2_ptr(a2) { }
    virtual BinaryFormula<Real,Op>* clone() const { return new BinaryFormula<Real,Op>(*this); }
    virtual VariableSet arguments() const { return join(arg1_ptr->arguments(),arg2_ptr->arguments()); }
    virtual Interval evaluate(const Valuation& x) const { return op(arg1_ptr->evaluate(x),arg2_ptr->evaluate(x)); }
    virtual ExpressionInterface* expression(const StateSpace& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os<<*arg1_ptr<<op<<*arg2_ptr; }
    Op op; FormulaPointer arg1_ptr; FormulaPointer arg2_ptr;
};

/*! \brief A simple formula in named variables.
 *  \details A %Formula differs from an Expression in three ways.
 *  Firstly, the independent variables are given string names, rather than an integer index.
 *  Secondly, formulae in different variables may be combined; the variables of the resulting formula
 *  are all variables occuring in all formulae.
 *  Thirdly, formulae may be manipulated symbolically, whereas expressions can only be manipulated numerically.
 *  \sa RealVariable, ExpressionInterface
 */
template<> class Formula<Real> {
  public:
    Formula(const int& c) : ptr(new ConstantFormula<Real,int>(c)) { }
    Formula(const double& c) : ptr(new ConstantFormula<Real,double>(c)) { }
    Formula(const Interval& c) : ptr(new ConstantFormula<Real,Interval>(c)) { }
    //! \brief The formula "v" in named variable \a v.
    Formula(const RealVariable& v) : ptr(new VariableFormula<Real>(v)) { }
    Formula(const FormulaInterface<Real>& e) : ptr(e.clone()) { }
    Formula(FormulaInterface<Real>* p) : ptr(p) { }
    Formula(shared_ptr<FormulaInterface<Real> > p) : ptr(p) { }
    //! \brief Create an expression \f$\R^n\rightarrow\R\f$ by assigning the variables of the formula the numbers given by \a spc.
    VariableSet arguments() const { return ptr->arguments(); }
    ExpressionInterface* expression(const StateSpace& spc) { return ptr->expression(spc); }
    Interval evaluate(const Valuation& x) const { return ptr->evaluate(x); }

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Formula<Real>& f);
  public:
    shared_ptr<FormulaInterface<Real> > ptr;
};

inline std::ostream& operator<<(std::ostream& os, const Formula<Real>& f) { return f.ptr->write(os); }

inline Assignment<RealVariable,RealFormula> RealVariable::operator=(const RealFormula& f) const {
    return Assignment<RealVariable,RealFormula>(*this,f); }

class AffineFormula
//    : public FormulaInterface
{
  public:
    explicit AffineFormula() : _b(0), _a() { }
    explicit AffineFormula(const int& c) : _b(c), _a() { }
    explicit AffineFormula(const double& c) : _b(c), _a() { }
    explicit AffineFormula(const Interval& c) : _b(c), _a() { }
    //! \brief The formula "v" in named variable \a v.
    explicit AffineFormula(const RealVariable& v) : _b(0.0), _a() { _a[v]=1.0; }

    AffineFormula* clone() const { return new AffineFormula(*this); }

    //! \brief Create an expression \f$\R^n\rightarrow\R\f$ by assigning the variables of the formula the numbers given by \a spc.
    AffineExpression* expression(const StateSpace& spc) const {
        Vector<Interval> a(spc.dimension());
        for(std::map<RealVariable,Interval>::const_iterator iter=_a.begin(); iter!=_a.end(); ++iter) {
            a[spc.index(iter->first)]+=iter->second;
        }
        return new AffineExpression(a,_b);
    }

    std::ostream& write(std::ostream& os) const;
    friend AffineFormula operator+(const AffineFormula&, const AffineFormula&);
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const AffineFormula& e);
  public:
    Interval _b;
    std::map<RealVariable,Interval> _a;
};

inline AffineFormula& operator*=(AffineFormula& x, const Interval& c) {
    x._b*=c;
    for(std::map<RealVariable,Interval>::iterator iter=x._a.begin(); iter!=x._a.end(); ++iter) {
        iter->second*=c;
    }
    return x;
}

inline AffineFormula& operator+=(AffineFormula& x, const AffineFormula& y) {
    x._b+=y._b;
    for(std::map<RealVariable,Interval>::const_iterator iter=y._a.begin(); iter!=y._a.end(); ++iter) {
        x._a[iter->first]+=iter->second;
    }
    return x;
}

inline AffineFormula& operator-=(AffineFormula& x, const AffineFormula& y) {
    x._b-=y._b;
    for(std::map<RealVariable,Interval>::const_iterator iter=y._a.begin(); iter!=y._a.end(); ++iter) {
        x._a[iter->first]-=iter->second;
    }
    return x;
}

inline AffineFormula operator-(const AffineFormula& x) {
    AffineFormula r; r-=x; return r; }

inline AffineFormula operator+(const AffineFormula& x1, const AffineFormula& x2) {
    AffineFormula r(x1); r+=x2; return r; }

inline AffineFormula operator-(const AffineFormula& x1, const AffineFormula& x2) {
    AffineFormula r(x1); r-=x2; return r; }

inline AffineFormula operator*(const Interval& c, const AffineFormula& x) {
    AffineFormula r(x); r*=c; return r; }

inline AffineFormula operator*(const AffineFormula& x, const Interval& c) {
    AffineFormula r(x); r*=c; return r; }

inline AffineFormula operator/(const AffineFormula& x, const Interval& c) {
    AffineFormula r(x); r*=(1.0/c); return r; }

inline std::ostream& operator<<(std::ostream& os, const AffineFormula& x) {
    if(x._b!=0) { os<<x._b; }
    for(std::map<RealVariable,Interval>::const_iterator iter=x._a.begin(); iter!=x._a.end(); ++iter) {
        if(iter->second!=0) {
            if(iter->second>0) { os << "+"; }
            os << iter->second << "*" << iter->first;
        }
    }
    return os;
}


inline ExpressionInterface* RealVariable::expression(const StateSpace& spc) const {
    return new ProjectionExpression(spc.dimension(),spc.index(*this));
}

template<class C>
inline ExpressionInterface* ConstantFormula<Real,C>::expression(const StateSpace& spc) const {
    return new ConstantExpression(spc.dimension(),Interval(this->val));
}

inline ExpressionInterface* VariableFormula<Real>::expression(const StateSpace& spc) const {
    return new ProjectionExpression(spc.dimension(),spc.index(this->var));
}

template<class Op>
ExpressionInterface* UnaryFormula<Real,Op>::expression(const StateSpace& spc) const {
    return new UnaryExpression<Op>(op,arg_ptr->expression(spc)); }

template<class Op>
ExpressionInterface* BinaryFormula<Real,Op>::expression(const StateSpace& spc) const {
    return new BinaryExpression<Op>(op,arg1_ptr->expression(spc),arg2_ptr->expression(spc)); }




//! \related Formula \brief .
template<class T> inline Formula<T> operator&&(Formula<T> e1, Formula<T> e2) {
    return Formula<T>(new BinaryFormula<T,And,T,T>(And(),e1.ptr,e2.ptr)); }
//! \related Formula \brief .
template<class T> inline Formula<T> operator||(Formula<T> e1, Formula<T> e2) {
    return Formula<T>(new BinaryFormula<T,Or>(Or(),e1.ptr,e2.ptr)); }
//! \related Formula \brief .
template<class T> inline Formula<T> operator!(Formula<T> e) {
    return Formula<T>(new UnaryFormula<T,Not>(Not(),e.ptr)); }

//! \related EnumeratedVariable \brief .
inline
Formula<bool>
operator==(const EnumeratedVariable& lhs, const EnumeratedValue& rhs) {
    return Formula<bool>(new BinaryFormula<bool,Equal,EnumeratedValue,EnumeratedValue>(Equal(),VariableFormula<EnumeratedValue>(lhs),ConstantFormula<EnumeratedValue>(rhs)));
}

//! \related EnumeratedVariable \brief .
inline
Formula<bool>
operator==(const EnumeratedVariable& lhs, const String& rhs) {
    return lhs==EnumeratedValue(rhs);
}

//! \related RealFormula \brief .
inline
Formula<tribool>
operator<=(const RealFormula& lhs, const RealFormula& rhs) {
    return Formula<tribool>(new BinaryFormula<tribool,Less,Real,Real>(Less(),lhs,rhs));
}

//! \related RealFormula \brief .
inline
Formula<tribool>
operator>=(const RealFormula& lhs, const RealFormula& rhs) {
    return Formula<tribool>(new BinaryFormula<tribool,Gtr,Real,Real>(Gtr(),lhs,rhs));
}

//! \related RealFormula \brief .
inline RealFormula operator+(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Real,Add>(Add(),e1.ptr,e2.ptr)); }
//! \related RealFormula \brief .
inline RealFormula operator-(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Real,Sub>(Sub(),e1.ptr,e2.ptr)); }
//! \related RealFormula \brief .
inline RealFormula operator*(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Real,Mul>(Mul(),e1.ptr,e2.ptr)); }
//! \related RealFormula \brief .
inline RealFormula operator/(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Real,Div>(Div(),e1.ptr,e2.ptr)); }


//! \related RealFormula \brief .
inline RealFormula exp(RealFormula e) {
    return RealFormula(new UnaryFormula<Real,Exp>(Exp(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula log(RealFormula e) {
    return RealFormula(new UnaryFormula<Real,Log>(Log(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula sin(RealFormula e) {
    return RealFormula(new UnaryFormula<Real,Sin>(Sin(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula cos(RealFormula e) {
    return RealFormula(new UnaryFormula<Real,Cos>(Cos(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula tan(RealFormula e) {
    return RealFormula(new UnaryFormula<Real,Tan>(Tan(),e.ptr)); }


inline bool operator==(const ContinuousPredicate& pred, bool truth) {
    ConstantFormula<tribool> const* f=dynamic_cast<ConstantFormula<tribool> const*>(&*pred.ptr);
    return (f && f->val==false);
}


template<class VAR> struct Next;

template<>
struct Next<EnumeratedVariable>
{
    EnumeratedVariable base;
    typedef Assignment<Next<EnumeratedVariable>,EnumeratedFormula> EnumeratedAssignment;
    friend Next<EnumeratedVariable> next(const EnumeratedVariable& v);
    EnumeratedAssignment operator=(const String& str) { return EnumeratedAssignment(*this,EnumeratedValue(str)); }
    EnumeratedAssignment operator=(const EnumeratedValue& val) { return EnumeratedAssignment(*this,val); }
    EnumeratedAssignment operator=(const EnumeratedVariable& var) { return EnumeratedAssignment(*this,var); }
    EnumeratedAssignment operator=(const EnumeratedFormula& fn) { return EnumeratedAssignment(*this,fn); }
    friend std::ostream& operator<<(std::ostream& os, const Next<EnumeratedVariable>& self) {
        return os << "next(" << self.base.name() << ")"; }
  private:
    Next(const EnumeratedVariable& v) : base(v) { }
};

template<>
struct Next<RealVariable>
{
    RealVariable base;
    typedef Assignment<Next<RealVariable>,RealFormula> RealAssignment;
    friend Next<RealVariable> next(const RealVariable& v);
    RealAssignment operator=(const double& x) { return RealAssignment(*this,RealFormula(x)); }
    RealAssignment operator=(const RealVariable& var) { return RealAssignment(*this,RealFormula(var)); }
    RealAssignment operator=(const RealFormula& fn) { return RealAssignment(*this,fn); }
    friend std::ostream& operator<<(std::ostream& os, const Next<RealVariable>& self) {
        return os << "next(" << self.base.name() << ")"; }
  private:
    Next(const RealVariable& v) : base(v) { }
};

inline Next<EnumeratedVariable> next(const EnumeratedVariable& v) { return Next<EnumeratedVariable>(v); }
inline Next<RealVariable> next(const RealVariable& v) { return Next<RealVariable>(v); }

template<class T> class DottedVariable;

template<class T>
struct DottedVariable : public Quantity<T>
{
    Variable<T> base;
    Assignment<DottedVariable<T>,Formula<T> > operator=(const Formula<T>& f) {
        return Assignment<DottedVariable<T>,Formula<T> >(*this,f); }
    Assignment<DottedVariable<T>,Formula<T> > operator=(const Variable<T>& v) {
        return Assignment<DottedVariable<T>,Formula<T> >(*this,Formula<T>(v)); }
    friend std::ostream& operator<<(std::ostream& os, const DottedVariable<T>& self) {
        return os << "dot(" << self.base.name() << ")"; }
    explicit DottedVariable(const Variable<T>& v) : Quantity<T>("dot("+v.name()+")"), base(v) { }
};

inline DottedVariable<Real> dot(const Variable<Real>& v) { return DottedVariable<Real>(v); }

typedef Assignment<EnumeratedVariable,EnumeratedFormula> EnumeratedAssignment;
typedef Assignment<Next<EnumeratedVariable>,EnumeratedFormula> EnumeratedUpdate;
typedef Assignment<RealVariable,RealFormula> RealAssignment;
typedef Assignment<Next<RealVariable>,RealFormula> RealUpdate;
typedef Assignment<DottedVariable<Real>,Formula<Real> > RealDynamic;







} // namespace Ariadne

#endif // ARIADNE_FORMULA_H
