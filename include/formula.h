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
#include "expression_interface.h"
#include "function_interface.h"
#include "function.h"

#include "macros.h"
#include "pointer.h"

#include "polynomial.h"

namespace Ariadne {





typedef std::string String;

class EnumeratedType;
class EnumeratedVariable;
class DiscreteValue;
class DiscreteAssignment;
class DiscreteUpdate;
class DiscretePredicate;
class DiscreteFormula;

class ArrayBuilderTag {};
static ArrayBuilderTag build_array;

template<class T>
inline std::vector<T> operator,(ArrayBuilderTag, const T& t) {
    std::vector<T> v; v.push_back(t); return v; }
template<class T>
inline std::vector<T> operator,(const T& t1, const T& t2) {
    std::vector<T> v; v.push_back(t1); v.push_back(t2); return v; }
template<class T>
inline std::vector<T> operator,(const std::vector<T>& v, const T& t) {
    std::vector<T> r(v); r.push_back(t); return r; }

inline std::vector<std::string> operator,(const std::string& s1, const char* s2) {
    std::vector<std::string> v; v.push_back(std::string(s1)); v.push_back(std::string(s2)); return v; }
inline std::vector<std::string> operator,(const std::vector<std::string>& v1, const char* s2) {
    std::vector<std::string> r(v1); r.push_back(std::string(s2)); return r; }

inline std::vector<std::string> operator,(const ArrayBuilderTag& a, const char* s) {
    std::vector<std::string> r; r.push_back(std::string(s)); return r; }


class SetBuilderTag {};
static SetBuilderTag build_set;
template<class T>
inline std::set<T> operator,(SetBuilderTag, const T& t) {
    std::set<T> s; s.insert(t); return s; }
template<class T>
inline std::set<T> operator,(const std::set<T>& s, const T& t) {
    std::set<T> r(s); r.insert(t); return r; }



/*! \brief An ASCII string. */
typedef std::string String;

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
struct Variable {
  protected:
    //! \brief Create an internal enumerated variable with the given \a name.
    Variable(const String& name) : _name_ptr(new std::string(name)) { }
  public:
    //! \brief The internal name of the variable.
    const String& name() const { return *this->_name_ptr; }
  private:
    shared_ptr<std::string> _name_ptr;
};


//! \brief A named variable which can take values of a given enumerated type.
//! \details \sa EnumeratedType, EnumeratedVariable, DiscretePredicate
struct EnumeratedVariable : public Variable {
  public:
    //! \brief Create an internal enumerated variable with the given \a name and \a type.
    EnumeratedVariable(const String& name, const EnumeratedType& type) : Variable(name), _type_ptr(&type) { }
    //! \brief The internal name of the variable.
    const String& name() const { return this->Variable::name(); }
    //! \brief The internal type of the variable.
    const EnumeratedType& type() const { return *this->_type_ptr; }
    //! \brief Tests if two %EnumeratedVariable objects represent the same internal variable.
    bool operator==(const EnumeratedVariable& other) const { return this->name()==other.name(); }
    //! \brief Tests if two %EnumeratedVariable objects represent different internal variables.
    bool operator!=(const EnumeratedVariable& other) const { return this->name()!=other.name(); }
    //! \brief Makes a discrete assignment object.
    DiscreteAssignment operator=(const DiscreteValue& val) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const EnumeratedVariable& self) {
        return os << self.name() << ":" << self.type().name(); }
  private:
    const EnumeratedType* _type_ptr;
};


struct NextEnumeratedVariable {
    explicit NextEnumeratedVariable(const EnumeratedVariable& v) : _base(v) { }
    DiscreteUpdate operator=(const DiscreteValue& val) const;
    DiscreteUpdate operator=(const EnumeratedVariable& var) const;
    DiscreteUpdate operator=(const std::string& val) const;
    EnumeratedVariable base() const { return this->_base; }
  private:
    EnumeratedVariable _base;
};

inline
NextEnumeratedVariable next(const EnumeratedVariable& v) {
    return NextEnumeratedVariable(v);
}


typedef EnumeratedType DiscreteType;
typedef EnumeratedVariable DiscreteVariable;
typedef NextEnumeratedVariable NextDiscreteVariable;


//! \brief A predicate over some discrete variables.
//! \details \sa DiscreteValue, EnumeratedVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscretePredicate {
    friend class DiscreteConstraint;
  public:
    DiscretePredicate() { }
    DiscretePredicate(bool val) { }
    // //! \brief Evaluate the predicate at the given state.
    // bool operator() (const DiscreteValuation&) const;

    //! \brief .
    friend DiscretePredicate operator==(const EnumeratedVariable&, const String&);
    //! \brief .
    friend DiscretePredicate operator!=(const EnumeratedVariable&, const String&);
    //! \brief .
    friend DiscretePredicate operator||(const DiscretePredicate&, const DiscretePredicate&);
    //! \brief .
    friend DiscretePredicate operator&&(const DiscretePredicate&, const DiscretePredicate&);
    //! \brief .
    friend DiscretePredicate operator!(const DiscretePredicate&);

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const DiscretePredicate&);
  private:
    std::vector< std::vector< std::pair<EnumeratedVariable, String> > > _cnf;
};

//! \brief A predicate over some discrete variables.
//! \details \sa DiscreteValue, EnumeratedVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscreteValue {
  public:
    //! \brief Construct a discrete value based on the string \a str.
    DiscreteValue(const std::string& str) : _value(str) { }
    //! \brief The string representation of the value.
    const std::string& str() const { return this->_value; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteValue& self) {
        return os << self._value; }
  private:
    std::string _value;
};

/*
inline
bool DiscretePredicate::operator() (const DiscreteValuation& state) const {
    for(uint i=0; i!=_cnf.size(); ++i) {
        bool value=false;
        for(uint j=0; j!=_cnf[i].size(); ++j) {
            if(state[_cnf[i][j].first]==_cnf[i][j].second) {
                value=true; break;
            }
        }
        if(value==false) { return false; }
    }
    return true;
}
*/

inline
DiscretePredicate operator==(const EnumeratedVariable& var, const String& val) {
    std::vector< std::pair<EnumeratedVariable, String> > clause(1u,std::make_pair(var,val));
    DiscretePredicate result; result._cnf.push_back(clause); return result;
}

inline
DiscretePredicate operator!=(const EnumeratedVariable&, const String&);

inline
DiscretePredicate operator||(const DiscretePredicate& dp1, const DiscretePredicate& dp2) {
    assert(dp1._cnf.size()==1);
    assert(dp1._cnf.size()==1);
    DiscretePredicate result;
    result._cnf.push_back(dp1._cnf[0]);
    for(uint j=0; j!=dp2._cnf[0].size(); ++j) {
        result._cnf[0].push_back(dp2._cnf[0][j]); }
    return result;
}

inline
DiscretePredicate operator&&(const DiscretePredicate& dp1, const DiscretePredicate& dp2) {
    DiscretePredicate result(dp1);
    for(uint i=0; i!=dp2._cnf.size(); ++i) {
        result._cnf.push_back(dp2._cnf[i]); }
    return result;
}

inline
DiscretePredicate operator!(const DiscretePredicate&);

inline
std::ostream& operator<<(std::ostream& os, const DiscretePredicate& self)
{
    if(self._cnf.empty()) { return os<<"true"; }
    for(uint i=0; i!=self._cnf.size(); ++i) {
        if(i!=0) { os << " & "; }
        for(uint j=0; j!=self._cnf[i].size(); ++j) {
            os << (j==0?"(":" | ") << self._cnf[i][j].first.name()<<"=="<<self._cnf[i][j].second;
        }
        os << ")";
    }
    return os;
}


struct DiscreteFormula {
    DiscreteFormula(const DiscreteValue& val) : _str(val.str()) { }
    DiscreteFormula(const DiscreteVariable& var) : _str(var.name()) { }
    friend std::ostream& operator<<(std::ostream& os, const DiscreteFormula& self) {
        return os << self._str; }
  private:
    std::string _str;
};

/*

//! \brief A state in an integer lattice.
//! \details \sa DiscreteConstraint
struct DiscreteState {
  public:
    DiscreteState(const std::vector<int>& values) : _values(values) { }

    uint size() const { return this->_values.size(); }
    uint operator[](uint i) const { return this->_values[i]; }
    bool operator==(const DiscreteState& other) const {
        assert(this->size()==other.size());
        for(uint i=0; i!=_values.size(); ++i) {
            if(this->_values[i]!=other._values[i]) { return false; } }
        return true; }
    bool operator!=(const DiscreteState& other) const {
        return  !(*this==other); }
    bool operator<(const DiscreteState& other) const {
        assert(this->size()==other.size());
        for(uint i=0; i!=_values.size(); ++i) {
            if(this->_values[i]!=other._values[i]) { return this->_values[i]<other._values[i]; } }
        return false; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteState& self) {
        for(uint i=0; i!=self._values.size(); ++i) { os << (i==0?'(':',') << "q["<<i<<"]="<<self._values[i]; } return os << ')'; }
  private:
    std::vector<int> _values;
};


//! \brief A finite state space, consisting of a list of enumerated variables.
//! \details \sa EnumeratedType, EnumeratedVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscreteSpace {
    friend class DiscreteState;
  public:
    //! \brief Create a space with named variables given by \a variables.
    DiscreteSpace(const std::vector<EnumeratedVariable>& variables) : _variables(variables) {
        for(uint i=0; i!=this->_variables.size(); ++i) {
            for(uint j=i+1; j!=this->_variables.size(); ++j) {
                assert(this->_variables[i].name()!=this->_variables[j].name()); } } }
    //! \brief The number of variables needed to describe a state in the space.
    uint size() const { return _variables.size(); }
    //! \brief The \a i<sup>th</sup> variable.
    const EnumeratedVariable& variable(uint i) const { return this->_variables.at(i); }
    //! \brief The index of the variable with name \a str.
    uint index(const String& str) const { for(uint i=0; i!=_variables.size(); ++i) {
        if(_variables[i].name()==str) { return i; } }
        assert(false);
    }
    //! \brief The index of the variable  \a var.
    uint index(const EnumeratedVariable& var) const { return this->index(var.name()); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteSpace& self) {
        for(uint i=0; i!=self._variables.size(); ++i) { os << (i==0?'(':',') << self._variables[i].name()<<":"<<self._variables[i].type().name(); } return os << ')'; }
  private:
    std::vector<EnumeratedVariable> _variables;
};



//! \brief A valuation of the enumerated variables of some discrete space.
//! \details \sa EnumeratedType, EnumeratedVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscreteValuation {
  public:
    //! \brief Assign \a values to the variables in \a space.
    DiscreteValuation(const DiscreteSpace& space, const std::vector<String>& values)
        : _space_ptr(&space), _values(values.size()) {
        assert(space.size()==values.size());
        for(uint i=0; i!=values.size(); ++i) {
            this->_values[i]=space.variable(i).type().index(values[i]); }
    }

    //! \brief The underlying discrete space.
    const DiscreteSpace& space() const { return *this->_space_ptr; }
    //! \brief The value of the \a i<sup>th</sup> variable.
    String operator[](uint i) const {
        return this->_space_ptr->variable(i).type().value(this->_values[i]); }
    //! \brief The value of the internal named variable \a var.
    String operator[](const EnumeratedVariable& var) const {
        for(uint i=0; i!=this->_space_ptr->size(); ++i) {
            if(this->_space_ptr->variable(i)==var) {
                return (*this)[i]; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<var<<" is not a member of space "<<this->space()); }
    //! \brief Test equality of two valuations on the same space.
    bool operator==(const DiscreteValuation& other) const {
        assert(this->_space_ptr==other._space_ptr);
        for(uint i=0; i!=_values.size(); ++i) {
            if(this->_values[i]!=other._values[i]) { return false; } }
        return true; }
    //! \brief Test inequality of two valuations on the same space.
    bool operator!=(const DiscreteValuation& other) const {
        return  !(*this==other); }
    //! \brief Lexicographic ordering of valuations on the same space.
    bool operator<(const DiscreteValuation& other) const {
        assert(this->_space_ptr==other._space_ptr);
        for(uint i=0; i!=_values.size(); ++i) {
            if(this->_values[i]!=other._values[i]) { return this->_values[i]<other._values[i]; } }
        return false; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteValuation& self) {
        for(uint i=0; i!=self._values.size(); ++i) { os << (i==0?'(':',') << self.space().variable(i).name()<<"="<<self[i]; } return os << ')'; }
  private:
    const DiscreteSpace* _space_ptr;
    std::vector<int> _values;
};

//! \brief A constraint over integer variables.
//! \details \sa DiscreteState
struct DiscreteConstraint {
  public:
    DiscreteConstraint(const DiscretePredicate& pred, const DiscreteSpace& spc);
    //! \brief Evaluate the predicate at the given state.
    bool operator() (const DiscreteState&) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const DiscreteConstraint&);
  private:
    std::vector< std::vector< std::pair<uint,uint> > > _cnf;
};

inline
DiscreteConstraint::DiscreteConstraint(const DiscretePredicate& pred, const DiscreteSpace& spc) {
    for(uint i=0; i!=pred._cnf.size(); ++i) {
        _cnf.push_back(std::vector<std::pair<uint,uint> >());
        for(uint j=0; j!=pred._cnf[i].size(); ++j) {
            EnumeratedVariable const& var=pred._cnf[i][j].first;
            String const& val=pred._cnf[i][j].second;
            this->_cnf[i].push_back(std::pair<uint,uint>(spc.index(var),var.type().index(val)));
        }
    }
}

inline
bool DiscreteConstraint::operator() (const DiscreteState& state) const {
    for(uint i=0; i!=_cnf.size(); ++i) {
        bool value=false;
        for(uint j=0; j!=_cnf[i].size(); ++j) {
            if(state[_cnf[i][j].first]==_cnf[i][j].second) {
                value=true; break;
            }
        }
        if(value==false) { return false; }
    }
    return true;
}

inline
std::ostream& operator<<(std::ostream& os, const DiscreteConstraint& self)
{
    for(uint i=0; i!=self._cnf.size(); ++i) {
        if(i!=0) { os << " & "; }
        for(uint j=0; j!=self._cnf[i].size(); ++j) {
            os << (j==0?"(":" | ") << "q["<<self._cnf[i][j].first<<"]=="<<self._cnf[i][j].second;
        }
        os << ")";
    }
    return os;
}

*/

struct DiscreteAssignment {
    DiscreteVariable lhs;
    DiscreteFormula rhs;
};

struct DiscreteUpdate {
    DiscreteVariable lhs;
    DiscreteFormula rhs;
};

inline
DiscreteAssignment EnumeratedVariable::operator=(const DiscreteValue& val) const {
    DiscreteAssignment a={*this,DiscreteFormula(val)}; return a; }

inline
DiscreteUpdate NextEnumeratedVariable::operator=(const DiscreteVariable& var) const {
    DiscreteUpdate a={this->base(),DiscreteFormula(var)}; return a; }

inline
DiscreteUpdate NextEnumeratedVariable::operator=(const DiscreteValue& val) const {
    DiscreteUpdate a={this->base(),DiscreteFormula(val)}; return a; }

inline
DiscreteUpdate NextEnumeratedVariable::operator=(const std::string& str) const {
    DiscreteUpdate a={this->base(),DiscreteFormula(DiscreteValue(str))}; return a; }


//! \brief A discrete event.
struct Event {
  public:
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

    friend EventSet operator!(const Event& e);
    friend EventSet operator!(const EventSet& s);

    friend std::ostream& operator<<(std::ostream& os, const EventSet& self) {
        if(self._is_complement) { os<<"!"; } return os << self._events; }
  private:
    std::set<Event> _events;
    bool _is_complement;
};

inline
EventSet operator!(const Event& e) {
    return !(EventSet(e)); }

inline
EventSet operator!(const EventSet& s) {
    EventSet r; r._events=s._events; r._is_complement=!s._is_complement; return r; }

inline
EventSet operator!(const std::set<Event>& s) {
    return !EventSet(s); }


class RealVariable;
class DottedVariable;
class NextVariable;
class Formula;
class Assignment;
class Comparison;
class Predicate;

struct GtrZero {}; struct LessZero {}; struct Gtr { }; struct Less { };

struct Add { template<class T> T operator()(const T& a1, const T& a2) const { return a1+a2; } };
struct Sub { template<class T> T operator()(const T& a1, const T& a2) const { return a1-a2; } };
struct Mul { template<class T> T operator()(const T& a1, const T& a2) const { return a1*a2; } };
struct Div { template<class T> T operator()(const T& a1, const T& a2) const { return a1/a2; } };

//struct Pow { template<class T, class N> T operator()(const T& a, const N& n) const { return Ariadne::pow(a,n); } };
struct Pow { Pow(int n) : _n(n) { } template<class T> T operator()(const T& a) const { return Ariadne::pow(a,_n); } private: int _n; };

struct Sqr { template<class T> T operator()(const T& a) const { return Ariadne::sqr(a); } };
struct Sqrt { template<class T> T operator()(const T& a) const { return Ariadne::sqrt(a); } };

struct Exp { template<class T> T operator()(const T& a) const { return Ariadne::exp(a); } };
struct Log { template<class T> T operator()(const T& a) const { return Ariadne::log(a); } };
struct Sin { template<class T> T operator()(const T& a) const { return Ariadne::sin(a); } };
struct Cos { template<class T> T operator()(const T& a) const { return Ariadne::cos(a); } };
struct Tan { template<class T> T operator()(const T& a) const { return Ariadne::tan(a); } };

class Variable;
class Space;
class FormulaInterface;

typedef shared_ptr<FormulaInterface> FormulaPointer;

struct FormulaInterface {
    virtual ~FormulaInterface() { }
    virtual FormulaInterface* clone() const = 0;
    virtual ExpressionInterface* expression(const Space& spc) const = 0;
    template<class Res> Res evaluate() const;
    virtual std::ostream& write(std::ostream& os) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const FormulaInterface& f) { return f.write(os); }



class Assignment;
class Dynamic;
class Update;

/*! \brief A named variable, suitable for use in formulae and in defining state spaces.
 *  Equality is performed using references by name, so different variables may have the same "name".
 *  \details \sa Space, Formula.
 */
struct RealVariable
{
  public:
    //! \brief Construct a new named variable with name \a name. If another variable has the same \c name,
    //! the two are \e not considered to be the same quantity.
    RealVariable(std::string name) : _name_ptr(new std::string(name)) { }
    //! \brief Copy constructer. Returns a C++ variable representing the same named variable.
    RealVariable(const RealVariable& v) : _name_ptr(v._name_ptr) { }
    //! \brief The name of the variable.
    const std::string& name() const { return *this->_name_ptr; }
    //! \brief A dynamically-allocated ProjectionExpression representing the variable.
    ExpressionInterface* expression(const Space& spc) const;
    //! \brief Equality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator==(const RealVariable& v) const { return this->_name_ptr==v._name_ptr; }
    //! \brief Inequality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator!=(const RealVariable& v) const { return !(*this==v); }

    bool operator<(const RealVariable& v) const { return *this->_name_ptr<*v._name_ptr; }

    Assignment operator=(const Formula& f) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const RealVariable& v);
  private:
    shared_ptr<const std::string> _name_ptr;
};


inline std::ostream& operator<<(std::ostream& os, const RealVariable& v) { return os << v.name(); }


struct DottedVariable
{
  public:
    friend DottedVariable dot(const RealVariable& v);
    Dynamic operator=(const Formula& f) const;
    RealVariable base() const { return _base_variable; }
  private:
    DottedVariable(const RealVariable& v) : _base_variable(v) { }
    RealVariable _base_variable;
};

inline DottedVariable dot(const RealVariable& v) {
    return DottedVariable(v);
}


struct NextVariable
{
  public:
    friend NextVariable next(const RealVariable& v);
    Update operator=(const RealVariable& v) const;
    Update operator=(const Formula& f) const;
    RealVariable base() const { return _base_variable; }
  private:
    explicit NextVariable(const RealVariable& v) : _base_variable(v) { }
    RealVariable _base_variable;
};

inline NextVariable next(const RealVariable& v) {
    return NextVariable(v);
}

/*! \brief A space defined as a list of named variables.
 *  \details \sa RealVariable
 */
struct Space
{
  public:
    //! \brief The trivial space \f$\R^0\f$.
    Space() : _variables() { }
    unsigned int size() const { return _variables.size(); }
    //! \brief The dimension of the space.
    unsigned int dimension() const { return _variables.size(); }
    //! \brief The \a i<sup>th</sup> named variable.
    const RealVariable& operator[](unsigned int i) const { return _variables.at(i); }
    const RealVariable& variable(unsigned int i) const { return _variables.at(i); }

    //! \brief The index of the named variable \a v.
    unsigned int index(const RealVariable& v) const {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(v==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"RealVariable "<<v<<" is not in the Space "<<*this);
        return _variables.size(); }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space& operator,(const RealVariable& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(_variables[i]!=v,"RealVariable "<<v<<" is already a variable of the Space "<<*this);
        }
        _variables.push_back(v); return *this; }
  public:
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const Space&);
  private:
    std::vector<RealVariable> _variables;
};

inline std::ostream& operator<<(std::ostream& os, const Space& spc) { return os << spc._variables; }

inline Space operator,(const RealVariable& v1, const RealVariable& v2) {
    Space r; r,v1,v2; return r; }



template<class C> struct ConstantFormula : public FormulaInterface {
  public:
    ConstantFormula(const C& c) : _value(c) { }
    virtual ConstantFormula<C>* clone() const { return new ConstantFormula<C>(*this); }
    virtual ExpressionInterface* expression(const Space& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os<<_value; }
  private:
    C _value;
};

struct VariableFormula : public FormulaInterface {
  public:
    VariableFormula(const RealVariable& var) : _var(var) { }
    virtual VariableFormula* clone() const { return new VariableFormula(*this); }
    virtual ExpressionInterface* expression(const Space& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os<<_var; }
  private:
    RealVariable _var;
};


template<class Op> struct UnaryFormula : public FormulaInterface {
    UnaryFormula(Op o, FormulaPointer a) : op(o), arg_ptr(a) { }
    UnaryFormula(Op o, const FormulaInterface& a) : op(o), arg_ptr(a.clone()) { }
    virtual UnaryFormula<Op>* clone() const { return new UnaryFormula<Op>(*this); }
    virtual ExpressionInterface* expression(const Space& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os<<op<<"("<<*arg_ptr<<")"; }
    //template<class Res> Res evaluate() const { return op(arg.evaluate<Res>()); }
    Op op; FormulaPointer arg_ptr;
};

template<class Op> struct BinaryFormula : public FormulaInterface {
    BinaryFormula(Op o, FormulaPointer a1, FormulaPointer a2) : op(o), arg1_ptr(a1), arg2_ptr(a2) { }
    virtual BinaryFormula<Op>* clone() const { return new BinaryFormula<Op>(*this); }
    virtual ExpressionInterface* expression(const Space& spc) const;
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
class Formula {
  public:
    Formula(const int& c) : ptr(new ConstantFormula<int>(c)) { }
    Formula(const double& c) : ptr(new ConstantFormula<double>(c)) { }
    Formula(const Interval& c) : ptr(new ConstantFormula<Interval>(c)) { }
    //! \brief The formula "v" in named variable \a v.
    Formula(const RealVariable& v) : ptr(new VariableFormula(v)) { }
    Formula(FormulaInterface* p) : ptr(p) { }
    Formula(shared_ptr<FormulaInterface> p) : ptr(p) { }
    //! \brief Create an expression \f$\R^n\rightarrow\R\f$ by assigning the variables of the formula the numbers given by \a spc.
    ExpressionInterface* expression(const Space& spc) { return ptr->expression(spc); }

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Formula& e);
  public:
    FormulaPointer ptr;
};

inline std::ostream& operator<<(std::ostream& os, const Formula& e) { return e.ptr->write(os); }


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
    AffineExpression* expression(const Space& spc) const {
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


inline ExpressionInterface* RealVariable::expression(const Space& spc) const {
    return new ProjectionExpression(spc.dimension(),spc.index(*this));
}

template<class C>
inline ExpressionInterface* ConstantFormula<C>::expression(const Space& spc) const {
    return new ConstantExpression(spc.dimension(),Interval(this->_value));
}

inline ExpressionInterface* VariableFormula::expression(const Space& spc) const {
    return new ProjectionExpression(spc.dimension(),spc.index(this->_var));
}

template<class Op>
ExpressionInterface* UnaryFormula<Op>::expression(const Space& spc) const {
    return new UnaryExpression<Op>(op,arg_ptr->expression(spc)); }

template<class Op>
ExpressionInterface* BinaryFormula<Op>::expression(const Space& spc) const {
    return new BinaryExpression<Op>(op,arg1_ptr->expression(spc),arg2_ptr->expression(spc)); }


inline FormulaPointer make_formula_pointer(const Float& c) { return FormulaPointer(new ConstantFormula<Float>(c)); }
inline FormulaPointer make_formula_pointer(const Interval& c) { return FormulaPointer(new ConstantFormula<Interval>(c)); }

inline Formula operator+(Formula e1, Formula e2) {
    return Formula(new BinaryFormula<Add>(Add(),e1.ptr,e2.ptr)); }
inline Formula operator-(Formula e1, Formula e2) {
    return Formula(new BinaryFormula<Sub>(Sub(),e1.ptr,e2.ptr)); }
inline Formula operator*(Formula e1, Formula e2) {
    return Formula(new BinaryFormula<Mul>(Mul(),e1.ptr,e2.ptr)); }
inline Formula operator/(Formula e1, Formula e2) {
    return Formula(new BinaryFormula<Div>(Div(),e1.ptr,e2.ptr)); }


/*
FormulaPointer operator+(FormulaPointer e1, double e2) { return e1+make_formula_pointer(e2); }
FormulaPointer operator-(FormulaPointer e1, double e2) { return e1-make_formula_pointer(e2); }
FormulaPointer operator*(FormulaPointer e1, double e2) { return e1*make_formula_pointer(e2); }
FormulaPointer operator/(FormulaPointer e1, double e2) { return e1/make_formula_pointer(e2); }

FormulaPointer operator+(FormulaPointer e1, Interval e2) { return e1+make_formula_pointer(e2); }
FormulaPointer operator-(FormulaPointer e1, Interval e2) { return e1-make_formula_pointer(e2); }
FormulaPointer operator*(FormulaPointer e1, Interval e2) { return e1*make_formula_pointer(e2); }
FormulaPointer operator/(FormulaPointer e1, Interval e2) { return e1/make_formula_pointer(e2); }

FormulaPointer operator+(double e1, FormulaPointer e2) { return make_formula_pointer(e1)+e2; }
FormulaPointer operator-(double e1, FormulaPointer e2) { return make_formula_pointer(e1)-e2; }
FormulaPointer operator*(double e1, FormulaPointer e2) { return make_formula_pointer(e1)*e2; }
FormulaPointer operator/(double e1, FormulaPointer e2) { return make_formula_pointer(e1)/e2; }

FormulaPointer operator+(Interval e1, FormulaPointer e2) { return make_formula_pointer(e1)+e2; }
FormulaPointer operator-(Interval e1, FormulaPointer e2) { return make_formula_pointer(e1)-e2; }
FormulaPointer operator*(Interval e1, FormulaPointer e2) { return make_formula_pointer(e1)*e2; }
FormulaPointer operator/(Interval e1, FormulaPointer e2) { return make_formula_pointer(e1)/e2; }
*/


inline Formula exp(Formula e) {
    return Formula(new UnaryFormula<Exp>(Exp(),e.ptr)); }
inline Formula log(Formula e) {
    return Formula(new UnaryFormula<Log>(Log(),e.ptr)); }
inline Formula sin(Formula e) {
    return Formula(new UnaryFormula<Sin>(Sin(),e.ptr)); }
inline Formula cos(Formula e) {
    return Formula(new UnaryFormula<Cos>(Cos(),e.ptr)); }
inline Formula tan(Formula e) {
    return Formula(new UnaryFormula<Tan>(Tan(),e.ptr)); }


/*
class Formula
{
  public:
    Formula(const FormulaInterface* f_ptr) : _ptr(const_cast<FormulaInterface*>(f_ptr)) { }
    Formula(const RealVariable& v) : _ptr(new VariableFormula(v)) { }
    shared_ptr<const ExpressionInterface> expression(const Space& spc) {
        return shared_ptr<const ExpressionInterface>(_ptr->expression(spc)); }
  private:
  public:
    Formula _ptr;
};

//inline std::ostream& operator<<(std::ostream& os, const Formula& f) { os<<f._ptr<<":"<<std::flush; return f._ptr->write(os); }
inline void dump(const Formula fptr) { fptr->write(std::cout); }


inline Formula make_constant(const int& c) { return Formula(new ConstantFormula<int>(c)); }
inline Formula make_constant(const double& c) { return Formula(new ConstantFormula<double>(c)); }
inline Formula make_constant(const Interval& c) { return Formula(new ConstantFormula<Interval>(c)); }

inline Formula operator+(const Formula& x1, const Formula& x2) {
    return Formula(new BinaryFormula<Add>(Add(),x1._ptr,x2._ptr)); }
inline Formula operator-(const Formula& x1, const Formula& x2) {
    return Formula(new BinaryFormula<Sub>(Sub(),x1._ptr,x2._ptr)); }
inline Formula operator*(const Formula& x1, const Formula& x2) {
    return Formula(new BinaryFormula<Mul>(Mul(),x1._ptr,x2._ptr)); }
inline Formula operator/(const Formula& x1, const Formula& x2) {
    return Formula(new BinaryFormula<Div>(Div(),x1._ptr,x2._ptr)); }

inline ConstantFormula<Interval> formula(const Interval& c) { return ConstantFormula<Interval>(c); }
inline VariableFormula<Interval> formula(const RealVariable& x) { return VariableFormula(c); }
inline FormulaInterface& formula(const FormulaInterface& e) { return e; }

template<class Op> make_binary_formula(Op op, const FormulaInterface& x1, const FormulaInterface& x2) {
    return BinaryFormula<Op>(op,
inline BinaryFormula<Add> operator+(const FormulaInterface& x1, const FormulaInterface& x2) {
    return make_binary_formula(Add(),x1,x2); }
inline BinaryFormula<Sub> operator-(const FormulaInterface& x1, const FormulaInterface& x2) {
    return make_binary_formula(Sub(),x1,x2); }
inline BinaryFormula<Mul> operator*(const FormulaInterface& x1, const FormulaInterface& x2) {
    return make_binary_formula(Mul(),x1,x2); }
inline BinaryFormula<Div> operator/(const FormulaInterface& x1, const FormulaInterface& x2) {
    return make_binary_formula(Div(),x1,x2); }
inline operator-(const FormulaInterface& x1, const FormulaInterface& x2) { return formula(x1)-formula(x2); }
inline operator*(const FormulaInterface& x1, const FormulaInterface& x2) { return formula(x1)*formula(x2); }
inline operator/(const FormulaInterface& x1, const FormulaInterface& x2) { return formula(x1)/formula(x2); }

inline BinaryFormula<Add> operator+(const RealVariable& x1, const RealVariable& x2) { return formula(x1)+formula(x2); }
inline BinaryFormula<Sub> operator-(const RealVariable& x1, const RealVariable& x2) { return formula(x1)-formula(x2); }
inline BinaryFormula<Mul> operator*(const RealVariable& x1, const RealVariable& x2) { return formula(x1)*formula(x2); }
inline BinaryFormula<Div> operator/(const RealVariable& x1, const RealVariable& x2) { return formula(x1)/formula(x2); }

inline operator+(const Interval& x1, const RealVariable& x2) { return formula(x1)+formula(x2); }
inline operator-(const Interval& x1, const RealVariable& x2) { return formula(x1)-formula(x2); }
inline operator*(const Interval& x1, const RealVariable& x2) { return formula(x1)*formula(x2); }
inline operator/(const Interval& x1, const RealVariable& x2) { return formula(x1)/formula(x2); }

inline operator+(const RealVariable& x1, const Interval& x2) { return formula(x1)+formula(x2); }
inline operator-(const RealVariable& x1, const Interval& x2) { return formula(x1)-formula(x2); }
inline operator*(const RealVariable& x1, const Interval& x2) { return formula(x1)*formula(x2); }
inline operator/(const RealVariable& x1, const Interval& x2) { return formula(x1)/formula(x2); }

inline operator+(const FormulaInterface& x1, const Interval& x2) { return formula(x1)+formula(x2); }
inline operator-(const FormulaInterface& x1, const Interval& x2) { return formula(x1)-formula(x2); }
inline operator*(const FormulaInterface& x1, const Interval& x2) { return formula(x1)*formula(x2); }
inline operator/(const FormulaInterface& x1, const Interval& x2) { return formula(x1)/formula(x2); }

inline operator+(const FormulaInterface& x1, const RealVariable& x2) { return formula(x1)+formula(x2); }
inline operator-(const FormulaInterface& x1, const RealVariable& x2) { return formula(x1)-formula(x2); }
inline operator*(const FormulaInterface& x1, const RealVariable& x2) { return formula(x1)*formula(x2); }
inline operator/(const FormulaInterface& x1, const RealVariable& x2) { return formula(x1)/formula(x2); }

inline operator+(const Interval& x1, const FormulaInterface& x2) { return formula(x1)+formula(x2); }
inline operator-(const Interval& x1, const FormulaInterface& x2) { return formula(x1)-formula(x2); }
inline operator*(const Interval& x1, const FormulaInterface& x2) { return formula(x1)*formula(x2); }
inline operator/(const Interval& x1, const FormulaInterface& x2) { return formula(x1)/formula(x2); }

inline operator+(const RealVariable& x1, const FormulaInterface& x2) { return formula(x1)+formula(x2); }
inline operator-(const RealVariable& x1, const FormulaInterface& x2) { return formula(x1)-formula(x2); }
inline operator*(const RealVariable& x1, const FormulaInterface& x2) { return formula(x1)*formula(x2); }
inline operator/(const RealVariable& x1, const FormulaInterface& x2) { return formula(x1)/formula(x2); }

inline UnaryFormula<Exp> exp(const RealVariable& x) {
    return UnaryFormula<Exp>(Exp(),VariableFormula(x)); }
inline UnaryFormula<Exp> exp(const FormulaInterface& x) {
    return UnaryFormula<Exp>(Exp(),x); }


Constant<int> make_constant(const int& x) { return Constant<int>(x); }
Constant<double> make_constant(const double& x) { return Constant<double>(x); }
Constant<Interval> make_constant(const Interval& x) { return Constant<Interval>(x); }

template<class Op> UnaryFormula<Op> make_unary_formula(const Op& op, const FormulaInterface& a) { return UnaryFormula<Op>(a.clone()); }




template<class Op> UnaryFormula<Op> make_unary_formula(const Op& op, const Formula& a) { return UnaryFormula<Op>(a); }
template<class Op> BinaryFormula<Op> make_binary_formula(const Op& op, const FormulaInterface& a1, const FormulaInterface& a2) { return BinaryFormula<Op>(a1,a2); }

template<class Op> Formula make_binary_formula(const Op& op, const Formula& a1, const Formula& a2) { return new BinaryFormula<Op>(a1,a2); }

inline Formula operator+(const Formula& arg1, const Formula& arg2) { return make_binary_formula(Add(),arg1,arg2); }
inline BinaryFormula<Sub> operator-(const Formula& arg1, const Formula& arg2) { return make_binary_formula(Sub(),arg1,arg2); }
inline BinaryFormula<Mul> operator*(const Formula& arg1, const Formula& arg2) { return make_binary_formula(Mul(),arg1,arg2); }
inline BinaryFormula<Div> operator/(const Formula& arg1, const Formula& arg2) { return make_binary_formula(Div(),arg1,arg2); }

inline BinaryFormula<Add> operator+(const Formula& arg1, const double& arg2) { return make_binary_formula(Add(),arg1,make_constant(arg2)); }
inline BinaryFormula<Sub> operator-(const Formula& arg1, const double& arg2) { return make_binary_formula(Sub(),arg1,make_constant(arg2)); }
inline BinaryFormula<Mul> operator*(const Formula& arg1, const double& arg2) { return make_binary_formula(Mul(),arg1,make_constant(arg2)); }
inline BinaryFormula<Div> operator/(const Formula& arg1, const double& arg2) { return make_binary_formula(Div(),arg1,make_constant(arg2)); }

inline UnaryFormula<Sqr> sqr(const Formula& arg) { return make_unary_formula(Sqr(),arg); }
inline UnaryFormula<Pow> pow(const Formula& arg, int n) { return make_unary_formula(Pow(n),arg); }
inline UnaryFormula<Sqrt> sqrt(const Formula& arg) { return make_unary_formula(Sqrt(),arg); }
inline UnaryFormula<Exp> exp(const Formula& arg) { return make_unary_formula(Exp(),arg); }
inline UnaryFormula<Log> log(const Formula& arg) { return make_unary_formula(Log(),arg); }
inline UnaryFormula<Sin> sin(const Formula& arg) { return make_unary_formula(Sin(),arg); }
inline UnaryFormula<Cos> cos(const Formula& arg) { return make_unary_formula(Cos(),arg); }
inline UnaryFormula<Tan> tan(const Formula& arg) { return make_unary_formula(Tan(),arg); }

*/

inline std::ostream& operator<<(std::ostream& os, const Add& v) { return os << "+"; }
inline std::ostream& operator<<(std::ostream& os, const Sub& v) { return os << "-"; }
inline std::ostream& operator<<(std::ostream& os, const Mul& v) { return os << "*"; }
inline std::ostream& operator<<(std::ostream& os, const Div& v) { return os << "/"; }

inline std::ostream& operator<<(std::ostream& os, const Pow& op) { return os << "pow"; }
inline std::ostream& operator<<(std::ostream& os, const Sqr& op) { return os << "sqr"; }
inline std::ostream& operator<<(std::ostream& os, const Sqrt& op) { return os << "sqrt"; }

inline std::ostream& operator<<(std::ostream& os, const Exp& op) { return os << "exp"; }
inline std::ostream& operator<<(std::ostream& os, const Log& op) { return os << "log"; }
inline std::ostream& operator<<(std::ostream& os, const Sin& op) { return os << "sin"; }
inline std::ostream& operator<<(std::ostream& os, const Cos& op) { return os << "cos"; }
inline std::ostream& operator<<(std::ostream& os, const Tan& op) { return os << "tan"; }


typedef Formula RealFormula;

struct Assignment {
    RealVariable lhs;
    RealFormula rhs;
};

struct Dynamic {
    RealVariable lhs;
    RealFormula rhs;
};

struct Update {
    RealVariable lhs;
    RealFormula rhs;
};




inline Assignment RealVariable::operator=(const RealFormula& f) const {
    Assignment a={*this,f}; return a;
}

inline Dynamic DottedVariable::operator=(const RealFormula& f) const {
    Dynamic d={this->base(),f}; return d;
}

inline Update NextVariable::operator=(const RealFormula& f) const {
    Update d={this->base(),f}; return d;
}

inline Update NextVariable::operator=(const RealVariable& v) const {
    Update d={this->base(),RealFormula(v)}; return d;
}


struct RealPredicate {
    Formula lhs;
    std::string cmp;
    Formula rhs;
};

inline
std::ostream& operator<<(std::ostream& os, const RealPredicate& pred) {
    return os << pred.lhs << pred.cmp << pred.rhs;
}

inline
RealPredicate operator<(const RealFormula& lhs, const RealFormula& rhs) {
    RealPredicate p={lhs,"<",rhs}; return p;
}

inline
RealPredicate operator<=(const RealFormula& lhs, const RealFormula& rhs) {
    RealPredicate p={lhs,"<=",rhs}; return p;
}

inline
RealPredicate operator>(const RealFormula& lhs, const RealFormula& rhs) {
    RealPredicate p={lhs,">",rhs}; return p;
}

inline
RealPredicate operator>=(const RealFormula& lhs, const RealFormula& rhs) {
    RealPredicate p={lhs,">=",rhs}; return p;
}

/*
class Reset
{
  public:
    std::vector<Polynomial<Interval> > _v;
    Reset(int rs, int as) : _v(rs,Polynomial<Interval>(as)) { assert(rs>0); assert(as>=0); }
    int result_size() const { return _v.size(); }
    int argument_size() const { return _v[0].argument_size(); }
    template<class F> Reset& operator,(Assignment<NextVariable,F> a) {
        assert(a.lhs.var.number()<result_size());
        _v[a.lhs.var.number()]=eval<Polynomial<Interval> >(argument_size(),a.rhs); return *this; }
};


class Dynamic
{
  public:
    std::vector<Polynomial<Interval> > _v;
    Dynamic(int s) : _v(s,Polynomial<Interval>(s)) { assert(s>0);}
    int size() const { return _v.size(); }
    int result_size() const { return _v.size(); }
    int argument_size() const { return _v[0].argument_size(); }
    template<class F> Dynamic& operator,(Assignment<DottedVariable,F> a) {
        assert(a.lhs.var.number()<result_size());
        _v[a.lhs.var.number()]=eval<Polynomial<Interval> >(argument_size(),a.rhs); return *this; }
};


class Guard
{
  public:
    Polynomial<Interval> _v;
    Guard(int as) : _v(Polynomial<Interval>(as)) { assert(as>0);}
    int argument_size() const { return _v.argument_size(); }
    template<class F1,class F2> Guard& operator,(Comparison<F1,Gtr,F2> c) {
        _v=eval<Polynomial<Interval> >(argument_size(),c.lhs-c.rhs); return *this; }
    template<class F1,class F2> Guard& operator,(Comparison<F1,Less,F2> c) {
        _v=eval<Polynomial<Interval> >(argument_size(),c.rhs-c.lhs); return *this; }
};


std::ostream& operator<<(std::ostream& os, const Reset& r) { os<<"Reset("<<r.result_size()<<","<<r.argument_size()<<")";
    for(unsigned int i=0; i!=r._v.size(); ++i) { os<<(i==0?"[":",")<<r._v[i]; } return os<<"]"; }
std::ostream& operator<<(std::ostream& os, const Dynamic& d) { os<<"Dynamic("<<d.size()<<")";
    for(unsigned int i=0; i!=d._v.size(); ++i) { os<<(i==0?"[":",")<<d._v[i]; } return os<<"]"; }
std::ostream& operator<<(std::ostream& os, const Guard& g) { os<<"Guard("<<g.argument_size()<<")";
    os<<"["<<g._v; return os<<"]"; }
*/

} // namespace Ariadne

#endif // ARIADNE_FORMULA_H
