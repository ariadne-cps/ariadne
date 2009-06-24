/***************************************************************************
 *            discrete_automaton.h
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


/*! \file discrete_automaton.h
 *  \brief Discrete formulae over variables
 */

#ifndef ARIADNE_DISCRETE_AUTOMATON_H
#define ARIADNE_DISCRETE_AUTOMATON_H

#include <cstdarg>
#include <iostream>
#include <string>
#include <vector>

#include "macros.h"
#include "pointer.h"


namespace Ariadne {





typedef std::string String;

class EnumeratedType;
class EnumeratedVariable;
class DiscreteValue;
class DiscreteAssignment;
class DiscreteUpdate;
class DiscretePredicate;
class DiscreteSpace;
class DiscreteState;

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

NextEnumeratedVariable next(const EnumeratedVariable& v) {
    return NextEnumeratedVariable(v);
}


typedef EnumeratedType DiscreteType;
typedef EnumeratedVariable DiscreteVariable;
typedef NextEnumeratedVariable NextDiscreteVariable;

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


//! \brief A predicate over some discrete variables.
//! \details \sa DiscreteValue, EnumeratedVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscretePredicate {
    friend class DiscreteConstraint;
  public:
    DiscretePredicate() { }
    DiscretePredicate(bool val) { }
    //! \brief Evaluate the predicate at the given state.
    bool operator() (const DiscreteValuation&) const;

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


struct DiscreteAssignment {
    DiscreteVariable lhs;
    DiscreteFormula rhs;
};

struct DiscreteUpdate {
    DiscreteVariable lhs;
    DiscreteFormula rhs;
};

DiscreteAssignment EnumeratedVariable::operator=(const DiscreteValue& val) const {
    DiscreteAssignment a={*this,DiscreteFormula(val)}; return a; }

DiscreteUpdate NextEnumeratedVariable::operator=(const DiscreteVariable& var) const {
    DiscreteUpdate a={this->base(),DiscreteFormula(var)}; return a; }

DiscreteUpdate NextEnumeratedVariable::operator=(const DiscreteValue& val) const {
    DiscreteUpdate a={this->base(),DiscreteFormula(val)}; return a; }

DiscreteUpdate NextEnumeratedVariable::operator=(const std::string& str) const {
    DiscreteUpdate a={this->base(),DiscreteFormula(DiscreteValue(str))}; return a; }


//! \brief A discrete event.
struct DiscreteEvent {
  public:
    DiscreteEvent(const std::string& name) { _names.push_back(name); this->_id=_names.size()-1; }
    DiscreteEvent(int n) { std::stringstream ss; ss<<"e"<<n; _names.push_back(ss.str()); this->_id=_names.size()-1; ; }
    int id() const { return this->_id; }
    const std::string& name() const { return _names[this->_id]; }
    bool operator==(const DiscreteEvent& other) const { return this->id()==other.id(); }
    bool operator!=(const DiscreteEvent& other) const { return this->id()!=other.id(); }
    bool operator<(const DiscreteEvent& other) const { return this->id()<other.id(); }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteEvent& self) { return os << self.name(); }
  private:
    uint _id;
  private:
    static std::vector<std::string> _names;
};

//! \brief A finite or cofinite set of discrete events.
struct DiscreteEventSet {
  public:
    DiscreteEventSet() : _events(), _is_complement(false) { }
    DiscreteEventSet(const std::set<DiscreteEvent>& s) : _events(s), _is_complement(false) { }
    DiscreteEventSet(const std::vector<DiscreteEvent>& v) : _events(v.begin(),v.end()), _is_complement(false) { }
    DiscreteEventSet(const DiscreteEvent& e) : _events(), _is_complement(false) { this-> _events.insert(e); }

    static DiscreteEventSet none() { DiscreteEventSet r; return r; }
    static DiscreteEventSet all() { DiscreteEventSet r; r._is_complement=true; return r; }

    friend DiscreteEventSet operator!(const DiscreteEvent& e);
    friend DiscreteEventSet operator!(const DiscreteEventSet& s);

    friend std::ostream& operator<<(std::ostream& os, const DiscreteEventSet& self) {
        if(self._is_complement) { os<<"!"; } return os << self._events; }
  private:
    std::set<DiscreteEvent> _events;
    bool _is_complement;
};

DiscreteEventSet operator!(const DiscreteEvent& e) {
    return !(DiscreteEventSet(e)); }
DiscreteEventSet operator!(const DiscreteEventSet& s) {
    DiscreteEventSet r; r._events=s._events; r._is_complement=!s._is_complement; return r; }
DiscreteEventSet operator!(const std::set<DiscreteEvent>& s) {
    return !DiscreteEventSet(s); }



} // namespace Ariadne

#endif // ARIADNE_DISCRETE_AUTOMATON_H
