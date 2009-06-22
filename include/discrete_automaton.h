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







class DiscreteType;
class DiscreteValue;
class DiscreteVariable;
class DiscreteSpace;
class DiscreteState;


template<class T>
inline std::vector<T> operator,(const T& t1, const T& t2) {
    std::vector<T> v; v.push_back(t1); v.push_back(t2); return v; }
template<class T>
inline std::vector<T> operator,(const std::vector<T>& v, const T& t) {
    std::vector<T> r(v); r.push_back(t); return r; }


//! \brief A value which can be taken by a discrete variable.
//! \details \sa DiscreteType, DiscreteValue, DiscreteVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscreteValue {
  public:
    DiscreteValue() { }
    DiscreteValue(const std::string& name);

    const std::string& name() const { return this->_name; }

    //! \brief .
    bool operator==(const DiscreteValue& other) const;
    //! \brief .
    bool operator!=(const DiscreteValue& other) const { return !(*this==other); }
    //! \brief .
    bool operator<(const DiscreteValue& other) const { return this->_name < other._name; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteValue& self) {
        return os<<self.name(); }
  private:
    std::string _name;
};

/*! \brief An internal enumerated type, suitable for use when a finite type is needed.
 *  \sa DiscreteType, DiscreteValue, DiscreteVariable, DiscreteSpace, DiscreteState, DiscretePredicate
 */
struct DiscreteType {
  public:
    DiscreteType(const std::string&, uint n, const std::string*);
    DiscreteType(const std::string& nm, const std::vector<DiscreteValue>& ary) : _name(nm), _elements() {
        for(uint i=0; i!=ary.size(); ++i) { _elements.push_back(ary[i].name()); } }
    bool has_value(const std::string& str) {
        for(uint i=0; i!=_elements.size(); ++i) { if(_elements[i]==str) { return true; } } return false; }
    DiscreteValue value(uint i) {
        return DiscreteValue(_elements.at(i)); }
    int id(const std::string&) const;
    bool operator==(const DiscreteType& other) const;
    bool operator!=(const DiscreteType& other) const;
    //! \brief .
    const std::string& name() const { return this->_name; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteType& self) {
        os << self.name() << "=";
        for(uint i=0; i!=self._elements.size(); ++i) { os<<(i==0?'{':',')<<self._elements[i]; } return os << '}'; }
  private:
    std::string _name;
    std::vector<std::string> _elements;

};


inline DiscreteValue::DiscreteValue(const std::string& name)
    : _name(name) { }

inline bool
DiscreteValue::operator==(const DiscreteValue& other) const {
   return this->_name == other._name;
}

/*
inline DiscreteValue::DiscreteValue(const DiscreteType& tp, const std::string& name)
    : _type_ptr(&tp), _value(tp.id(name)) { }

inline bool
DiscreteValue::operator==(const DiscreteValue& other) const {
    ARIADNE_ASSERT_MSG(this->type()==other.type(),"Type mismatch: type of "<<this->name()<<" is "<<this->type().name()<<"; type of "<<other.name()<<" is "<<other.type().name()<<".");
    return this->_value == other._value;
}
*/


//! \brief A named variable which can take values of a given discrete type.
//! \details \sa DiscreteType, DiscreteValue, DiscreteVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscreteVariable {
    DiscreteVariable(std::string name, DiscreteType type) : _name(name), _type(type) { }
  public:
    //! \brief .
    const DiscreteType& type() const { return this->_type; }
    //! \brief .
    const std::string& name() const { return this->_name; }
    //! \brief .
    bool operator==(const DiscreteVariable& other) const { return this->_name==other._name; }
    bool operator!=(const DiscreteVariable& other) const { return this->_name!=other._name; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteVariable& self) {
        return os << self.name() << ":" << self.type().name(); }
  private:
    std::string _name;
    DiscreteType _type;
};

//! \brief A finite state space, consisting of valuations for a set of discrete variables..
//! \details \sa DiscreteType, DiscreteValue, DiscreteVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscreteSpace {
    friend class DiscreteState;
  public:
    DiscreteSpace(const std::vector<DiscreteVariable>& variables) : _variables(variables) { }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteSpace& self) {
        for(uint i=0; i!=self._variables.size(); ++i) { os << (i==0?'(':',') << self._variables[i].name(); } return os << ')'; };
  private:
    std::vector<DiscreteVariable> _variables;
};



//! \brief A predicate over some discrete variables.
//! \details \sa DiscreteType, DiscreteValue, DiscreteVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscreteState {
  public:
    DiscreteState(const DiscreteSpace& spc) : _variables(spc._variables), _values(_variables.size()) { }

    DiscreteSpace space() const { return reinterpret_cast<DiscreteSpace const&>(this->_variables); }
    DiscreteValue operator[](const DiscreteVariable& var) const {
        for(uint i=0; i!=this->_variables.size(); ++i) {
            if(_variables[i]==var) { return this->_values[i]; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<var<<" is not a member of space "<<this->space()); }
    DiscreteValue& operator[](const DiscreteVariable& var) {
        for(uint i=0; i!=this->_variables.size(); ++i) {
            if(_variables[i]==var) { return this->_values[i]; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<var<<" is not a member of space "<<this->space()); }
    bool operator==(const DiscreteState& other) const {
        for(uint i=0; i!=_values.size(); ++i) { if(this->_values[i]!=other._values[i]) { return false; } } return true; }
    bool operator!=(const DiscreteState& other) const {
        return  !(*this==other); }
    bool operator<(const DiscreteState& other) const {
        for(uint i=0; i!=_values.size(); ++i) {
            if(this->_values[i]!=other._values[i]) { return this->_values[i]<other._values[i]; } }
        return false; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const DiscreteState& self) {
        for(uint i=0; i!=self._values.size(); ++i) { os << (i==0?'(':',') << self._values[i]; } return os << ')'; }
  private:
    std::vector<DiscreteVariable> _variables;
    std::vector<DiscreteValue> _values;
};

//! \brief A predicate over some discrete variables.
//! \details \sa DiscreteValue, DiscreteVariable, DiscreteSpace, DiscreteState, DiscretePredicate
struct DiscretePredicate {
  public:
    //! \brief Evaluate the predicate at the given state.
    bool operator() (const DiscreteState&) const;

    //! \brief .
    friend DiscretePredicate operator==(const DiscreteVariable&, const DiscreteValue&);
    //! \brief .
    friend DiscretePredicate operator!=(const DiscreteVariable&, const DiscreteValue&);
    //! \brief .
    friend DiscretePredicate operator||(const DiscretePredicate&, const DiscretePredicate&);
    //! \brief .
    friend DiscretePredicate operator&&(const DiscretePredicate&, const DiscretePredicate&);
    //! \brief .
    friend DiscretePredicate operator!(const DiscretePredicate&);

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const DiscretePredicate&);
  private:
    std::vector< std::vector< std::pair<DiscreteVariable, DiscreteValue> > > _cnf;
};

inline
bool DiscretePredicate::operator() (const DiscreteState& state) const {
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
DiscretePredicate operator==(const DiscreteVariable& var, const DiscreteValue& val) {
    std::vector< std::pair<DiscreteVariable, DiscreteValue> > clause(1u,std::make_pair(var,val));
    DiscretePredicate result; result._cnf.push_back(clause); return result;
}

inline
DiscretePredicate operator!=(const DiscreteVariable&, const DiscreteValue&);

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
    for(uint i=0; i!=self._cnf.size(); ++i) {
        if(i!=0) { os << " & "; }
        for(uint j=0; j!=self._cnf[i].size(); ++j) {
            os << (j==0?"(":" | ") << self._cnf[i][j].first.name()<<"=="<<self._cnf[i][j].second.name();
        }
        os << ")";
    }
    return os;
}


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
  private:
    uint _id;
  private:
    static std::vector<std::string> _names;
};


} // namespace Ariadne

#endif // ARIADNE_DISCRETE_AUTOMATON_H
