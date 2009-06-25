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

class Variable;
class Space;

class EnumeratedType;
class EnumeratedVariable;
class EnumeratedValue;
class EnumeratedFormula;

class RealVariable;
class RealFormula;

// Functionality for building arrays or std::vector from comma-separated lists
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
inline std::vector<std::string> operator,(ArrayBuilderTag, const char* s) {
    std::vector<std::string> v; v.push_back(std::string(s)); return v; }




/*! \brief An ASCII string. */
typedef std::string String;

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
    friend EventSet operator!(const std::set<Event>& s);

    friend std::ostream& operator<<(std::ostream& os, const EventSet& self);
  private:
    std::set<Event> _events;
    bool _is_complement;
};

inline EventSet operator!(const Event& e) { return !(EventSet(e)); }
inline EventSet operator!(const EventSet& s) { EventSet r; r._events=s._events; r._is_complement=!s._is_complement; return r; };
inline EventSet operator!(const std::set<Event>& s) { return !EventSet(s); }

inline std::ostream& operator<<(std::ostream& os, const EventSet& self) {
    if(self._is_complement) { os<<"!"; } return os << self._events; }

template<class LHS, class RHS>
struct Assignment
{
    Assignment(const LHS& l, const RHS& r) : lhs(l), rhs(r) { }
    LHS lhs; RHS rhs;
};


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
    //! \brief Equality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator==(const Variable& v) const { return this->_name_ptr==v._name_ptr; }
    //! \brief Inequality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator!=(const Variable& v) const { return !(*this==v); }
    //! \brief Less-than comparison orders %Variable objects lexicographically by name.
    bool operator<(const Variable& v) const { return *this->_name_ptr<*v._name_ptr; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Variable& v) {
        return os << v.name(); }
  private:
    shared_ptr<std::string> _name_ptr;
};


//! \brief A named discrete variable.
struct DiscreteVariable : public Variable {
    DiscreteVariable(const std::string& name) : Variable(name) { }
};


//! \brief A named variable which can take values of a given enumerated type.
//! \details \sa EnumeratedType, EnumeratedVariable, DiscretePredicate
struct EnumeratedVariable : public DiscreteVariable {
  public:
    //! \brief Create an internal enumerated variable with the given \a name and \a type.
    EnumeratedVariable(const String& name, const EnumeratedType& type) : DiscreteVariable(name), _type_ptr(&type) { }
    //! \brief The internal type of the variable.
    const EnumeratedType& type() const { return *this->_type_ptr; }
    //! \brief Makes a discrete assignment object.
    Assignment<EnumeratedVariable,EnumeratedFormula> operator=(const EnumeratedValue& val) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const EnumeratedVariable& v) {
        return os << v.name(); }
        //return os << v.name() << ":" << v.type().name(); }
  private:
    const EnumeratedType* _type_ptr;
};


//! \brief A named integer variable, suitable for use in formulae and in defining the discrete state.
//! Equality is performed using references by name, so different variables may have the same "name".
//! \details \sa Space, Formula.
struct IntegerVariable : public DiscreteVariable
{
  public:
    //! \brief Construct a new named integer variable with name \a name. If another variable has the same \c name,
    //! the two are \e not considered to be the same quantity.
    IntegerVariable(std::string name) : DiscreteVariable(name) { }
    //! \brief Copy constructer. Returns a C++ variable representing the same named variable.
    IntegerVariable(const IntegerVariable& v) : DiscreteVariable(v) { }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const IntegerVariable& v) {
        return os << v.name(); }
        //return os << v.name() << ":Integer"; }
  private:
    shared_ptr<const std::string> _name_ptr;
};

//! \brief A named continuous variable.
struct ContinuousVariable : public Variable {
    ContinuousVariable(const std::string& name) : Variable(name) { }
};


//! \brief A named real variable, suitable for use in formulae and in defining state spaces.
//! Equality is performed using references by name, so different variables may have the same "name".
//! \details \sa Space, Formula.
struct RealVariable : public ContinuousVariable
{
  public:
    //! \brief Construct a new named variable with name \a name. If another variable has the same \c name,
    //! the two are \e not considered to be the same quantity.
    RealVariable(std::string name) : ContinuousVariable(name) { }
    //! \brief Copy constructer. Returns a C++ variable representing the same named variable.
    RealVariable(const RealVariable& v) : ContinuousVariable(v) { }
    //! \brief A dynamically-allocated ProjectionExpression representing the variable.
    ExpressionInterface* expression(const Space& spc) const;

    //! \brief Create an assignment object representing the variable's defining equation.
    Assignment<RealVariable,RealFormula> operator=(const RealFormula& f) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const RealVariable& v) {
        return os << v.name(); }
        //return os << v.name() << ":Real"; }
  private:
    shared_ptr<const std::string> _name_ptr;
};


/*! \brief A space defined as a list of named variables.
 *  \details \sa Variable
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
    const Variable& operator[](unsigned int i) const { return _variables.at(i); }
    const Variable& variable(unsigned int i) const { return _variables.at(i); }

    //! \brief The index of the named variable \a v.
    unsigned int index(const Variable& v) const {
        for(uint i=0; i!=_variables.size(); ++i) {
            if(v==_variables[i]) { return i; } }
        ARIADNE_ASSERT_MSG(false,"Variable "<<v<<" is not in the Space "<<*this);
        return _variables.size(); }
    //! \brief Append the named variable \a v to the variables defining the space.
    Space& operator,(const Variable& v) {
        for(uint i=0; i!=_variables.size(); ++i) {
            ARIADNE_ASSERT_MSG(_variables[i]!=v,"Variable "<<v<<" is already a variable of the Space "<<*this);
        }
        _variables.push_back(v); return *this; }
  public:
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const Space&);
  private:
    std::vector<Variable> _variables;
};

inline std::ostream& operator<<(std::ostream& os, const Space& spc) { return os << spc._variables; }

inline Space operator,(const Variable& v1, const Variable& v2) {
    Space r; r,v1,v2; return r; }



//! \brief A predicate over some discrete variables.
//! \details \sa EnumeratedValue, EnumeratedVariable, DiscretePredicate
struct EnumeratedValue {
  public:
    //! \brief Construct a discrete value based on the string \a str.
    EnumeratedValue(const String& str) : _value(str) { }
    //! \brief The string representation of the value.
    operator const String& () const { return this->_value; }
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const EnumeratedValue& self) {
        return os << self._value; }
  private:
    std::string _value;
};


//! \brief A formula returning an EnumeratedVariable.
struct EnumeratedFormula {
    EnumeratedFormula(const EnumeratedValue& val) : _str(val) { }
    EnumeratedFormula(const EnumeratedVariable& var) : _str(var.name()) { }
    friend std::ostream& operator<<(std::ostream& os, const EnumeratedFormula& self) {
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



struct GtrZero {}; struct LessZero {}; struct Gtr { }; struct Less { }; struct Equal { };

struct And { template<class T> T operator()(const T& a1, const T& a2) const { return a1 && a2; } };
struct Or { template<class T> T operator()(const T& a1, const T& a2) const { return a1 || a2; } };
struct Not { template<class T> T operator()(const T& a) const { return !a; } };


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
class RealFormula {
  public:
    RealFormula(const int& c) : ptr(new ConstantFormula<int>(c)) { }
    RealFormula(const double& c) : ptr(new ConstantFormula<double>(c)) { }
    RealFormula(const Interval& c) : ptr(new ConstantFormula<Interval>(c)) { }
    //! \brief The formula "v" in named variable \a v.
    RealFormula(const RealVariable& v) : ptr(new VariableFormula(v)) { }
    RealFormula(FormulaInterface* p) : ptr(p) { }
    RealFormula(shared_ptr<FormulaInterface> p) : ptr(p) { }
    //! \brief Create an expression \f$\R^n\rightarrow\R\f$ by assigning the variables of the formula the numbers given by \a spc.
    ExpressionInterface* expression(const Space& spc) { return ptr->expression(spc); }

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const RealFormula& e);
  public:
    FormulaPointer ptr;
};

inline std::ostream& operator<<(std::ostream& os, const RealFormula& e) { return e.ptr->write(os); }


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

//! \related RealFormula \brief .
inline RealFormula operator+(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Add>(Add(),e1.ptr,e2.ptr)); }
//! \related RealFormula \brief .
inline RealFormula operator-(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Sub>(Sub(),e1.ptr,e2.ptr)); }
//! \related RealFormula \brief .
inline RealFormula operator*(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Mul>(Mul(),e1.ptr,e2.ptr)); }
//! \related RealFormula \brief .
inline RealFormula operator/(RealFormula e1, RealFormula e2) {
    return RealFormula(new BinaryFormula<Div>(Div(),e1.ptr,e2.ptr)); }


//! \related RealFormula \brief .
inline RealFormula exp(RealFormula e) {
    return RealFormula(new UnaryFormula<Exp>(Exp(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula log(RealFormula e) {
    return RealFormula(new UnaryFormula<Log>(Log(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula sin(RealFormula e) {
    return RealFormula(new UnaryFormula<Sin>(Sin(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula cos(RealFormula e) {
    return RealFormula(new UnaryFormula<Cos>(Cos(),e.ptr)); }
//! \related RealFormula \brief .
inline RealFormula tan(RealFormula e) {
    return RealFormula(new UnaryFormula<Tan>(Tan(),e.ptr)); }


inline std::ostream& operator<<(std::ostream& os, const Less& v) { return os << "<="; }
inline std::ostream& operator<<(std::ostream& os, const Gtr& v) { return os << ">="; }
inline std::ostream& operator<<(std::ostream& os, const Equal& v) { return os << "=="; }

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

template<class CVAR>
struct Dotted
{
    CVAR base;
    Assignment<Dotted<RealVariable>,RealFormula> operator=(const RealFormula& f) {
        return Assignment<Dotted<RealVariable>,RealFormula>(*this,f); }
    Assignment<Dotted<RealVariable>,RealFormula> operator=(const RealVariable& v) {
        return Assignment<Dotted<RealVariable>,RealFormula>(*this,RealFormula(v)); }
    friend Dotted<RealVariable> dot(const RealVariable& v);
    friend std::ostream& operator<<(std::ostream& os, const Dotted<CVAR>& self) {
        return os << "next(" << self.base.name() << ")"; }
  private:
    Dotted(const CVAR& v) : base(v) { }
};

inline Dotted<RealVariable> dot(const RealVariable& v) { return Dotted<RealVariable>(v); }

typedef Assignment<EnumeratedVariable,EnumeratedFormula> EnumeratedAssignment;
typedef Assignment<Next<EnumeratedVariable>,EnumeratedFormula> EnumeratedUpdate;
typedef Assignment<RealVariable,RealFormula> RealAssignment;
typedef Assignment<Next<RealVariable>,RealFormula> RealUpdate;
typedef Assignment<Dotted<RealVariable>,RealFormula> RealDynamic;


template<class T> class PredicateInterface {
  public:
    virtual PredicateInterface<T>* clone() const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
};

template<class T> inline
std::ostream& operator<<(std::ostream& os, const PredicateInterface<T>& pred) {
    return pred.write(os);
}


template<class T> struct ConstantPredicate : public PredicateInterface<T> {
  public:
    ConstantPredicate(const T& c) : _value(c) { }
    virtual ConstantPredicate<T>* clone() const { return new ConstantPredicate<T>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os<<_value; }
  private:
    T _value;
};

template<class T, class LHS, class Cmp, class RHS> struct ComparisonPredicate : public PredicateInterface<T> {
  public:
    ComparisonPredicate(const LHS& l, const Cmp& c, const RHS& r) : lhs_ptr(new LHS(l)), cmp(c), rhs_ptr(new RHS(r)) { }
    virtual ComparisonPredicate<T,LHS,Cmp,RHS>* clone() const { return new ComparisonPredicate<T,LHS,Cmp,RHS>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os<<*lhs_ptr<<cmp<<*rhs_ptr; }
  private:
    shared_ptr<LHS> lhs_ptr; Cmp cmp; shared_ptr<RHS> rhs_ptr;
};


template<class T, class Op> struct UnaryPredicate : public PredicateInterface<T> {
    typedef shared_ptr< PredicateInterface<T> > PredicatePointer;
    UnaryPredicate(Op o, PredicatePointer a) : op(o), arg_ptr(a) { }
    UnaryPredicate(Op o, const PredicateInterface<T>& a) : op(o), arg_ptr(a.clone()) { }
    virtual UnaryPredicate<T,Op>* clone() const { return new UnaryFormula<Op>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os<<op<<"("<<*arg_ptr<<")"; }
    //template<class Res> Res evaluate() const { return op(arg.evaluate<Res>()); }
    Op op; PredicatePointer arg_ptr;
};

template<class T, class Op> struct BinaryPredicate : public PredicateInterface<T> {
    typedef shared_ptr< PredicateInterface<T> > PredicatePointer;
    BinaryPredicate(Op o, PredicatePointer a1, PredicatePointer a2) : op(o), arg1_ptr(a1), arg2_ptr(a2) { }
    virtual BinaryPredicate<T,Op>* clone() const { return new BinaryFormula<Op>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os<<*arg1_ptr<<op<<*arg2_ptr; }
    Op op; PredicatePointer arg1_ptr; PredicatePointer arg2_ptr;
};

template<class T> class Predicate;
template<class T> std::ostream& operator<<(std::ostream&, const Predicate<T>&);

//! \brief .
template<class T> class Predicate {
  public:
    //! \brief .
    Predicate<T>(const T& c) : _ptr(new ConstantPredicate<T>(c)) { }
    //! \brief .
    Predicate<T>(PredicateInterface<T>* p) : _ptr(p) { }
    //! \brief .
    Predicate<T>(shared_ptr< PredicateInterface<T> > p) : _ptr(p) { }

    //! \brief Write to an output stream.
    friend std::ostream& operator<< <>(std::ostream& os, const Predicate<T>& self);
  private:
    shared_ptr< PredicateInterface<T> > _ptr;
};

//! \related Predicate<T> \brief .
template<class T> inline Predicate<T> operator&&(Predicate<T> e1, Predicate<T> e2) {
    return Predicate<T>(new BinaryPredicate<T,And>(And(),e1.ptr,e2.ptr)); }
//! \related Predicate<T> \brief .
template<class T> inline Predicate<T> operator||(Predicate<T> e1, Predicate<T> e2) {
    return Predicate<T>(new BinaryPredicate<T,Or>(Or(),e1.ptr,e2.ptr)); }
//! \related Predicate<T> \brief .
template<class T> inline Predicate<T> operator!(Predicate<T> e) {
    return Predicate<T>(new UnaryPredicate<T,Or>(Or(),e.ptr)); }

template<class T> inline std::ostream& operator<<(std::ostream& os, const Predicate<T>& self) {
    return os << *self._ptr; }

//! \related EnumeratedVariable \brief .
inline
Predicate<bool>
operator==(const EnumeratedVariable& lhs, const EnumeratedValue& rhs) {
    return Predicate<bool>(new ComparisonPredicate<bool,EnumeratedVariable,Equal,EnumeratedValue>(lhs,Equal(),rhs));
}

//! \related EnumeratedVariable \brief .
inline
Predicate<bool>
operator==(const EnumeratedVariable& lhs, const String& rhs) {
    return Predicate<bool>(new ComparisonPredicate<bool,EnumeratedVariable,Equal,EnumeratedValue>(lhs,Equal(),EnumeratedValue(rhs)));
}

//! \related RealFormula \brief .
inline
Predicate<tribool>
operator<=(const RealFormula& lhs, const RealFormula& rhs) {
    return Predicate<tribool>(new ComparisonPredicate<tribool,RealFormula,Less,RealFormula>(lhs,Less(),rhs));
}

//! \related RealFormula \brief .
inline
Predicate<tribool>
operator>=(const RealFormula& lhs, const RealFormula& rhs) {
    return Predicate<tribool>(new ComparisonPredicate<tribool,RealFormula,Gtr,RealFormula>(rhs,Gtr(),lhs));
}

typedef Predicate<bool> DiscretePredicate;
typedef Predicate<tribool> ContinuousPredicate;

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
