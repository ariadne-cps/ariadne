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

class Variable;
class Formula;
class Assignment;
class Comparison;

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



/*! \brief A named variable, suitable for use in formulae and in defining state spaces.
 *  Equality is performed using references by name, so different variables may have the same "name".
 *  \details \sa Space, Formula.
 */
struct Variable
{
  public:
    //! \brief Construct a new named variable with name \a name. If another variable has the same \c name,
    //! the two are \e not considered to be the same quantity.
    Variable(std::string name) : _name_ptr(new std::string(name)) { }
    //! \brief Copy constructer. Returns a C++ variable representing the same named variable.
    Variable(const Variable& v) : _name_ptr(v._name_ptr) { }
    //! \brief The name of the variable.
    const std::string& name() const { return *this->_name_ptr; }
    //! \brief A dynamically-allocated ProjectionExpression representing the variable.
    ExpressionInterface* expression(const Space& spc) const;
    //! \brief Equality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator==(const Variable& v) const { return this->_name_ptr==v._name_ptr; }
    //! \brief Inequality operator. Note that variables with the same name might not be considered equal (the same).
    bool operator!=(const Variable& v) const { return !(*this==v); }

    bool operator<(const Variable& v) const { return *this->_name_ptr<*v._name_ptr; }

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Variable& v);
  private:
    shared_ptr<const std::string> _name_ptr;
};

inline std::ostream& operator<<(std::ostream& os, const Variable& v) { return os << v.name(); }



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
    VariableFormula(const Variable& var) : _var(var) { }
    virtual VariableFormula* clone() const { return new VariableFormula(*this); }
    virtual ExpressionInterface* expression(const Space& spc) const;
    virtual std::ostream& write(std::ostream& os) const { return os<<_var; }
  private:
    Variable _var;
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
 *  \sa Variable, ExpressionInterface
 */
class Formula {
  public:
    Formula(const int& c) : ptr(new ConstantFormula<int>(c)) { }
    Formula(const double& c) : ptr(new ConstantFormula<double>(c)) { }
    Formula(const Interval& c) : ptr(new ConstantFormula<Interval>(c)) { }
    //! \brief The formula "v" in named variable \a v.
    Formula(const Variable& v) : ptr(new VariableFormula(v)) { }
    Formula(FormulaInterface* p) : ptr(p) { }
    Formula(shared_ptr<FormulaInterface> p) : ptr(p) { }
    //! \brief Create an expression \f$\R^n\rightarrow\R\f$ by assigning the variables of the formula the numbers given by \a spc.
    ExpressionInterface* expression(const Space& spc) { return ptr->expression(spc); }

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Formula& e);
  public:
    FormulaPointer ptr;
};

std::ostream& operator<<(std::ostream& os, const Formula& e) { return e.ptr->write(os); }


class AffineFormula
//    : public FormulaInterface
{
  public:
    explicit AffineFormula() : _b(0), _a() { }
    explicit AffineFormula(const int& c) : _b(c), _a() { }
    explicit AffineFormula(const double& c) : _b(c), _a() { }
    explicit AffineFormula(const Interval& c) : _b(c), _a() { }
    //! \brief The formula "v" in named variable \a v.
    explicit AffineFormula(const Variable& v) : _b(0.0), _a() { _a[v]=1.0; }

    AffineFormula* clone() const { return new AffineFormula(*this); }

    //! \brief Create an expression \f$\R^n\rightarrow\R\f$ by assigning the variables of the formula the numbers given by \a spc.
    AffineExpression* expression(const Space& spc) const {
        Vector<Interval> a(spc.dimension());
        for(std::map<Variable,Interval>::const_iterator iter=_a.begin(); iter!=_a.end(); ++iter) {
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
    std::map<Variable,Interval> _a;
};

inline AffineFormula& operator*=(AffineFormula& x, const Interval& c) {
    x._b*=c;
    for(std::map<Variable,Interval>::iterator iter=x._a.begin(); iter!=x._a.end(); ++iter) {
        iter->second*=c;
    }
    return x;
}

inline AffineFormula& operator+=(AffineFormula& x, const AffineFormula& y) {
    x._b+=y._b;
    for(std::map<Variable,Interval>::const_iterator iter=y._a.begin(); iter!=y._a.end(); ++iter) {
        x._a[iter->first]+=iter->second;
    }
    return x;
}

inline AffineFormula& operator-=(AffineFormula& x, const AffineFormula& y) {
    x._b-=y._b;
    for(std::map<Variable,Interval>::const_iterator iter=y._a.begin(); iter!=y._a.end(); ++iter) {
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
    for(std::map<Variable,Interval>::const_iterator iter=x._a.begin(); iter!=x._a.end(); ++iter) {
        if(iter->second!=0) {
            if(iter->second>0) { os << "+"; }
            os << iter->second << "*" << iter->first;
        }
    }
    return os;
}


ExpressionInterface* Variable::expression(const Space& spc) const {
    return new ProjectionExpression(spc.dimension(),spc.index(*this));
}

template<class C>
ExpressionInterface* ConstantFormula<C>::expression(const Space& spc) const {
    return new ConstantExpression(spc.dimension(),Interval(this->_value));
}

ExpressionInterface* VariableFormula::expression(const Space& spc) const {
    return new ProjectionExpression(spc.dimension(),spc.index(this->_var));
}

template<class Op>
ExpressionInterface* UnaryFormula<Op>::expression(const Space& spc) const {
    return new UnaryExpression<Op>(op,arg_ptr->expression(spc)); }

template<class Op>
ExpressionInterface* BinaryFormula<Op>::expression(const Space& spc) const {
    return new BinaryExpression<Op>(op,arg1_ptr->expression(spc),arg2_ptr->expression(spc)); }


FormulaPointer make_formula_pointer(const Float& c) { return FormulaPointer(new ConstantFormula<Float>(c)); }
FormulaPointer make_formula_pointer(const Interval& c) { return FormulaPointer(new ConstantFormula<Interval>(c)); }

Formula operator+(Formula e1, Formula e2) {
    return Formula(new BinaryFormula<Add>(Add(),e1.ptr,e2.ptr)); }
Formula operator-(Formula e1, Formula e2) {
    return Formula(new BinaryFormula<Sub>(Sub(),e1.ptr,e2.ptr)); }
Formula operator*(Formula e1, Formula e2) {
    return Formula(new BinaryFormula<Mul>(Mul(),e1.ptr,e2.ptr)); }
Formula operator/(Formula e1, Formula e2) {
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


Formula exp(Formula e) {
    return Formula(new UnaryFormula<Exp>(Exp(),e.ptr)); }
Formula log(Formula e) {
    return Formula(new UnaryFormula<Log>(Log(),e.ptr)); }
Formula sin(Formula e) {
    return Formula(new UnaryFormula<Sin>(Sin(),e.ptr)); }
Formula cos(Formula e) {
    return Formula(new UnaryFormula<Cos>(Cos(),e.ptr)); }
Formula tan(Formula e) {
    return Formula(new UnaryFormula<Tan>(Tan(),e.ptr)); }



/*
class Formula
{
  public:
    Formula(const FormulaInterface* f_ptr) : _ptr(const_cast<FormulaInterface*>(f_ptr)) { }
    Formula(const Variable& v) : _ptr(new VariableFormula(v)) { }
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
inline VariableFormula<Interval> formula(const Variable& x) { return VariableFormula(c); }
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

inline BinaryFormula<Add> operator+(const Variable& x1, const Variable& x2) { return formula(x1)+formula(x2); }
inline BinaryFormula<Sub> operator-(const Variable& x1, const Variable& x2) { return formula(x1)-formula(x2); }
inline BinaryFormula<Mul> operator*(const Variable& x1, const Variable& x2) { return formula(x1)*formula(x2); }
inline BinaryFormula<Div> operator/(const Variable& x1, const Variable& x2) { return formula(x1)/formula(x2); }

inline operator+(const Interval& x1, const Variable& x2) { return formula(x1)+formula(x2); }
inline operator-(const Interval& x1, const Variable& x2) { return formula(x1)-formula(x2); }
inline operator*(const Interval& x1, const Variable& x2) { return formula(x1)*formula(x2); }
inline operator/(const Interval& x1, const Variable& x2) { return formula(x1)/formula(x2); }

inline operator+(const Variable& x1, const Interval& x2) { return formula(x1)+formula(x2); }
inline operator-(const Variable& x1, const Interval& x2) { return formula(x1)-formula(x2); }
inline operator*(const Variable& x1, const Interval& x2) { return formula(x1)*formula(x2); }
inline operator/(const Variable& x1, const Interval& x2) { return formula(x1)/formula(x2); }

inline operator+(const FormulaInterface& x1, const Interval& x2) { return formula(x1)+formula(x2); }
inline operator-(const FormulaInterface& x1, const Interval& x2) { return formula(x1)-formula(x2); }
inline operator*(const FormulaInterface& x1, const Interval& x2) { return formula(x1)*formula(x2); }
inline operator/(const FormulaInterface& x1, const Interval& x2) { return formula(x1)/formula(x2); }

inline operator+(const FormulaInterface& x1, const Variable& x2) { return formula(x1)+formula(x2); }
inline operator-(const FormulaInterface& x1, const Variable& x2) { return formula(x1)-formula(x2); }
inline operator*(const FormulaInterface& x1, const Variable& x2) { return formula(x1)*formula(x2); }
inline operator/(const FormulaInterface& x1, const Variable& x2) { return formula(x1)/formula(x2); }

inline operator+(const Interval& x1, const FormulaInterface& x2) { return formula(x1)+formula(x2); }
inline operator-(const Interval& x1, const FormulaInterface& x2) { return formula(x1)-formula(x2); }
inline operator*(const Interval& x1, const FormulaInterface& x2) { return formula(x1)*formula(x2); }
inline operator/(const Interval& x1, const FormulaInterface& x2) { return formula(x1)/formula(x2); }

inline operator+(const Variable& x1, const FormulaInterface& x2) { return formula(x1)+formula(x2); }
inline operator-(const Variable& x1, const FormulaInterface& x2) { return formula(x1)-formula(x2); }
inline operator*(const Variable& x1, const FormulaInterface& x2) { return formula(x1)*formula(x2); }
inline operator/(const Variable& x1, const FormulaInterface& x2) { return formula(x1)/formula(x2); }

inline UnaryFormula<Exp> exp(const Variable& x) {
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
