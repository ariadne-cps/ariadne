/***************************************************************************
 *            symbolic_function.h
 *
 *  Copyright 2008-12  Pieter Collins
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

/*! \file symbolic_function.h
 *  \brief Symbolic functions
 */

#ifndef ARIADNE_SYMBOLIC_FUNCTION_H
#define ARIADNE_SYMBOLIC_FUNCTION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function_interface.h"

#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/container.h"
#include "utility/metaprogramming.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"

#include "function/function_mixin.h"

namespace Ariadne {

//------------------------ Formula functions  -----------------------------------//

//! A function defined by a formula
template<class X>
struct ScalarFormulaFunction
    : ScalarFunctionMixin<ScalarFormulaFunction<X>,InformationTag<X>>
{
    typedef InformationTag<X> P;
    Nat _argument_size;
    Formula<X> _formula;

    ScalarFormulaFunction(Nat as, const Formula<X>& f) : _argument_size(as), _formula(f) { }
    operator Formula<X>() const { return _formula; }

    virtual Nat argument_size() const { return _argument_size; }
    virtual ScalarFunctionInterface<P>* _derivative(uint j) const { return new ScalarFormulaFunction<X>(_argument_size,Ariadne::derivative(_formula,j)); }
    virtual std::ostream& write(std::ostream& os) const { return os << this->_formula; }
    virtual std::ostream& repr(std::ostream& os) const { return os << "FormulaFunction("<<this->_argument_size<<","<<this->_formula<<")"; }
    template<class Y> void _compute(Y& r, const Vector<Y>& x) const { r=Ariadne::evaluate(_formula,x); }
};

inline EffectiveScalarFunction function(Nat n, Formula<Real> f) { return EffectiveScalarFunction(new ScalarFormulaFunction<Real>(n,f)); }

typedef ScalarFormulaFunction<Real> RealScalarFormulaFunction;

//! A vector function defined by formulae
template<class X>
struct VectorFormulaFunction
    : VectorFunctionMixin<VectorFormulaFunction<X>,InformationTag<X>>
{
    Nat _argument_size;
    Vector< Formula<X> > _formulae;

    VectorFormulaFunction(uint as, const List< Formula<X> >& f) : _argument_size(as), _formulae(f) { }
    VectorFormulaFunction(uint as, const Vector< Formula<X> >& f) : _argument_size(as), _formulae(f) { }

    virtual Nat result_size() const { return this->_formulae.size(); }
    virtual Nat argument_size() const { return this->_argument_size; }
    virtual ScalarFormulaFunction<X>* _get(uint i) const { return new ScalarFormulaFunction<X>(this->_argument_size,this->_formulae[i]); }
    virtual std::ostream& write(std::ostream& os) const { return os << this->_formulae; }
    virtual std::ostream& repr(std::ostream& os) const { return os << "VectorFormulaFunction("<<this->result_size()<<","<<this->argument_size()<<","<<this->_formulae<<")"; }
    template<class Y> void _compute(Vector<Y>& r, const Vector<Y>& x) const { r=Ariadne::cached_evaluate(this->_formulae,x); }
};



//------------------------ Arithmetic scalar functions  -----------------------------------//


//! A constant function f(x)=c
template<class X>
struct ConstantFunction
    : ScalarFunctionMixin<ConstantFunction<X>,InformationTag<X>>
{
    typedef InformationTag<X> P;
  public:
    Nat _argument_size;
    X _value;

    ConstantFunction(uint as, const X& c) : _argument_size(as), _value(c) { }
    operator X() const { return _value; }

    virtual Nat argument_size() const { return _argument_size; }
    virtual ScalarFunctionInterface<P>* _derivative(uint j) const { return new ConstantFunction<P>(_argument_size,P(0)); }
    virtual std::ostream& repr(std::ostream& os) const { return os << this->_value; }
    virtual std::ostream& write(std::ostream& os) const { return os << "CF[R"<<this->_argument_size<<"]("<<_value<<")"; }
    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=x.zero_element()+_value; }
};

//! A coordinate function \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=x_i\f$.
template<class P>
struct CoordinateFunction
    : ScalarFunctionMixin<CoordinateFunction<P>,P>
{
    typedef CanonicalNumberType<P> X;

    Nat _argument_size;
    Nat _index;

    CoordinateFunction(uint as, uint i) : _argument_size(as), _index(i) { }
    Nat index() const { return _index; }

    virtual Nat argument_size() const { return _argument_size; }
    virtual ScalarFunctionInterface<P>* _derivative(uint j) const {
        if(j==_index) { return new ConstantFunction<X>(_argument_size,X(1)); }
        else { return new ConstantFunction<X>(_argument_size,X(0)); } }
    virtual std::ostream& repr(std::ostream& os) const { return os << "x"<<this->_index; }
    virtual std::ostream& write(std::ostream& os) const { return os << "IF[R"<<this->_argument_size<<"](x"<<this->_index<<")"; }
    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=x[_index]; }
};


//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
template<class P>
struct UnaryFunction
    : ScalarFunctionMixin< UnaryFunction<P>, P >
{
    typedef CanonicalNumberType<P> X;
  public:
    UnaryFunction(const OperatorCode& op, const ScalarFunction<P>& arg)
        : _op(op), _arg(arg) { }
    virtual UnaryFunction<P>* clone() const { return new UnaryFunction<P>(*this); }
    virtual Nat argument_size() const {
        return this->_arg.argument_size(); }

    virtual ScalarFunctionInterface<P>* _derivative(uint j) const {
        return static_cast<const ScalarFunctionInterface<P>&>(this->derivative(j))._clone();
    }

    virtual ScalarFunction<P> derivative(uint j) const {
        switch(_op) {
            case POS: return _arg.derivative(j);
            case NEG: return -_arg.derivative(j);
            case REC: return -_arg.derivative(j)/sqr(_arg);
            case SQR: return 2*_arg.derivative(j)*_arg;
            case SQRT: return _arg.derivative(j)/(2*sqrt(_arg));
            case EXP: return _arg*_arg.derivative(j);
            case LOG: return _arg.derivative(j)/_arg;
            case SIN: return _arg.derivative(j)*cos(_arg);
            case COS: return -_arg.derivative(j)*sin(_arg);
            default: ARIADNE_FAIL_MSG("Unknown unary function "<<this->_op);
        }
    }

    virtual std::ostream& repr(std::ostream& os) const {
        return os << "UF[R" << this->argument_size() << "](" << *this << ")"; }
    virtual std::ostream& write(std::ostream& os) const {
        return os << _op << '(' << _arg << ')'; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=Ariadne::compute(_op,_arg.evaluate(x)); }

    OperatorCode _op;
    ScalarFunction<P> _arg;
};


template<class P> ScalarFunction<P> sqr(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(SQR,f)); }
template<class P> ScalarFunction<P> sqrt(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(SQRT,f)); }
template<class P> ScalarFunction<P> sin(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(SIN,f)); }
template<class P> ScalarFunction<P> cos(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(COS,f)); }

template<class P>
struct BinaryFunction
    : ScalarFunctionMixin< BinaryFunction<P>, P >
{
    typedef CanonicalNumberType<P> X;
  public:
    BinaryFunction(OperatorCode op, const ScalarFunction<P>& arg1, const ScalarFunction<P>& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) { ARIADNE_ASSERT_MSG(arg1.argument_size()==arg2.argument_size(),"op='"<<op<<"', arg1="<<arg1<<", arg2="<<arg2); }
    virtual BinaryFunction<P>* clone() const { return new BinaryFunction<P>(*this); }
    virtual Nat argument_size() const {
        return this->_arg1.argument_size(); }

    virtual ScalarFunctionInterface<P>* _derivative(uint j) const {
        return static_cast<const ScalarFunctionInterface<P>&>(this->derivative(j))._clone();
    }

    virtual ScalarFunction<P> derivative(uint j) const {
        switch(_op) {
            case ADD:
                return _arg1.derivative(j)+_arg2.derivative(j);
            case SUB:
                return _arg1.derivative(j)-_arg2.derivative(j);
            case MUL:
                return _arg1.derivative(j)*_arg2+_arg1*_arg2.derivative(j);
            case DIV:
                if(dynamic_cast<const ConstantFunction<X>*>(_arg2.raw_pointer())) {
                    return _arg1.derivative(j)/_arg2;
                } else {
                    return _arg1.derivative(j)/_arg2-_arg2.derivative(j)*_arg1/sqr(_arg2);
                }
            default: ARIADNE_FAIL_MSG("Unknown binary function "<<this->_op);
        }
    }

    virtual std::ostream& repr(std::ostream& os) const {
        return os << "BF[R" << this->argument_size() << "](" << *this << ")"; }
    virtual std::ostream& write(std::ostream& os) const {
        if(_op==ADD || _op==SUB) { return os << '(' << _arg1 << symbol(_op) << _arg2 << ')'; }
        else { return os << _arg1 << symbol(_op) << _arg2; } }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=Ariadne::compute(_op,_arg1.evaluate(x),_arg2.evaluate(x)); }

    OperatorCode _op;
    ScalarFunction<P> _arg1;
    ScalarFunction<P> _arg2;
};


// \brief The power function \f$(x,n)\mapsto x^n\f$.
template<class P>
class PowerFunction
    : public ScalarFunctionMixin< PowerFunction<P>, P >
{
    typedef CanonicalNumberType<P> X;
  public:
    PowerFunction(OperatorCode op, const ScalarFunction<P>& arg1, const Int& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {  }
    virtual PowerFunction<P>* clone() const { return new PowerFunction<P>(*this); }
    virtual Nat argument_size() const {
        return this->_arg1.argument_size(); }

    virtual ScalarFunctionInterface<P>* _derivative(uint j) const {
        ARIADNE_NOT_IMPLEMENTED;
    }

    virtual ScalarFunction<P> derivative(uint j) const {
        if(_arg2==0) { return ScalarFunction<P>::constant(this->argument_size(),X(0)); }
        if(_arg2==1) { return _arg1.derivative(j); }
        if(_arg2==2) { return 2*_arg1.derivative(j)*_arg1; }
        return _arg2*_arg1.derivative(j)*pow(_arg1,_arg2-1);
    }

    virtual std::ostream& repr(std::ostream& os) const {
        return os << "PF["<<this->argument_size()<<"]("<< *this <<")"; }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "pow(" << _arg1 << "," << _arg2 << ")"; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=pow(_arg1.evaluate(x),_arg2); }

    OperatorCode _op;
    ScalarFunction<P> _arg1;
    Int _arg2;
};

typedef ConstantFunction<Real> RealConstantFunction;
typedef ConstantFunction<EffectiveNumber> EffectiveConstantFunction;
typedef CoordinateFunction<EffectiveTag> EffectiveCoordinateFunction;
typedef UnaryFunction<EffectiveTag> EffectiveUnaryFunction;
typedef BinaryFunction<EffectiveTag> EffectiveBinaryFunction;
typedef PowerFunction<EffectiveTag> EffectivePowerFunction;


//------------------------ Results of functional operations  -----------------------------------//

template<class P>
struct ScalarEmbeddedFunction
    : ScalarFunctionMixin<ScalarEmbeddedFunction<P>,P>
{
    ScalarEmbeddedFunction(Nat as1, const ScalarFunction<P>& f2, Nat as3)
        : _as1(as1), _f2(f2), _as3(as3) { }
    virtual Nat argument_size() const { return _as1+_f2.argument_size()+_as3; }
    virtual ScalarFunctionInterface<P>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { return os << "ScalarEmbeddedFunction( as1="<<_as1<<", f2="<<_f2<<", as3="<<_as3<<" )"; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        Vector<XX> px=project(x,Range(_as1,_as1+_f2.argument_size())); r=_f2.evaluate(px); }

    Nat _as1;
    ScalarFunction<P> _f2;
    Nat _as3;
};


template<class P>
struct VectorEmbeddedFunction
    : VectorFunctionMixin<VectorEmbeddedFunction<P>,P>
{
    VectorEmbeddedFunction(Nat as1, const VectorFunction<P>& f2, Nat as3)
        : _as1(as1), _f2(f2), _as3(as3) { }
    virtual Nat result_size() const { return _f2.result_size(); }
    virtual Nat argument_size() const { return _as1+_f2.argument_size()+_as3; }
    virtual ScalarFunctionInterface<P>* _get(uint j) const { return new ScalarEmbeddedFunction<P>(_as1,_f2.get(j),_as3); }
    virtual std::ostream& write(std::ostream& os) const { return os << "VectorEmbeddedFunction( as1="<<_as1<<", f2="<<_f2<<", as3="<<_as3<<" )"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        Vector<XX> px=project(x,Range(_as1,_as1+_f2.argument_size())); r=_f2.evaluate(px); }

    Nat _as1;
    VectorFunction<P> _f2;
    Nat _as3;
};


template<class P>
struct ScalarComposedFunction
    : ScalarFunctionMixin<ScalarComposedFunction<P>,P>
{
    ScalarComposedFunction(const ScalarFunction<P>& f, const VectorFunction<P>& g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual Nat argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<P>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { return os << "ScalarComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    ScalarFunction<P> _f;
    VectorFunction<P> _g;
};


template<class P>
struct VectorComposedFunction
    :  VectorFunctionMixin<VectorComposedFunction<P>,P>
{
    VectorComposedFunction(VectorFunction<P> f, VectorFunction<P> g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual Nat result_size() const { return _f.result_size(); }
    virtual Nat argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<P>* _get(uint i) const { return new ScalarComposedFunction<P>(_f[i],_g); }
    virtual std::ostream& write(std::ostream& os) const { return os << "ComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    VectorFunction<P> _f;
    VectorFunction<P> _g;
};


template<class P>
struct JoinedFunction
    : VectorFunctionMixin<JoinedFunction<P>,P>
{
    JoinedFunction(VectorFunction<P> f1, VectorFunction<P> f2)
        : _f1(f1), _f2(f2) { ARIADNE_ASSERT(f1.argument_size()==f2.argument_size()); }
    virtual Nat result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual Nat argument_size() const { return _f1.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os << "JoinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }
    virtual ScalarFunctionInterface<P>* _get(uint i) const { return (i<_f1.result_size()) ? _f1.raw_pointer()->_get(i) : _f2.raw_pointer()->_get(i-_f1.result_size()); }
    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        r=join(_f1.evaluate(x),_f2.evaluate(x)); }

    VectorFunction<P> _f1;
    VectorFunction<P> _f2;
};


template<class P>
class CombinedFunction
    : VectorFunctionMixin<CombinedFunction<P>,P>
{
    CombinedFunction(VectorFunction<P> f1, VectorFunction<P> f2)
        : _f1(f1), _f2(f2) { }
    virtual Nat result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual Nat argument_size() const { return _f1.argument_size()+_f2.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os << "CombinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        return r=combine(_f1.evaluate(project(x,range(0,_f1.argument_size()))),
                         _f2.evaluate(project(x,range(_f1.argument_size(),this->argument_size())))); }

    VectorFunction<P> _f1;
    VectorFunction<P> _f2;
};


// A Lie deriviative \f$\nabla g\cdot f\f$.
template<class P>
struct LieDerivativeFunction
    : ScalarFunctionMixin<LieDerivativeFunction<P>,P>
{
    //! \brief Construct the identity function in dimension \a n.
    LieDerivativeFunction(const ScalarFunction<P>& g, const VectorFunction<P>& f) {
        ARIADNE_ASSERT(g.argument_size()==f.argument_size());
        ARIADNE_ASSERT(f.result_size()==f.argument_size());
        _g=g; for(uint j=0; j!=g.argument_size(); ++j) { _dg[j]=g.derivative(j); } _f=f; }
    Nat argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<P>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    std::ostream& write(std::ostream& os) const { return os << "LieDerivative( g="<<_g<<", f="<<_f<<" )"; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        //const Vector<R> fx=_f.evaluate(x); r=0; for(uint i=0; i!=_dg.size(); ++i) { r+=fx[i]+_dg[i].evaluate(x); } }
        Vector<XX> fx=_f.evaluate(x);
        r=0;
        for(uint i=0; i!=_dg.size(); ++i) {
            r+=fx[i]+_dg[i].evaluate(x);
        }
    }

    ScalarFunction<P> _g;
    List< ScalarFunction<P> > _dg;
    VectorFunction<P> _f;
};




} // namespace Ariadne

#endif
