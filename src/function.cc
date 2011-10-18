/***************************************************************************
 *            function.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "numeric.h"
#include "differential.h"
#include "function.h"
#include "operators.h"
#include "formula.h"
#include "algebra.h"

#include "function_mixin.h"
#include "function_mixin.tcc"

#include "function_model.h"

#include "taylor_function.h"


namespace Ariadne {

template<class T> inline std::string str(const T& t) {
    std::stringstream ss; ss << t; return ss.str(); }

typedef uint Nat;
typedef int Int;

// Functions for converting from Expression classes to Function classes via Formula
template<class X> class Expression;
template<class X> class Variable;
template<class X> class Space;
Nat dimension(const Space<Real>& spc);
Nat len(const List< Variable<Real> >& vars);
Formula<Real> formula(const Expression<Real>& expr, const Space<Real>& spc);
Formula<Real> formula(const Expression<Real>& expr, const List< Variable<Real> >& vars);
ScalarFunction<Real> make_function(const Expression<Real>& expr, const Space<Real>& spc) {
    return ScalarFunction<Real>(dimension(spc),formula(expr,spc)); }
ScalarFunction<Real> make_function(const Expression<Real>& expr, const List< Variable<Real> >& vars) {
    return ScalarFunction<Real>(len(vars),formula(expr,vars)); }



//------------------------ FunctionInterface forwarded functions  -----------------------------------//

Vector<Float> ScalarFunctionInterface<Float>::gradient(const Vector<Float>& x) const {
    return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
Vector<Interval> ScalarFunctionInterface<Interval>::gradient(const Vector<Interval>& x) const {
    return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }
Matrix<Float> VectorFunctionInterface<Float>::jacobian(const Vector<Float>& x) const {
    return this->evaluate(Differential<Float>::variables(1u,x)).jacobian(); }
Matrix<Interval> VectorFunctionInterface<Interval>::jacobian(const Vector<Interval>& x) const {
    return this->evaluate(Differential<Interval>::variables(1u,x)).jacobian(); }

//------------------------ Formula functions  -----------------------------------//

//! A function defined by a formula
template<class X>
struct ScalarFormulaFunction
    : ScalarFunctionMixin<ScalarFormulaFunction<X>,X>
{
    Nat _argument_size;
    Formula<X> _formula;

    ScalarFormulaFunction(Nat as, const Formula<X>& f) : _argument_size(as), _formula(f) { }
    operator Formula<X>() const { return _formula; }

    virtual Nat argument_size() const { return _argument_size; }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const { return new ScalarFormulaFunction<X>(_argument_size,Ariadne::derivative(_formula,j)); }
    virtual std::ostream& write(std::ostream& os) const { return os << this->_formula; }
    virtual std::ostream& repr(std::ostream& os) const { return os << "FF["<<this->_argument_size<<"]("<<this->_formula<<")"; }
    template<class Y> void _compute(Y& r, const Vector<Y>& x) const { r=Ariadne::evaluate(_formula,x); }
};

inline ScalarFunction<Real> function(Nat n, Formula<Real> f) { return new ScalarFormulaFunction<Real>(n,f); }

typedef ScalarFormulaFunction<Real> RealScalarFormulaFunction;

//! A vector function defined by formulae
template<class X>
struct VectorFormulaFunction
    : VectorFunctionMixin<VectorFormulaFunction<X>,X>
{
    Nat _argument_size;
    Vector< Formula<X> > _formulae;

    VectorFormulaFunction(uint as, const List< Formula<X> >& f) : _argument_size(as), _formulae(f) { }
    VectorFormulaFunction(uint as, const Vector< Formula<X> >& f) : _argument_size(as), _formulae(f) { }

    virtual Nat result_size() const { return this->_formulae.size(); }
    virtual Nat argument_size() const { return this->_argument_size; }
    virtual ScalarFormulaFunction<X>* _get(uint i) const { return new ScalarFormulaFunction<X>(this->_argument_size,this->_formulae[i]); }
    virtual std::ostream& write(std::ostream& os) const { return os << "FF[R"<<this->argument_size()<<"->R"<<this->result_size()<<"]("<<this->_formulae<<")"; }
    template<class Y> void _compute(Vector<Y>& r, const Vector<Y>& x) const { r=Ariadne::cached_evaluate(this->_formulae,x); }
};



//------------------------ Arithmetic scalar functions  -----------------------------------//


//! A constant function f(x)=c
template<class X>
struct ConstantFunction
    : ScalarFunctionMixin<ConstantFunction<X>,X>
{
    Nat _argument_size;
    X _value;

    ConstantFunction(uint as, const X& c) : _argument_size(as), _value(c) { }
    operator X() const { return _value; }

    virtual Nat argument_size() const { return _argument_size; }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const { return new ConstantFunction<X>(_argument_size,0.0); }
    virtual std::ostream& repr(std::ostream& os) const { return os << this->_value; }
    virtual std::ostream& write(std::ostream& os) const { return os << "CF["<<this->_argument_size<<"]("<<_value<<")"; }
    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=x.zero_element()+_value; }
};

//! A coordinate function \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=x_i\f$.
template<class X>
struct CoordinateFunction
    : ScalarFunctionMixin<CoordinateFunction<X>,X>
{
    Nat _argument_size;
    Nat _index;

    CoordinateFunction(uint as, uint i) : _argument_size(as), _index(i) { }
    Nat index() const { return _index; }

    virtual Nat argument_size() const { return _argument_size; }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const {
        if(j==_index) { return new CoordinateFunction<X>(_argument_size,1.0); }
        else { return new CoordinateFunction<X>(_argument_size,0.0); } }
    virtual std::ostream& repr(std::ostream& os) const { return os << "x"<<this->_index; }
    virtual std::ostream& write(std::ostream& os) const { return os << "IF["<<this->_argument_size<<"](x"<<this->_index<<")"; }
    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=x[_index]; }
};


//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
template<class X>
struct UnaryFunction
    : ScalarFunctionMixin< UnaryFunction<X>, X >
{
  public:
    UnaryFunction(const OperatorCode& op, const ScalarFunction<X>& arg)
        : _op(op), _arg(arg) { }
    virtual UnaryFunction<X>* clone() const { return new UnaryFunction<X>(*this); }
    virtual Nat argument_size() const {
        return this->_arg.argument_size(); }

    virtual ScalarFunctionInterface<X>* _derivative(uint j) const {
        return static_cast<const ScalarFunctionInterface<X>&>(this->derivative(j))._clone();
    }

    virtual ScalarFunction<X> derivative(uint j) const {
        switch(_op) {
            case POS: return _arg.derivative(j);
            case NEG: return -_arg.derivative(j);
            case REC: return -_arg.derivative(j)/sqr(_arg);
            case SQR: return 2.0*_arg.derivative(j)*_arg;
            case SQRT: return _arg.derivative(j)/(2.0*sqrt(_arg));
            case EXP: return _arg*_arg.derivative(j);
            case LOG: return _arg.derivative(j)/_arg;
            case SIN: return _arg.derivative(j)*cos(_arg);
            case COS: return -_arg.derivative(j)*sin(_arg);
            default: ARIADNE_FAIL_MSG("Unknown unary function "<<this->_op);
        }
    }

    virtual std::ostream& repr(std::ostream& os) const {
        return os << "UF[" << this->argument_size() << "](" << *this << ")"; }
    virtual std::ostream& write(std::ostream& os) const {
        return os << _op << '(' << _arg << ')'; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=Ariadne::compute(_op,_arg.evaluate(x)); }

    OperatorCode _op;
    ScalarFunction<X> _arg;
};


template<class X> ScalarFunction<X> sqr(const ScalarFunction<X>& f) {
    return new UnaryFunction<X>(SQR,f); }
template<class X> ScalarFunction<X> sqrt(const ScalarFunction<X>& f) {
    return new UnaryFunction<X>(SQRT,f); }
template<class X> ScalarFunction<X> sin(const ScalarFunction<X>& f) {
    return new UnaryFunction<X>(SIN,f); }
template<class X> ScalarFunction<X> cos(const ScalarFunction<X>& f) {
    return new UnaryFunction<X>(COS,f); }

template<class X>
struct BinaryFunction
    : ScalarFunctionMixin< BinaryFunction<X>, X >
{
  public:
    BinaryFunction(OperatorCode op, const ScalarFunction<X>& arg1, const ScalarFunction<X>& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) { ARIADNE_ASSERT_MSG(arg1.argument_size()==arg2.argument_size(),"op='"<<op<<"', arg1="<<arg1<<", arg2="<<arg2); }
    virtual BinaryFunction<X>* clone() const { return new BinaryFunction<X>(*this); }
    virtual Nat argument_size() const {
        return this->_arg1.argument_size(); }

    virtual ScalarFunctionInterface<X>* _derivative(uint j) const {
        return static_cast<const ScalarFunctionInterface<X>&>(this->derivative(j))._clone();
    }

    virtual ScalarFunction<X> derivative(uint j) const {
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
        return os << "BF[" << this->argument_size() << "](" << *this << ")"; }
    virtual std::ostream& write(std::ostream& os) const {
        if(_op==ADD || _op==SUB) { return os << '(' << _arg1 << symbol(_op) << _arg2 << ')'; }
        else { return os << _arg1 << symbol(_op) << _arg2; } }
        
    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=Ariadne::compute(_op,_arg1.evaluate(x),_arg2.evaluate(x)); }

    OperatorCode _op;
    ScalarFunction<X> _arg1;
    ScalarFunction<X> _arg2;
};


// \brief The power function \f$(x,n)\mapsto x^n\f$.
template<class X>
class PowerFunction
    : public ScalarFunctionMixin< PowerFunction<X>, X >
{
  public:
    PowerFunction(OperatorCode op, const ScalarFunction<X>& arg1, const Int& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {  }
    virtual PowerFunction<X>* clone() const { return new PowerFunction<X>(*this); }
    virtual Nat argument_size() const {
        return this->_arg1.argument_size(); }

    virtual ScalarFunctionInterface<X>* _derivative(uint j) const {
        ARIADNE_NOT_IMPLEMENTED;
    }

    virtual ScalarFunction<X> derivative(uint j) const {
        if(_arg2==0) { return ScalarFunction<X>::constant(this->argument_size(),0.0); }
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
    ScalarFunction<X> _arg1;
    Int _arg2;
};

typedef ConstantFunction<Real> RealConstantFunction;
typedef CoordinateFunction<Real> RealCoordinateFunction;
typedef UnaryFunction<Real> RealUnaryFunction;
typedef BinaryFunction<Real> RealBinaryFunction;
typedef PowerFunction<Real> RealPowerFunction;









//------------------------ Vector of Scalar functions  -----------------------------------//

template<class X>
struct VectorOfScalarFunction
    : VectorFunctionMixin<VectorOfScalarFunction<X>,X>
{
    VectorOfScalarFunction(uint rs, uint as)
        : _as(as), _vec(rs,ScalarFunction<X>(as)) { }
    VectorOfScalarFunction(uint rs, const ScalarFunction<X>& f)
        : _as(f.argument_size()), _vec(rs,f) { }

    void set(uint i, const ScalarFunction<X>& f) {
        if(this->argument_size()==0u) { this->_as=f.argument_size(); }
        ARIADNE_ASSERT(f.argument_size()==this->argument_size());
        this->_vec[i]=f; }
    ScalarFunction<X> get(uint i) const {
        return this->_vec[i]; }

    virtual uint result_size() const {
        return _vec.size(); }
    virtual uint argument_size() const {
        return _as; }

    virtual ScalarFunctionInterface<X>* _get(uint i) const { return static_cast<const ScalarFunctionInterface<X>&>(this->_vec[i])._clone(); }

    ScalarFunction<X> operator[](uint i) const {
        return this->_vec[i]; }

    virtual std::ostream& write(std::ostream& os) const {
        os << "VoSF[" << this->result_size() << "," << this->argument_size() << "][";
        for(uint i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->repr(os); }
        return os << "]"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        r=Vector<XX>(this->_vec.size());
        for(uint i=0; i!=r.size(); ++i) {
            r[i]=_vec[i].evaluate(x); } }

    uint _as;
    Vector< ScalarFunction<X> > _vec;

};


template<class X>
struct FunctionElement
    : ScalarFunctionMixin<FunctionElement<X>,X>
{
    FunctionElement(const VectorFunction<X>& f, Nat i)
        : _f(f), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }

    virtual Nat argument_size() const { return _f.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os<<_f<<"["<<_i<<"]"; }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=this->_f.evaluate(x)[_i]; }

    VectorFunction<X> _f;
    Nat _i;
};



//------------------------ Results of functional operations  -----------------------------------//

template<class X>
struct ScalarEmbeddedFunction
    : ScalarFunctionMixin<ScalarEmbeddedFunction<X>,X>
{
    ScalarEmbeddedFunction(Nat as1, const ScalarFunction<X>& f2, Nat as3)
        : _as1(as1), _f2(f2), _as3(as3) { }
    virtual Nat argument_size() const { return _as1+_f2.argument_size()+_as3; }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { return os << "ScalarEmbeddedFunction( as1="<<_as1<<", f2="<<_f2<<", as3="<<_as3<<" )"; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        Vector<XX> px=project(x,Range(_as1,_as1+_f2.argument_size())); r=_f2.evaluate(px); }

    Nat _as1;
    ScalarFunction<X> _f2;
    Nat _as3;
};


template<class X>
struct VectorEmbeddedFunction
    : VectorFunctionMixin<VectorEmbeddedFunction<X>,X>
{
    VectorEmbeddedFunction(Nat as1, const VectorFunction<X>& f2, Nat as3)
        : _as1(as1), _f2(f2), _as3(as3) { }
    virtual Nat result_size() const { return _f2.result_size(); }
    virtual Nat argument_size() const { return _as1+_f2.argument_size()+_as3; }
    virtual ScalarFunctionInterface<X>* _get(uint j) const { return new ScalarEmbeddedFunction<X>(_as1,_f2.get(j),_as3); }
    virtual std::ostream& write(std::ostream& os) const { return os << "VectorEmbeddedFunction( as1="<<_as1<<", f2="<<_f2<<", as3="<<_as3<<" )"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        Vector<XX> px=project(x,Range(_as1,_as1+_f2.argument_size())); r=_f2.evaluate(px); }

    Nat _as1;
    VectorFunction<X> _f2;
    Nat _as3;
};


template<class X>
struct ScalarComposedFunction
    : ScalarFunctionMixin<ScalarComposedFunction<X>,X>
{
    ScalarComposedFunction(const ScalarFunction<X>& f, const VectorFunction<X>& g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual Nat argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { return os << "ScalarComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    ScalarFunction<X> _f;
    VectorFunction<X> _g;
};


template<class X>
struct VectorComposedFunction
    :  VectorFunctionMixin<VectorComposedFunction<X>,X>
{
    VectorComposedFunction(VectorFunction<X> f, VectorFunction<X> g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual Nat result_size() const { return _f.result_size(); }
    virtual Nat argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<X>* _get(uint i) const { return new ScalarComposedFunction<X>(_f[i],_g); }
    virtual std::ostream& write(std::ostream& os) const { return os << "ComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    VectorFunction<X> _f;
    VectorFunction<X> _g;
};


template<class X>
struct JoinedFunction
    : VectorFunctionMixin<JoinedFunction<X>,X>
{
    JoinedFunction(VectorFunction<X> f1, VectorFunction<X> f2)
        : _f1(f1), _f2(f2) { ARIADNE_ASSERT(f1.argument_size()==f2.argument_size()); }
    virtual Nat result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual Nat argument_size() const { return _f1.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os << "JoinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }
    virtual ScalarFunctionInterface<X>* _get(uint i) const { return (i<_f1.result_size()) ? _f1.raw_pointer()->_get(i) : _f2.raw_pointer()->_get(i-_f1.result_size()); }
    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        r=join(_f1.evaluate(x),_f2.evaluate(x)); }

    VectorFunction<X> _f1;
    VectorFunction<X> _f2;
};


template<class X>
class CombinedFunction
    : VectorFunctionMixin<CombinedFunction<X>,X>
{
    CombinedFunction(VectorFunction<X> f1, VectorFunction<X> f2)
        : _f1(f1), _f2(f2) { }
    virtual Nat result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual Nat argument_size() const { return _f1.argument_size()+_f2.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os << "CombinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        return r=combine(_f1.evaluate(project(x,range(0,_f1.argument_size()))),
                         _f2.evaluate(project(x,range(_f1.argument_size(),this->argument_size())))); }

    VectorFunction<X> _f1;
    VectorFunction<X> _f2;
};


// A Lie deriviative \f$\nabla g\cdot f\f$.
template<class X>
struct LieDerivativeFunction
    : ScalarFunctionMixin<LieDerivativeFunction<X>,X>
{
    //! \brief Construct the identity function in dimension \a n.
    LieDerivativeFunction(const ScalarFunction<X>& g, const VectorFunction<X>& f) {
        ARIADNE_ASSERT(g.argument_size()==f.argument_size());
        ARIADNE_ASSERT(f.result_size()==f.argument_size());
        _g=g; for(uint j=0; j!=g.argument_size(); ++j) { _dg[j]=g.derivative(j); } _f=f; }
    Nat argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    std::ostream& write(std::ostream& os) const { return os << "LieDerivative( g="<<_g<<", f="<<_f<<" )"; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        //const Vector<R> fx=_f.evaluate(x); r=0; for(uint i=0; i!=_dg.size(); ++i) { r+=fx[i]+_dg[i].evaluate(x); } }
        Vector<XX> fx=_f.evaluate(x);
        r=0;
        for(uint i=0; i!=_dg.size(); ++i) {
            r+=fx[i]+_dg[i].evaluate(x);
        }
    }

    ScalarFunction<X> _g;
    List< ScalarFunction<X> > _dg;
    VectorFunction<X> _f;
};




//------------------------ Conversion to a formula -----------------------------------//

Formula<Real> formula(const ScalarFunction<Real>& f) {
    const ScalarFunctionInterface<Real>* fptr=f.raw_pointer();
    const ScalarFormulaFunction<Real>* ff=dynamic_cast<const ScalarFormulaFunction<Real>*>(fptr);
    if(ff) { return ff->_formula; }
    const RealConstantFunction* cf=dynamic_cast<const RealConstantFunction*>(fptr);
    const RealCoordinateFunction* indf=dynamic_cast<const RealCoordinateFunction*>(fptr);
    const RealUnaryFunction* uf=dynamic_cast<const RealUnaryFunction*>(fptr);
    const RealBinaryFunction* bf=dynamic_cast<const RealBinaryFunction*>(fptr);
    const RealPowerFunction* pf=dynamic_cast<const RealPowerFunction*>(fptr);
    if(cf) { return Formula<Real>::constant(cf->_value); }
    if(indf) { return Formula<Real>::coordinate(indf->_index); }
    if(uf) { return make_formula(uf->_op,formula(uf->_arg)); }
    if(bf) { return make_formula(bf->_op,formula(bf->_arg1),formula(bf->_arg2)); }
    if(pf) { return make_formula(pf->_op,formula(pf->_arg1),pf->_arg2); }
    ARIADNE_FAIL_MSG("Cannot compute formula for function "<<f<<"\n");
}

Vector< Formula<Real> > formula(const VectorFunction<Real>& f) {
    const VectorFunctionInterface<Real>& fi=f;
    const VectorFormulaFunction<Real>* ff;
    const VectorOfScalarFunction<Real>* vf;
    if( (vf=dynamic_cast<const VectorOfScalarFunction<Real>*>(&fi)) ) {
        Vector< Formula<Real> > r(vf->result_size());
        for(uint i=0; i!=r.size(); ++i) { r[i]=formula((*vf)[i]); }
        return r;
    } else if( (ff=dynamic_cast<const VectorFormulaFunction<Real>*>(&fi)) ) {
        return ff->_formulae;
    } else {
        ARIADNE_FAIL_MSG("Cannot compute formula for function "<<f<<"\n");
    }
}





//------------------------ ScalarFunction -----------------------------------//

template<class X> ScalarFunction<X>::ScalarFunction()
    : _ptr(new ScalarFormulaFunction<X>(0u,Formula<X>::zero()))
{
}

template<class X> ScalarFunction<X>::ScalarFunction(Nat n)
    : _ptr(new ScalarFormulaFunction<X>(n,Formula<X>::zero()))
{
}

template<class X> ScalarFunction<X>::ScalarFunction(Nat n, Formula<X> e)
    : _ptr(new ScalarFormulaFunction<X>(n,e))
{
}

template<class X> ScalarFunction<X> ScalarFunction<X>::zero(Nat n)
{
    return new ScalarFormulaFunction<X>(n,Formula<X>::zero());
}

template<class X> ScalarFunction<X> ScalarFunction<X>::constant(Nat n, X c)
{
    return new ScalarFormulaFunction<X>(n,Formula<X>::constant(c));
}

template<class X> ScalarFunction<X> ScalarFunction<X>::coordinate(Nat n, Nat j)
{
    return new ScalarFormulaFunction<X>(n,Formula<X>::coordinate(j));
}

template<class X> List< ScalarFunction<X> > ScalarFunction<X>::coordinates(Nat n)
{
    List< ScalarFunction<X> > lsf; lsf.reserve(n);
    for(Nat i=0; i!=n; ++i) { lsf.append(ScalarFunction<X>::coordinate(n,i)); }
    return lsf;
}

template class ScalarFunction<Float>;
template class ScalarFunction<Interval>;
template class ScalarFunction<Real>;



//------------------------ Scalar arithmetic operators -----------------------------------//





template<class OP> inline
RealScalarFunction make_unary_function(OP op, const RealScalarFunction& f) {
    return RealScalarFunction(new RealUnaryFunction(op.code(),f)); }

template<class OP> inline
RealScalarFunction make_binary_function(OP op, const RealScalarFunction& f1, const RealScalarFunction& f2) {
    return RealScalarFunction(new RealBinaryFunction(op.code(),f1,f2)); }

template<class OP> inline
RealScalarFunction make_binary_function(OP op, const RealScalarFunction& f1, const Int& n2) {
    return RealScalarFunction(new RealPowerFunction(op.code(),f1,n2)); }


RealScalarFunction operator+(const RealScalarFunction& f)
{
    return f;
}

RealScalarFunction operator-(const RealScalarFunction& f)
{
    const RealScalarFormulaFunction* ff=dynamic_cast<const RealScalarFormulaFunction*>(f.raw_pointer());
    if(ff) { return function(ff->_argument_size,-ff->_formula); }
    return make_unary_function(Neg(),f);
}

RealScalarFunction operator+(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula+e2->_formula);
    }
    return make_binary_function(Add(),f1,f2);
}

RealScalarFunction operator-(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula-e2->_formula);
    }
    return make_binary_function(Sub(),f1,f2);
}

RealScalarFunction operator*(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula*e2->_formula);
    }
    return make_binary_function(Mul(),f1,f2);
}

RealScalarFunction operator/(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula/e2->_formula);
    }
    return make_binary_function(Div(),f1,f2);
}


RealScalarFunction operator+(const RealScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula+s2); }
    return f1+RealScalarFunction::constant(f1.argument_size(),s2);
}

RealScalarFunction operator-(const RealScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula-s2); }
    return f1-RealScalarFunction::constant(f1.argument_size(),s2);
}

RealScalarFunction operator*(const RealScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula*s2); }
    return f1*RealScalarFunction::constant(f1.argument_size(),s2);
}

RealScalarFunction operator/(const RealScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula/s2); }
    return f1/RealScalarFunction::constant(f1.argument_size(),s2);
}

RealScalarFunction operator/(const RealScalarFunction& f1, int s2)
{
    return f1/Real(s2);
}

RealScalarFunction operator+(const Real& s1, const RealScalarFunction& f2)
{
    return f2+s1;
}

RealScalarFunction operator-(const Real& s1, const RealScalarFunction& f2)
{
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e2) { return function(e2->_argument_size,s1-e2->_formula); }
    return RealScalarFunction::constant(f2.argument_size(),s1)-f2;
}

RealScalarFunction operator*(const Real& s1, const RealScalarFunction& f2)
{
    return f2*s1;
}

RealScalarFunction operator/(const Real& s1, const RealScalarFunction& f2)
{
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e2) { return function(e2->_argument_size,s1/e2->_formula); }
    return RealScalarFunction::constant(f2.argument_size(),s1)/f2;
}

RealScalarFunction pow(const RealScalarFunction& f, Nat m)
{
    const RealScalarFormulaFunction* e=dynamic_cast<const RealScalarFormulaFunction*>(f.raw_pointer());
    if(e) { return function(e->_argument_size,pow(e->_formula,m)); }
    return make_binary_function(Pow(),f,m);
}


RealScalarFunction pow(const RealScalarFunction& f, Int n)
{
    const RealScalarFormulaFunction* e=dynamic_cast<const RealScalarFormulaFunction*>(f.raw_pointer());
    if(e) { return function(e->_argument_size, pow(e->_formula,n)); }
    return make_binary_function(Pow(),f,n);
}

RealScalarFunction neg(const RealScalarFunction& f) {
    return make_unary_function(Neg(),f); }

RealScalarFunction rec(const RealScalarFunction& f) {
    return make_unary_function(Rec(),f); }

RealScalarFunction sqr(const RealScalarFunction& f) {
    return make_unary_function(Sqr(),f); }

RealScalarFunction sqrt(const RealScalarFunction& f) {
    return make_unary_function(Sqrt(),f); }

RealScalarFunction exp(const RealScalarFunction& f) {
    return make_unary_function(Exp(),f); }

RealScalarFunction log(const RealScalarFunction& f) {
    return make_unary_function(Log(),f); }

RealScalarFunction sin(const RealScalarFunction& f) {
    return make_unary_function(Sin(),f); }

RealScalarFunction cos(const RealScalarFunction& f) {
    return make_unary_function(Cos(),f); }

RealScalarFunction tan(const RealScalarFunction& f) {
    return make_unary_function(Tan(),f); }






//------------------------ Vector Function ----------------------------------//


template<class X> VectorFunction<X>::VectorFunction()
    : _ptr(new VectorOfScalarFunction<X>(0u,ScalarFunction<X>()))
{
}

template<class X> VectorFunction<X>::VectorFunction(Nat rs, Nat as)
    : _ptr(new VectorOfScalarFunction<X>(rs,ScalarFunction<X>::zero(as)))
{
}

template<class X> VectorFunction<X>::VectorFunction(const List< ScalarFunction<X> >& lsf)
{
    ARIADNE_ASSERT(lsf.size()>0);
    Nat as=lsf[0].argument_size();
    VectorOfScalarFunction<X>* new_ptr=new VectorOfScalarFunction<X>(lsf.size(),as);
    for(uint i=0; i!=lsf.size(); ++i) {
        new_ptr->set(i,lsf[i]);
    }
    this->_ptr=shared_ptr< const VectorFunctionInterface<X> >(new_ptr);
}

template<class X> VectorFunction<X>::VectorFunction(Nat as, const List< Formula<X> >& le)
{
    ARIADNE_ASSERT(le.size()>0);
    VectorOfScalarFunction<X>* new_ptr=new VectorOfScalarFunction<X>(le.size(),as);
    for(uint i=0; i!=le.size(); ++i) {
        new_ptr->set(i,ScalarFormulaFunction<X>(as,le[i]));
    }
    this->_ptr=shared_ptr< const VectorFunctionInterface<X> >(new_ptr);
}

template<class X> VectorFunction<X>::VectorFunction(Nat rs, ScalarFunction<X> sf)
    : _ptr(new VectorOfScalarFunction<X>(rs,sf))
{
}

template<class X> VectorFunction<X> VectorFunction<X>::zeros(Nat rs, Nat as)
{
    VectorOfScalarFunction<X>* res = new VectorOfScalarFunction<X>(rs,as);
    for(uint i=0; i!=rs; ++i) {
        res->_vec[i]=ScalarFunction<X>::zero(as);
    }
    return res;
}

template<class X> VectorFunction<X> VectorFunction<X>::identity(Nat n)
{
    VectorOfScalarFunction<X>* res = new VectorOfScalarFunction<X>(n,n);
    for(uint i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<X>::coordinate(n,i);
    }
    return res;
}

template<class X> Void VectorFunction<X>::set(Nat i, ScalarFunction<X> sf)
{
    if(!this->_ptr.unique()) {
        ARIADNE_WARN("VectorFunction<X>::set(Nat, ScalarFunction<X>): Cloning shared pointer.\n");
        this->_ptr=shared_ptr< const VectorFunctionInterface<X> >(this->_ptr->_clone());
    }
    VectorFunctionInterface<X>* vf=const_cast<VectorFunctionInterface<X>*>(this->raw_pointer());
    VectorOfScalarFunction<X>* vsf=dynamic_cast<VectorOfScalarFunction<X>*>(vf);
    if(vsf) {
        vsf->set(i,sf);
    } else {
        ARIADNE_THROW(std::runtime_error,"Void VectorFunction<X>::set(Nat i, ScalarFunction<X> sf)","Cannot set element of "<<*vf<<"\n");
    }
}

template class VectorFunction<Float>;
template class VectorFunction<Interval>;
template class VectorFunction<Real>;



//------------------------ Vector function operators -------------------------------//

RealVectorFunction operator*(const RealScalarFunction& f, const Vector<Real>& e) {
    for(uint i=0; i!=e.size(); ++i) { ARIADNE_ASSERT(e[i]==Real(0) || e[i]==Real(1)); }
    RealVectorFunction r(e.size(),f.argument_size());
    for(uint i=0; i!=e.size(); ++i) {
        if(e[i]==Real(1)) { r.set(i,f); }
    }
    return r;
}

RealVectorFunction operator+(const RealVectorFunction& f1, const RealVectorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    RealVectorFunction r(f1.result_size(),f1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,f1[i]+f2[i]);
    }
    return r;
}

RealVectorFunction operator-(const RealVectorFunction& f1, const RealVectorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    RealVectorFunction r(f1.result_size(),f1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,f1[i]-f2[i]);
    }
    return r;
}

VectorFunction<Real> operator*(const VectorFunction<Real>& vf, const ScalarFunction<Real>& sf) {
    ARIADNE_ASSERT(vf.argument_size()==sf.argument_size());
    RealVectorFunction r(vf.result_size(),vf.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,vf[i]*sf);
    }
    return r;
}

VectorFunction<Real> operator*(const ScalarFunction<Real>& sf, const VectorFunction<Real>& vf) {
    ARIADNE_ASSERT(sf.argument_size()==vf.argument_size());
    RealVectorFunction r(vf.result_size(),vf.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,sf*vf[i]);
    }
    return r;
}



RealVectorFunction join(const RealScalarFunction& f1, const RealScalarFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    RealVectorFunction r(2,f1.argument_size());
    r.set(0,f1);
    r.set(1,f2);
    return r;
}

RealVectorFunction join(const RealVectorFunction& f1, const RealScalarFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    RealVectorFunction r(f1.result_size()+1u,f1.argument_size());
    for(uint i=0u; i!=f1.result_size(); ++i) {
        r.set(i,f1.get(i));
    }
    r.set(f1.result_size(),f2);
    return r;
}

RealVectorFunction join(const RealScalarFunction& f1, const RealVectorFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    RealVectorFunction r(f2.result_size()+1u,f1.argument_size());
    r.set(0u,f1);
    for(uint i=0u; i!=f2.result_size(); ++i) {
        r.set(i+1u,f2.get(i));
    }
    return r;
}

RealVectorFunction join(const RealVectorFunction& f1, const RealVectorFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    RealVectorFunction r(f1.result_size()+f2.result_size(),f1.argument_size());
    for(uint i=0u; i!=f1.result_size(); ++i) {
        r.set(i,f1.get(i));
    }
    for(uint i=0u; i!=f2.result_size(); ++i) {
        r.set(i+f1.result_size(),f2.get(i));
    }
    return r;
}

RealScalarFunction embed(Nat as1, const RealScalarFunction& f2, Nat as3) {
    return new ScalarEmbeddedFunction<Real>(as1,f2,as3);
}

RealVectorFunction embed(Nat as1, const RealVectorFunction& f2, Nat as3) {
    return new VectorEmbeddedFunction<Real>(as1,f2,as3);
}

RealScalarFunction compose(const RealScalarFunction& f, const RealVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return new ScalarComposedFunction<Real>(f,g);
}

RealVectorFunction compose(const RealVectorFunction& f, const RealVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return new VectorComposedFunction<Real>(f,g);
}

RealScalarFunction lie_derivative(const RealScalarFunction& g, const RealVectorFunction& f) {
    ARIADNE_ASSERT_MSG(g.argument_size()==f.result_size(),"f="<<f<<", g="<<g<<"\n");
    ARIADNE_ASSERT_MSG(f.result_size()==f.argument_size(),"f="<<f<<", g="<<g<<"\n");
    ARIADNE_ASSERT_MSG(f.result_size()>0,"f="<<f<<", g="<<g<<"\n");

    try {
        RealScalarFunction r=g.derivative(0)*f[0];
        for(uint i=1; i!=g.argument_size(); ++i) {
            r=r+g.derivative(i)*f[i];
        }
        return r;
    }
    catch(...) {
        ARIADNE_FAIL_MSG("Failed to compute Lie derivative of "<<g<<" under vector field "<<f<<"\n");
    }
}



//------------------------ Interval function operators -------------------------------//

IntervalScalarFunction operator-(IntervalScalarFunction const& f1, IntervalScalarFunction const& f2) {
    shared_ptr<IntervalScalarFunctionModelInterface const> f1p=boost::dynamic_pointer_cast<IntervalScalarFunctionModelInterface const>(f1.managed_pointer());
    shared_ptr<IntervalScalarFunctionModelInterface const> f2p=boost::dynamic_pointer_cast<IntervalScalarFunctionModelInterface const>(f2.managed_pointer());
    if(f1p && f2p) {
        return IntervalScalarFunctionModel(*f1p) - IntervalScalarFunctionModel(*f2p);
    }
    return new BinaryFunction<Interval>(SUB,f1,f2);
}

IntervalScalarFunction compose(const IntervalScalarFunction& f, const IntervalVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return compose(f,IntervalVectorFunctionModel(dynamic_cast<IntervalVectorFunctionModelInterface const&>(*g.raw_pointer())));
}

IntervalVectorFunction compose(const IntervalVectorFunction& f, const IntervalVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return compose(f,IntervalVectorFunctionModel(dynamic_cast<IntervalVectorFunctionModelInterface const&>(*g.raw_pointer())));
}

IntervalVectorFunction join(IntervalVectorFunction const& f1, const IntervalVectorFunction& f2) {
    if(dynamic_cast<VectorTaylorFunction const*>(f1.raw_pointer()) && dynamic_cast<VectorTaylorFunction const*>(f2.raw_pointer())) {
        return join(dynamic_cast<VectorTaylorFunction const&>(*f1.raw_pointer()),dynamic_cast<VectorTaylorFunction const&>(*f2.raw_pointer()));
    }
    VectorOfScalarFunction<Interval> r(f1.result_size()+f2.result_size(),f1.argument_size());
    for(uint i=0; i!=f1.result_size(); ++i) { r[i]=f1[i]; }
    for(uint i=0; i!=f2.result_size(); ++i) { r[i+f1.result_size()]=f2[i]; }
    return r;
}









} // namespace Ariadne
