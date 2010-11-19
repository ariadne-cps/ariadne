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
#include "function.h"
#include "function_template.h"
#include "polynomial.h"
#include "differential.h"
#include "taylor_model.h"
#include "operators.h"
#include "variables.h"
#include "expression.h"
#include "assignment.h"
#include "space.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

template<class T> inline std::string str(const T& t) {
    std::stringstream ss; ss << t; return ss.str(); }

typedef uint SizeType;
typedef uint Nat;
typedef int Int;

inline std::ostream& operator<<(std::ostream& os, const ScalarFunctionInterface<Void>* f) { return f->repr(os); }





struct ScalarExpressionFunctionBody
    : ScalarFunctionTemplate<ScalarExpressionFunctionBody>
{
    Space<Real> _space;
    Expression<Real> _expression;

    ScalarExpressionFunctionBody(const Expression<Real>& e, const Space<Real>& s)
        : _space(s), _expression(e) { }
    virtual SizeType argument_size() const { return _space.size(); }
    virtual RealScalarFunction derivative(uint j) const { return RealScalarFunction(Ariadne::derivative(_expression,_space[j]),_space); }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "EF["<<this->argument_size()<<"]( space="<<this->_space<<", expression="<<this->_expression<<" )"; }

    template<class X> void _compute(X& r, const Vector<X>& x) const {
        Map<ExtendedRealVariable,X> valuation;
        for(uint i=0; i!=x.size(); ++i) { valuation[_space.variable(i)]=x[i]; }
        r=Ariadne::evaluate(_expression,valuation); }
};

struct VectorExpressionFunctionBody
    : VectorFunctionTemplate<VectorExpressionFunctionBody>
{
    VectorExpressionFunctionBody(const List<ExtendedRealVariable>& rv,
                             const List< Assignment<ExtendedRealVariable,RealExpression> >& eq,
                             const List<RealVariable>& av)
        : _result_variables(rv), _argument_variables(av), _assignments(eq) { }

    List<ExtendedRealVariable> _result_variables;
    List<RealVariable> _argument_variables;
    List< Assignment<ExtendedRealVariable,RealExpression> > _assignments;

    virtual uint result_size() const { return this->_result_variables.size(); }
    virtual uint argument_size() const { return this->_argument_variables.size(); }

    virtual VectorExpressionFunctionBody* clone() const { return new VectorExpressionFunctionBody(*this); }

    template<class X> void _compute(Vector<X>& result, const Vector<X>& x) const {
        Map<ExtendedRealVariable,X> values;
        for(uint i=0; i!=_argument_variables.size(); ++i) {
            values[_argument_variables[i]]=x[i];
        }
        for(uint i=0; i!=_assignments.size(); ++i) {
            values[_assignments[i].lhs]=Ariadne::evaluate(_assignments[i].rhs,values);
        }
        for(uint i=0; i!=result.size(); ++i) {
            result[i]=values[_result_variables[i]];
        }
    }

    virtual RealScalarFunction operator[](uint i) const {
        if(_result_variables.size()==_assignments.size()) {
            for(uint i=0; i!=_assignments.size(); ++i) {
                if(_assignments[i].lhs==_result_variables[i]) {
                    return RealScalarFunction(_assignments[i].rhs,_argument_variables);
                }
                ARIADNE_THROW(std::runtime_error,"VectorExpressionFunction::operator[]",
                              "Variable "<<_argument_variables[i]<<" has no assignment");
            }
        }
        ARIADNE_FAIL_MSG("Can only compute index of a function with no dependencies");
    }

    std::ostream& write(std::ostream& os) const {
        return os << "EF["<<this->result_size()<<","<<this->argument_size()<<"]"
                  <<"( expression=" << this->_assignments << ")";
    }

};





//! A constant function f(x)=c
struct ScalarConstantFunctionBody
    : ScalarFunctionTemplate<ScalarConstantFunctionBody>
{
    SizeType _argument_size;
    Real _value;

    ScalarConstantFunctionBody(uint as, const Real& c) : _argument_size(as), _value(c) { }
    operator Real() const { return _value; }

    virtual SizeType argument_size() const { return _argument_size; }
    virtual RealScalarFunction derivative(uint j) const { return RealScalarFunction::constant(_argument_size,0.0); }
    virtual std::ostream& repr(std::ostream& os) const { return os << this->_value; }
    virtual std::ostream& write(std::ostream& os) const { return os << "CF["<<this->_argument_size<<"]("<<Float(this->_value)<<")"; }
    template<class X> void _compute(X& r, const Vector<X>& x) const {
        if(x.size()==0) { r=_value; } else { r=x[0]*0+_value; } }
};


//! A coordinate function \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=x_i\f$.
struct ScalarCoordinateFunctionBody
    : ScalarFunctionTemplate<ScalarCoordinateFunctionBody>
{
    SizeType _argument_size;
    SizeType _index;

    ScalarCoordinateFunctionBody(uint as, uint i) : _argument_size(as), _index(i) { }
    SizeType index() const { return _index; }

    virtual SizeType argument_size() const { return _argument_size; }
    virtual RealScalarFunction derivative(uint j) const {
        if(j==_index) { return RealScalarFunction::constant(_argument_size,1.0); }
        else { return RealScalarFunction::constant(_argument_size,0.0); } }
    virtual std::ostream& repr(std::ostream& os) const { return os << "x"<<this->_index; }
    virtual std::ostream& write(std::ostream& os) const { return os << "F["<<this->_argument_size<<"](x"<<this->_index<<")"; }
    template<class X> void _compute(X& r, const Vector<X>& x) const { r=x[_index]; }
};


//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
template<class Op> struct UnaryFunctionBody
    : ScalarFunctionTemplate< UnaryFunctionBody<Op> >
{
  public:
    UnaryFunctionBody(const Op& op, const RealScalarFunction& arg)
        : _op(op), _arg(arg) { }
    virtual UnaryFunctionBody* clone() const { return new UnaryFunctionBody<Op>(*this); }
    virtual SizeType argument_size() const {
        return this->_arg.argument_size(); }

    virtual RealScalarFunction derivative(uint j) const {
        switch(_op.code()) {
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
        return os << _op << "(" << _arg._raw_pointer() << ")"; }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "F["<<this->argument_size()<<"](" << this <<")"; }

    template<class R, class A> void _compute(R& r, const A& x) const { r=_op(_arg.evaluate(x)); }

    Op _op;
    RealScalarFunction _arg;
};



template<class Op> struct BinaryFunctionBody
    : ScalarFunctionTemplate< BinaryFunctionBody<Op> >
{
  public:
    BinaryFunctionBody(const Op& op, const RealScalarFunction& arg1, const RealScalarFunction& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) { ARIADNE_ASSERT_MSG(arg1.argument_size()==arg2.argument_size(),"op='"<<op<<"', arg1="<<arg1<<", arg2="<<arg2); }
    virtual BinaryFunctionBody<Op>* clone() const { return new BinaryFunctionBody<Op>(*this); }
    virtual SizeType argument_size() const {
        return this->_arg1.argument_size(); }

    virtual RealScalarFunction derivative(uint j) const {
        switch(_op.code()) {
            case ADD:
                return _arg1.derivative(j)+_arg2.derivative(j);
            case SUB:
                return _arg1.derivative(j)-_arg2.derivative(j);
            case MUL:
                return _arg1.derivative(j)*_arg2+_arg1*_arg2.derivative(j);
            case DIV:
                if(dynamic_cast<const ScalarConstantFunctionBody*>(_arg2.pointer())) {
                    return _arg1.derivative(j)/_arg2;
                } else {
                    return _arg1.derivative(j)/_arg2-_arg2.derivative(j)*_arg1/sqr(_arg2);
                }
            default: ARIADNE_FAIL_MSG("Unknown binary function "<<this->_op);
        }
    }

    virtual std::ostream& repr(std::ostream& os) const {
        return os << _arg1._raw_pointer() << symbol(_op.code()) << _arg2._raw_pointer(); }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "F["<<this->argument_size()<<"]("<< this <<")"; }

    template<class R, class A> void _compute(R& r, const A& x) const { r=_op(_arg1.evaluate(x),_arg2.evaluate(x)); }

    Op _op;
    RealScalarFunction _arg1;
    RealScalarFunction _arg2;
};


// \brief The power function \f$(x,n)\mapsto x^n\f$.
template<> struct BinaryFunctionBody<Pow>
    : ScalarFunctionTemplate< BinaryFunctionBody<Pow> >
{
  public:
    BinaryFunctionBody(const Pow& op, const RealScalarFunction& arg1, const Int& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {  }
    virtual BinaryFunctionBody<Pow>* clone() const { return new BinaryFunctionBody<Pow>(*this); }
    virtual SizeType argument_size() const {
        return this->_arg1.argument_size(); }

    virtual RealScalarFunction derivative(uint j) const {
        if(_arg2==0) { return RealScalarFunction::constant(this->argument_size(),0.0); }
        if(_arg2==1) { return _arg1.derivative(j); }
        if(_arg2==2) { return 2*_arg1.derivative(j)*_arg1; }
        return _arg2*_arg1.derivative(j)*pow(_arg1,_arg2-1);
    }

    virtual std::ostream& repr(std::ostream& os) const {
        return os << "pow(" << _arg1._raw_pointer() << "," << _arg2 << ")"; }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "F["<<this->argument_size()<<"]("<< this <<")"; }

    template<class R, class A> void _compute(R& r, const A& x) const { r=pow(_arg1.evaluate(x),_arg2); }

    Pow _op;
    RealScalarFunction _arg1;
    Int _arg2;
};






//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
struct ScalarAffineFunctionBody
    : ScalarFunctionTemplate<ScalarAffineFunctionBody>
{
    explicit ScalarAffineFunctionBody(uint n) : _affine(n) { }
    ScalarAffineFunctionBody(const Affine<Real>& aff) : _affine(aff) { }
    ScalarAffineFunctionBody(const Vector<Real>& g, const Real& c) : _affine(g,c) { }

    virtual SizeType argument_size() const { return _affine.argument_size(); }
    virtual RealScalarFunction derivative(uint j) const {
        return RealScalarFunction::constant(_affine.argument_size(),_affine.gradient(j)); }
    virtual std::ostream& repr(std::ostream& os) const { return os << Affine<Float>(_affine); }
    virtual std::ostream& write(std::ostream& os) const { return os << "AF["<<this->argument_size()<<"]("<<Affine<Float>(_affine)<<")"; }
    template<class X> void _compute(X& r, const Vector<X>& x) const {
        r=x[0]*0+_affine.value(); for(uint j=0; j!=argument_size(); ++j) { r+=_affine.gradient(j)*x[j]; } }
    Affine<Real> _affine;
};


struct ScalarPolynomialFunctionBody
    : ScalarFunctionTemplate<ScalarPolynomialFunctionBody>
{
    explicit ScalarPolynomialFunctionBody(uint n) : _polynomial(n) { }
    ScalarPolynomialFunctionBody(const Polynomial<Real>& p) : _polynomial(p) { }

    virtual SizeType argument_size() const { return _polynomial.argument_size(); }
    virtual RealScalarFunction derivative(uint j) const { return RealScalarFunction(Ariadne::derivative(_polynomial,j)); }
    virtual std::ostream& repr(std::ostream& os) const {
        return os << Polynomial<Float>(_polynomial); }
    virtual std::ostream& write(std::ostream& os) const { return
        os << "PF["<<this->argument_size()<<"]("<<Polynomial<Float>(_polynomial)<<")"; }
    template<class X> void _compute(X& r, const Vector<X>& x) const {
        r=Ariadne::evaluate(_polynomial,x); }
    Polynomial<Real> _polynomial;
};

struct ScalarExpansionFunctionBody
: ScalarFunctionTemplate<ScalarExpansionFunctionBody>
{
    explicit ScalarExpansionFunctionBody(uint n) : _expansion(n) { }
    ScalarExpansionFunctionBody(const Expansion<Real>& p) : _expansion(p) { }

    virtual SizeType argument_size() const { return _expansion.argument_size(); }
    virtual RealScalarFunction derivative(uint j) const { return RealScalarFunction(Ariadne::derivative(Polynomial<Real>(_expansion),j)); }
    virtual std::ostream& repr(std::ostream& os) const {
        return os << Expansion<Float>(_expansion); }
    virtual std::ostream& write(std::ostream& os) const { return
        os << "EF["<<this->argument_size()<<"]("<<Polynomial<Float>(_expansion)<<")"; }
    template<class X> void _compute(X& r, const Vector<X>& x) const {
        r=x[0]*0; r=Ariadne::evaluate(_expansion,x); }
        Expansion<Real> _expansion;
};

inline void add_error(Float& x, const Float& e) { }
inline void add_error(Interval& x, const Float& e) { x+=Interval(-e,+e); }
inline void add_error(Differential<Float>& x, const Float& e) { }
inline void add_error(Differential<Interval>& x, const Float& e) { x+=Interval(-e,+e); }
inline void add_error(TaylorModel<Float>& x, const Float& e) { }
inline void add_error(TaylorModel<Interval>& x, const Float& e) { x+=Interval(-e,+e); }

/*

class ScalarTaylorFunction;

//! An Taylor function expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=p(s^{-1}(x)\f$.
struct ScalarTaylorFunctionBody
    : ScalarFunctionTemplate<ScalarTaylorFunctionBody>
{
    explicit ScalarTaylorFunctionBody(uint n) : _domain(n,Interval(-1,+1)), _polynomial(n), _error(0.0) { }
    ScalarTaylorFunctionBody(const ScalarTaylorFunction& tm) : _domain(tm._domain), _polynomial(tm._polynomial), _error(tm._error) { }

    virtual SizeType argument_size() const { return _polynomial.argument_size(); }
    virtual RealScalarFunction derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { return os << _polynomial; }
    virtual std::ostream& repr(std::ostream& os) const { return os << "TF["<<this->argument_size()<<"]("<<_polynomial<<")"; }
    template<class X> void _compute(X& r, const Vector<X>& x) const {
        Vector<X> u=x; for(uint i=0; i!=n; ++i) { u[i]/=_domain[i].radius(); u[i]-=_domain[i].centre(); }
        r=evaluate(_polynomial,u); add_error(r,_error); }

    Vector<Interval> _domain;
    Polynomial<Float> _polynomial;
    Float _error;
};

*/



struct VectorOfScalarFunctionBody
    : VectorFunctionTemplate<VectorOfScalarFunctionBody>
{
    VectorOfScalarFunctionBody(uint rs, uint as)
        : _as(as), _vec(rs,RealScalarFunction(as)) { }
    VectorOfScalarFunctionBody(uint rs, const RealScalarFunction& f)
        : _as(f.argument_size()), _vec(rs,f) { }

    void set(uint i, const RealScalarFunction& f) {
        if(this->argument_size()==0u) { this->_as=f.argument_size(); }
        ARIADNE_ASSERT(f.argument_size()==this->argument_size());
        this->_vec[i]=f; }
    RealScalarFunction get(uint i) const {
        return this->_vec[i]; }

    virtual uint result_size() const {
        return _vec.size(); }
    virtual uint argument_size() const {
        return _as; }

    virtual RealScalarFunction operator[](uint i) const {
        return RealScalarFunction(this->_vec[i]); }

    virtual std::ostream& write(std::ostream& os) const {
        os << "F[" << this->result_size() << "," << this->argument_size() << "][";
        for(uint i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i]._raw_pointer()->repr(os); }
        return os << "]"; }

    template<class X> void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=Vector<X>(this->_vec.size());
        for(uint i=0; i!=r.size(); ++i) {
            r[i]=_vec[i].evaluate(x); } }

    uint _as;
    Vector<RealScalarFunction> _vec;

};





//! A constant function \f$ x\mapsto c\f$ from \f$R^m\f$ to \f$R^n\f$.
struct VectorConstantFunctionBody
    : VectorFunctionTemplate<VectorConstantFunctionBody>
{
    VectorConstantFunctionBody(uint rs, double c0, ...) : _as(), _c(rs) {
        ARIADNE_ASSERT(rs>0); va_list args; va_start(args,c0);
        _c[0]=c0; for(size_t i=1; i!=rs; ++i) { _c[i]=va_arg(args,double); } _as=va_arg(args,int);
        va_end(args);
    }

    VectorConstantFunctionBody(uint rs, Interval c0, ...) : _as(), _c(rs) {
        ARIADNE_ASSERT(rs>0); double l,u; va_list args; va_start(args,c0);
        _c[0]=c0; for(size_t i=1; i!=rs; ++i) { l=va_arg(args,double); u=va_arg(args,double); _c[i]=Interval(l,u); } _as=va_arg(args,int);
        va_end(args);
    }

    VectorConstantFunctionBody(const Vector<Float>& c, uint as) : _as(as), _c(c) { }
    VectorConstantFunctionBody(const Vector<Interval>& c, uint as) : _as(as), _c(c) { }
    VectorConstantFunctionBody(const Vector<Real>& c, uint as) : _as(as), _c(c) { }

    const Vector<Real>& c() const { return _c; }

    SizeType result_size() const { return _c.size(); }
    SizeType argument_size() const { return _as; }
    RealScalarFunction operator[](uint i) const { return RealScalarFunction::constant(this->_as,this->_c[i]); }
    std::ostream& write(std::ostream& os) const {
        return os << "VectorConstantFunctionBody( argument_size=" << this->argument_size()
                  << ", c=" << this->c() << " )"; }

    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_c[i]; } }

    uint _as;
    Vector<Real> _c;
};


// An affine function \f$x\mapsto Ax+b\f$.
struct VectorAffineFunctionBody
    : public VectorFunctionTemplate<VectorAffineFunctionBody>
{

    VectorAffineFunctionBody(const Matrix<Real>& A, const Vector<Real>& b)
        : _A(A), _b(b) { ARIADNE_ASSERT(A.row_size()==b.size()); }
    VectorAffineFunctionBody(uint rs, double b0, ...) : _A(), _b(rs) {
        ARIADNE_ASSERT(rs>0);
        va_list args; va_start(args,b0);
        _b[0]=b0; for(size_t i=1; i!=rs; ++i) { _b[i]=va_arg(args,double); }
        uint as=va_arg(args,int); ARIADNE_ASSERT(as>0); _A=Matrix<Real>(rs,as);
        for(size_t i=0; i!=rs; ++i) { for(size_t j=0; i!=as; ++j) { _A[i][j]=va_arg(args,double); } }
        va_end(args);
    }

    const Matrix<Real>& A() const { return _A; }
    const Vector<Real>& b() const { return _b; }

    virtual SizeType result_size() const { return _A.row_size(); }
    virtual SizeType argument_size() const { return _A.column_size(); }

    template<class X> void _compute(Vector<X>& r, const Vector<X>& x) const { r=prod(_A,x)+_b; }

    virtual Matrix<Real> jacobian(const Vector<Real>& x) const { return this->_A; }

    virtual RealScalarFunction operator[](unsigned int i) const  {
        ARIADNE_ASSERT(i<this->result_size());
        Vector<Real> g(this->argument_size()); for(unsigned int j=0; j!=g.size(); ++j) { g[j]=this->_A[i][j]; }
        return ScalarAffineFunction(g,_b[i]); }

    virtual std::ostream& write(std::ostream& os) const;

    Matrix<Real> _A; Vector<Real> _b;
};

inline std::ostream& VectorAffineFunctionBody::write(std::ostream& os) const {
    const VectorAffineFunctionBody& f=*this;
    os << "AF[" << f.result_size() <<"," << f.argument_size() << "]";
    for(uint i=0; i!=f.result_size(); ++i) {
        os<<(i==0?"( ":"; ");
        if(f.b()[i]!=0) { os<<Float(f._b[i]); }
        for(uint j=0; j!=f.argument_size(); ++j) {
            if(f._A[i][j]!=0) {
                if(f._A[i][j]>0) { os<<"+"; } else { os<<"-"; }
                if(abs(f._A[i][j])!=1) { os<<abs(Float(f._A[i][j])); }
                //ss<<char('x'+j);
                os<<"x"<<j;
            }
        }
    }
    os<<" )";
    return os;
}



//! A polynomial function.
struct VectorPolynomialFunctionBody
{
    explicit VectorPolynomialFunctionBody(uint rs, uint as) : _polynomials(rs,Polynomial<Real>(as)) { }
    VectorPolynomialFunctionBody(const Vector< Polynomial<Real> >& p) : _polynomials(p) { }

    virtual SizeType result_size() const { return _polynomials.size(); }
    virtual SizeType argument_size() const { return _polynomials[0].argument_size(); }
    virtual RealScalarFunction operator[](uint i) const { return RealScalarFunction(this->_polynomials[i]); }
    virtual std::ostream& write(std::ostream& os) const { return os << _polynomials; }
    template<class R, class A> void _compute(R& r, const A& x) const { r=Ariadne::evaluate(_polynomials,x); }

    Vector< Polynomial<Real> > _polynomials;
};




// \brief A projection function \f$ x'_i= x_{p(i)}\f$.
struct ProjectionFunctionBody
    : VectorFunctionTemplate<ProjectionFunctionBody>
{

    // \brief Construct the identity function in dimension \a n.
    ProjectionFunctionBody(uint n) : _as(n), _p(n) {
        for(uint i=0; i!=n; ++i) { _p[i]=i; } }
    // \brief Construct the projection functions \f$f_i(x)=x_{i+k}\f$ for \f$i=0,\ldots,m-1\f$. Precondition: \f$m+k\leq n\f$.
    ProjectionFunctionBody(uint m, uint n, uint k) : _as(n), _p(m) {
        ARIADNE_ASSERT(m+k<=n); for(uint j=0; j!=m; ++j) { _p[j]=k+j; } }
    // \brief Construct the projection function  with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,m-1\f$.
    ProjectionFunctionBody(uint m, uint n, const array<uint>& p) : _as(n), _p(p) {
        ARIADNE_ASSERT(p.size()==m); for(uint i=0; i!=_p.size(); ++i) { ARIADNE_ASSERT(p[i]<n); } }
    // \brief Construct the projection function with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,|p|-1\f$.
    ProjectionFunctionBody(const array<uint>& p, uint n) : _as(n), _p(p) {
        for(uint i=0; i!=_p.size(); ++i) { ARIADNE_ASSERT(p[i]<n); } }
    ProjectionFunctionBody(const Range& rng, uint as) : _as(as), _p(rng.size()) {
        ARIADNE_ASSERT(rng.start()+rng.size()<=as);
        for(uint i=0; i!=_p.size(); ++i) { _p[i]=rng.start()+i; } }
    SizeType p(unsigned int j) const { return this->_p[j]; }
    SizeType result_size() const { return this->_p.size(); }
    SizeType argument_size() const { return this->_as; }
    RealScalarFunction operator[](uint i) const { return RealScalarFunction::coordinate(_as,_p[i]); }
    std::ostream& write(std::ostream& os) const {
        return os << "ProjectionFunctionBody( p=" << this->_p << " )"; }

    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=x[_p[i]]; } }

    uint _as;
    array<uint> _p;
};




//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
struct IdentityFunctionBody
    : VectorFunctionTemplate<IdentityFunctionBody>
{
    IdentityFunctionBody(uint n) : _n(n) { }
    SizeType result_size() const { return this->_n; }
    SizeType argument_size() const { return this->_n; }
    RealScalarFunction operator[](uint i) const { return RealScalarFunction::coordinate(_n,i); }
    std::ostream& write(std::ostream& os) const {
        return os << "IdentityFunctionBody( size=" << this->result_size() << " )"; }

    template<class R, class A> void _compute(R& r, const A& x) const { r=x; }

    uint _n;
};

//! \brief The unscaling function \f$ x_i\mapsto (x_i-c_i)/r_i\f$ in \f$\R^n\f$.
struct VectorUnscalingFunctionBody
: VectorFunctionTemplate<VectorUnscalingFunctionBody>
{
    VectorUnscalingFunctionBody(const Vector<Interval>& bx) : _bx(bx) { }
    SizeType result_size() const { return this->_bx.size(); }
    SizeType argument_size() const { return this->_bx.size(); }
    RealScalarFunction operator[](uint i) const { return (RealScalarFunction::coordinate(_bx.size(),i)-Real(med_ivl(_bx[i])))/Real(rad_ivl(_bx[i])); }
    std::ostream& write(std::ostream& os) const {
        return os << "VectorUnscalingFunction( domain=" << this->_bx << " )"; }

    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) {
            r[i]=(x[i]-Real(med_ivl(_bx[i])))/Real(rad_ivl(_bx[i])); } }

    Vector<Interval> _bx;
};



// A Lie deriviative \f$\nabla g\cdot f\f$.
struct LieDerivativeFunctionBody
    : ScalarFunctionTemplate<LieDerivativeFunctionBody>
{
    //! \brief Construct the identity function in dimension \a n.
    LieDerivativeFunctionBody(const RealScalarFunction& g, const RealVectorFunction& f) {
        ARIADNE_ASSERT(g.argument_size()==f.argument_size());
        ARIADNE_ASSERT(f.result_size()==f.argument_size());
        _g=g; for(uint j=0; j!=g.argument_size(); ++j) { _dg[j]=g.derivative(j); } _f=f; }
    SizeType argument_size() const { return _g.argument_size(); }
    std::ostream& write(std::ostream& os) const { return os << "LieDerivative( g="<<_g<<", f="<<_f<<" )"; }

    template<class R, class A> void _compute(R& r, const A& x) const {
        Vector<R> fx=_f.evaluate(x); r=0; for(uint i=0; i!=_dg.size(); ++i) { r+=fx[i]+_dg[i].evaluate(x); } }

    RealScalarFunction _g;
    List<RealScalarFunction> _dg;
    RealVectorFunction _f;
};







struct FunctionElementBody
    : ScalarFunctionTemplate<FunctionElementBody>
{
    FunctionElementBody(const RealVectorFunction& f, SizeType i)
        : _f(f), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }

    virtual SizeType argument_size() const { return _f.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os<<_f<<"["<<_i<<"]"; }

    template<class Y> inline void _compute(Y& r, const Vector<Y>& x) const {
        r=this->_f.evaluate(x)[_i]; }

    RealVectorFunction _f;
    SizeType _i;
};





struct ScalarComposedFunctionBody
    : ScalarFunctionTemplate<ScalarComposedFunctionBody>
{
    ScalarComposedFunctionBody(const RealScalarFunction& f, const RealVectorFunction& g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual SizeType argument_size() const { return _g.argument_size(); }
    virtual RealScalarFunction derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { return os << "ScalarComposedFunctionBody( f="<<_f<<", g="<<_g<<" )"; }

    template<class X> inline void _compute(X& r, const Vector<X>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    RealScalarFunction _f;
    RealVectorFunction _g;
};


struct VectorComposedFunctionBody
    :  VectorFunctionTemplate<VectorComposedFunctionBody>
{
    VectorComposedFunctionBody(RealVectorFunction f, RealVectorFunction g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual SizeType result_size() const { return _f.result_size(); }
    virtual SizeType argument_size() const { return _g.argument_size(); }
    virtual RealScalarFunction operator[](uint i) const { return compose(_f[i],_g); }
    virtual std::ostream& write(std::ostream& os) const { return os << "ComposedFunctionBody( f="<<_f<<", g="<<_g<<" )"; }

    template<class X> inline void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    RealVectorFunction _f;
    RealVectorFunction _g;
};

struct JoinedFunctionBody
    : VectorFunctionTemplate<JoinedFunctionBody>
{
    JoinedFunctionBody(RealVectorFunction f1, RealVectorFunction f2)
        : _f1(f1), _f2(f2) { ARIADNE_ASSERT(f1.argument_size()==f2.argument_size()); }
    virtual SizeType result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual SizeType argument_size() const { return _f1.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os << "JoinedFunctionBody( f1="<<_f1<<", f2="<<_f2<<" )"; }

    template<class X> inline void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=join(_f1.evaluate(x),_f2.evaluate(x)); }

    RealVectorFunction _f1;
    RealVectorFunction _f2;
};

class CombinedFunctionBody
    : VectorFunctionTemplate<JoinedFunctionBody>
{
    CombinedFunctionBody(RealVectorFunction f1, RealVectorFunction f2)
        : _f1(f1), _f2(f2) { }
    virtual SizeType result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual SizeType argument_size() const { return _f1.argument_size()+_f2.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os << "CombinedFunctionBody( f1="<<_f1<<", f2="<<_f2<<" )"; }

    template<class X> inline void _compute(Vector<X>& r, const Vector<X>& x) const {
        return r=combine(_f1.evaluate(project(x,range(0,_f1.argument_size()))),
                         _f2.evaluate(project(x,range(_f1.argument_size(),this->argument_size())))); }

    RealVectorFunction _f1;
    RealVectorFunction _f2;
};






//------------------------ RealScalarFunction -----------------------------------//

RealScalarFunction RealScalarFunction::constant(Nat n, double c)
{
    return RealScalarFunction(new ScalarConstantFunctionBody(n,Real(c)));

    Polynomial<Interval> p(n);
    p[MultiIndex::zero(n)]=c;
    return RealScalarFunction(p);
}

RealScalarFunction RealScalarFunction::constant(Nat n, Real c)
{
    return RealScalarFunction(new ScalarConstantFunctionBody(n,c));

    Polynomial<Interval> p(n);
    p[MultiIndex::zero(n)]=Interval(c);
    return RealScalarFunction(p);
}

RealScalarFunction RealScalarFunction::variable(Nat n, Nat j)
{
    ARIADNE_DEPRECATED("RealScalarFunction::variable","Use RealScalarFunction::coordinate instead");
    //return RealScalarFunction(new ScalarCoordinateFunctionBody(n,j));
    Polynomial<Interval> p(n);
    p[MultiIndex::unit(n,j)]=1.0;
    return RealScalarFunction(p);
}

RealScalarFunction RealScalarFunction::coordinate(Nat n, Nat j)
{
    return RealScalarFunction(new ScalarCoordinateFunctionBody(n,j));

    Polynomial<Interval> p(n);
    p[MultiIndex::unit(n,j)]=1.0;
    return RealScalarFunction(p);
}

List<RealScalarFunction> RealScalarFunction::coordinates(Nat n)
{
    List<RealScalarFunction> result(n);
    for(Nat i=0; i!=n; ++i) {
        result[i]=RealScalarFunction::coordinate(n,i);
    }
    return result;
}

RealScalarFunction::ScalarFunction(Nat n)
    : _ptr(new ScalarPolynomialFunctionBody(Polynomial<Interval>::constant(n,Interval(0.0))))
{
}

RealScalarFunction::ScalarFunction(const Polynomial<Real>& p)
    : _ptr(new ScalarPolynomialFunctionBody(Polynomial<Interval>(p)))
{
}

RealScalarFunction::ScalarFunction(const Expression<Real>& e, const Space<Real>& s)
    : _ptr(new ScalarExpressionFunctionBody(e,s))
{
}


RealScalarFunction::ScalarFunction(const Expression<Real>& e, const List< Variable<Real> >& s)
    : _ptr(new ScalarExpressionFunctionBody(e,s))
{
}

RealScalarFunction::ScalarFunction(const Expression<tribool>& e, const List< Variable<Real> >& s)
{
    *this = RealScalarFunction(new ScalarExpressionFunctionBody(indicator(e),s));
}

RealScalarFunction::ScalarFunction(const Expansion<Float>& e)
{
    *this = RealScalarFunction(new ScalarExpansionFunctionBody(e));
}




RealScalarFunction RealScalarFunction::derivative(Nat j) const
{
    return this->pointer()->derivative(j);
}


Polynomial<Real> RealScalarFunction::polynomial() const
{
    const ScalarPolynomialFunctionBody* p=dynamic_cast<const ScalarPolynomialFunctionBody*>(this->pointer());
    if(p) { return p->_polynomial; }
    const ScalarExpressionFunctionBody* e=dynamic_cast<const ScalarExpressionFunctionBody*>(this->pointer());
    if(e) { return Ariadne::polynomial<Real>(e->_expression,e->_space); }
    ARIADNE_THROW(std::runtime_error,"RealScalarFunction::polynomial()","FunctionBody "<<*this<<" is not a polynomial.");
}


template<class Op> inline
RealScalarFunction make_unary_function(Op op, const RealScalarFunction& f) {
    return RealScalarFunction(new UnaryFunctionBody<Op>(op,f)); }

template<class Op> inline
RealScalarFunction make_binary_function(Op op, const RealScalarFunction& f1, const RealScalarFunction& f2) {
    return RealScalarFunction(new BinaryFunctionBody<Op>(op,f1,f2)); }

inline
RealScalarFunction make_binary_function(Pow op, const RealScalarFunction& f1, const Int& n2) {
    return RealScalarFunction(new BinaryFunctionBody<Pow>(op,f1,n2)); }


RealScalarFunction operator+(const RealScalarFunction& f)
{
    //std::cerr<<"+f1, f1="<<f1<<"\n";
    const ScalarPolynomialFunctionBody* p=dynamic_cast<const ScalarPolynomialFunctionBody*>(f.pointer());
    if(p) { return RealScalarFunction( +p->_polynomial ); }
    return f;
}

RealScalarFunction operator-(const RealScalarFunction& f)
{
    //std::cerr<<"-f1, f1="<<f1<<"\n";
    const ScalarPolynomialFunctionBody* p=dynamic_cast<const ScalarPolynomialFunctionBody*>(f.pointer());
    if(p) { return RealScalarFunction( -p->_polynomial ); }
    const ScalarExpressionFunctionBody* e=dynamic_cast<const ScalarExpressionFunctionBody*>(f.pointer());
    if(e) {
        return RealScalarFunction(-e->_expression,e->_space);
    }
    return make_unary_function(Neg(),f);
}

RealScalarFunction operator+(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    //std::cerr<<"f1+f2, f1="<<f1<<", f2="<<f2<<"\n";
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    const ScalarPolynomialFunctionBody* p2=dynamic_cast<const ScalarPolynomialFunctionBody*>(f2.pointer());
    if(p1 && p2) {
        return RealScalarFunction(p1->_polynomial + p2->_polynomial );
    }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    const ScalarExpressionFunctionBody* e2=dynamic_cast<const ScalarExpressionFunctionBody*>(f2.pointer());
    if(e1 && e2 && e1->_space==e2->_space) {
        return RealScalarFunction(e1->_expression+e2->_expression,e1->_space);
    }
    return make_binary_function(Add(),f1,f2);
}

RealScalarFunction operator-(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    //std::cerr<<"f1-f2, f1="<<f1<<", f2="<<f2<<"\n";
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    const ScalarPolynomialFunctionBody* p2=dynamic_cast<const ScalarPolynomialFunctionBody*>(f2.pointer());
    if(p1 && p2) {
        return RealScalarFunction(p1->_polynomial - p2->_polynomial );
    }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    const ScalarExpressionFunctionBody* e2=dynamic_cast<const ScalarExpressionFunctionBody*>(f2.pointer());
    if(e1 && e2 && e1->_space==e2->_space) {
        return RealScalarFunction(e1->_expression-e2->_expression,e1->_space);
    }
    return make_binary_function(Sub(),f1,f2);
}

RealScalarFunction operator*(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    //std::cerr<<"f1*f2, f1="<<f1<<", f2="<<f2<<"\n";
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    const ScalarPolynomialFunctionBody* p2=dynamic_cast<const ScalarPolynomialFunctionBody*>(f2.pointer());
    if(p1 && p2) {
        return RealScalarFunction(p1->_polynomial * p2->_polynomial );
    }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    const ScalarExpressionFunctionBody* e2=dynamic_cast<const ScalarExpressionFunctionBody*>(f2.pointer());
    if(e1 && e2 && e1->_space==e2->_space) {
        return RealScalarFunction(e1->_expression*e2->_expression,e1->_space);
    }
    return make_binary_function(Mul(),f1,f2);
}

RealScalarFunction operator/(const RealScalarFunction& f1, const RealScalarFunction& f2)
{
    //std::cerr<<"f1/f2, f1="<<f1<<", f2="<<f2<<"\n";
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    const ScalarPolynomialFunctionBody* p2=dynamic_cast<const ScalarPolynomialFunctionBody*>(f2.pointer());
    if(p1 && p2) {
        if(p2->_polynomial.degree()==0) {
            return RealScalarFunction(p1->_polynomial/p2->_polynomial);
        }
    }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    const ScalarExpressionFunctionBody* e2=dynamic_cast<const ScalarExpressionFunctionBody*>(f2.pointer());
    if(e1 && e2 && e1->_space==e2->_space) {
        return RealScalarFunction(e1->_expression/e2->_expression,e1->_space);
    }
    return make_binary_function(Div(),f1,f2);
}


RealScalarFunction operator+(const RealScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    if(p1) { return RealScalarFunction(p1->_polynomial + static_cast<Real>(s2) ); }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    if(e1) { return RealScalarFunction(e1->_expression+s2,e1->_space); }
    return f1+RealScalarFunction::constant(f1.argument_size(),s2);
}

RealScalarFunction operator-(const RealScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    if(p1) { return RealScalarFunction(p1->_polynomial - static_cast<Real>(s2) ); }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    if(e1) { return RealScalarFunction(e1->_expression-s2,e1->_space); }
    return f1-RealScalarFunction::constant(f1.argument_size(),s2);
}

RealScalarFunction operator*(const RealScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    if(p1) { return RealScalarFunction(p1->_polynomial * static_cast<Real>(s2) ); }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    if(e1) { return RealScalarFunction(e1->_expression*s2,e1->_space); }
    return f1*RealScalarFunction::constant(f1.argument_size(),s2);
}

RealScalarFunction operator/(const RealScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunctionBody* p1=dynamic_cast<const ScalarPolynomialFunctionBody*>(f1.pointer());
    if(p1) { return RealScalarFunction(p1->_polynomial / static_cast<Real>(s2) ); }
    const ScalarExpressionFunctionBody* e1=dynamic_cast<const ScalarExpressionFunctionBody*>(f1.pointer());
    if(e1) { return RealScalarFunction(e1->_expression/s2,e1->_space); }
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
    const ScalarPolynomialFunctionBody* p2=dynamic_cast<const ScalarPolynomialFunctionBody*>(f2.pointer());
    if(p2) { return RealScalarFunction( static_cast<Real>(s1) - p2->_polynomial ); }
    const ScalarExpressionFunctionBody* e2=dynamic_cast<const ScalarExpressionFunctionBody*>(f2.pointer());
    if(e2) { return RealScalarFunction(s1-e2->_expression,e2->_space); }
    return RealScalarFunction::constant(f2.argument_size(),s1)-f2;
}

RealScalarFunction operator*(const Real& s1, const RealScalarFunction& f2)
{
    return f2*s1;
}

RealScalarFunction operator/(const Real& s1, const RealScalarFunction& f2)
{
    const ScalarExpressionFunctionBody* e2=dynamic_cast<const ScalarExpressionFunctionBody*>(f2.pointer());
    if(e2) { return RealScalarFunction(s1/e2->_expression,e2->_space); }
    return RealScalarFunction::constant(f2.argument_size(),s1)/f2;
}

RealScalarFunction pow(const RealScalarFunction& f, Nat m)
{
    const ScalarPolynomialFunctionBody* p=dynamic_cast<const ScalarPolynomialFunctionBody*>(f.pointer());
    if(p) { return RealScalarFunction(pow(p->_polynomial,m)); }
    const ScalarExpressionFunctionBody* e=dynamic_cast<const ScalarExpressionFunctionBody*>(f.pointer());
    if(e) { return RealScalarFunction(pow(e->_expression,m),e->_space); }
    return make_binary_function(Pow(),f,m);
}


RealScalarFunction pow(const RealScalarFunction& f, Int n)
{
    const ScalarPolynomialFunctionBody* p=dynamic_cast<const ScalarPolynomialFunctionBody*>(f.pointer());
    if(p && n>=0) { return RealScalarFunction(pow(p->_polynomial,Nat(n))); }
    const ScalarExpressionFunctionBody* e=dynamic_cast<const ScalarExpressionFunctionBody*>(f.pointer());
    if(e) { return RealScalarFunction(pow(e->_expression,n),e->_space); }
    return make_binary_function(Pow(),f,n);
}

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



RealScalarFunction embed(const RealScalarFunction& f, uint k) {
    const ScalarPolynomialFunctionBody* p=dynamic_cast<const ScalarPolynomialFunctionBody*>(f.pointer());
    if(p) {
        return RealScalarFunction(new ScalarPolynomialFunctionBody(embed(0u,p->_polynomial,k)));
    }
    const ScalarExpressionFunctionBody* e=dynamic_cast<const ScalarExpressionFunctionBody*>(f.pointer());
    if(e) {
        Space<Real> new_space=e->_space;
        for(uint i=0; i!=k; ++i) { new_space.append(Variable<Real>("x"+str(k))); }
        return new ScalarExpressionFunctionBody(e->_expression,new_space);
    }
    ARIADNE_FAIL_MSG("Cannot embed function "<<f<<" in a different space.");
}




//------------------------ Vector Function ----------------------------------//

RealVectorFunction::VectorFunction(VectorFunctionInterface<Real>* fptr)
    : _ptr(fptr)
{
}

RealVectorFunction::VectorFunction()
    : _ptr(new VectorOfScalarFunctionBody(0u,RealScalarFunction()))
{
}

RealVectorFunction::VectorFunction(Nat rs, Nat as)
    : _ptr(new VectorOfScalarFunctionBody(rs,RealScalarFunction::constant(as,Real(0.0))))
{
}

RealVectorFunction::VectorFunction(Nat rs, const RealScalarFunction& sf)
    : _ptr(new VectorOfScalarFunctionBody(rs,sf))
{
}

RealVectorFunction::VectorFunction(const List<RealScalarFunction>& lsf)
{
    ARIADNE_ASSERT(lsf.size()>0);
    this->_ptr=shared_ptr<RealVectorFunctionInterface>(new VectorOfScalarFunctionBody(lsf.size(),lsf[0].argument_size()));
    VectorOfScalarFunctionBody* vec = static_cast<VectorOfScalarFunctionBody*>(this->_ptr.operator->());
    for(uint i=0; i!=lsf.size(); ++i) {
        vec->set(i,lsf[i]);
    }
}

RealVectorFunction::VectorFunction(const Vector< Polynomial<Real> >& p)
{
    ARIADNE_ASSERT(p.size()>0);
    this->_ptr=shared_ptr<RealVectorFunctionInterface>(new VectorOfScalarFunctionBody(p.size(),p[0].argument_size()));
    VectorOfScalarFunctionBody* vec = static_cast<VectorOfScalarFunctionBody*>(this->_ptr.operator->());
    for(uint i=0; i!=p.size(); ++i) {
        vec->set(i,RealScalarFunction(p[i]));
    }
}

RealVectorFunction::VectorFunction(const List< Expression<Real> >& e, const List< Variable<Real> >& s)
    : _ptr(new VectorOfScalarFunctionBody(e.size(),s.size()))
{
    VectorOfScalarFunctionBody* vec = static_cast<VectorOfScalarFunctionBody*>(this->_ptr.operator->());
    for(uint i=0; i!=e.size(); ++i) {
        vec->set(i,RealScalarFunction(e[i],s));
    }
}

/*
RealVectorFunction::RealVectorFunction(const Space<Real>& rs, const Map<RealVariable,RealExpression>& e, const Space<Real>& as)
    : _ptr(new VectorOfScalarFunctionBody(rs.size(),as.size()))
{
    VectorOfScalarFunctionBody* vec = static_cast<VectorOfScalarFunctionBody*>(this->_ptr.operator->());
    for(uint i=0; i!=rs.size(); ++i) {
        vec->set(i,RealScalarFunction(e[rs[i]],as));
    }
}
*/






RealVectorFunction::VectorFunction(const List<ExtendedRealVariable>& rv,
                               const List<ExtendedRealAssignment>& eq,
                               const List<RealVariable>& av)
    : _ptr(new VectorExpressionFunctionBody(rv,eq,av))
{
}

RealVectorFunction::VectorFunction(const List<RealVariable>& rv,
                               const List<RealAssignment>& eq,
                               const List<RealVariable>& av)
    : _ptr(new VectorExpressionFunctionBody(List<ExtendedRealVariable>(rv),List<ExtendedRealAssignment>(eq),av))
{
}

RealVectorFunction::VectorFunction(const List<DottedRealVariable>& rv,
                               const List<DottedRealAssignment>& eq,
                               const List<RealVariable>& av)
    : _ptr(new VectorExpressionFunctionBody(List<ExtendedRealVariable>(rv),List<ExtendedRealAssignment>(eq),av))
{
}

RealVectorFunction::VectorFunction(const List<PrimedRealVariable>& rv,
                               const List<PrimedRealAssignment>& eq,
                               const List<RealVariable>& av)
    : _ptr(new VectorExpressionFunctionBody(List<ExtendedRealVariable>(rv),List<ExtendedRealAssignment>(eq),av))
{
}




RealVectorFunction RealVectorFunction::constant(const Vector<Real>& c, Nat as)
{
    const Nat rs=c.size();
    VectorOfScalarFunctionBody* res = new VectorOfScalarFunctionBody(rs,as);
    for(uint i=0; i!=rs; ++i) {
        res->_vec[i]=RealScalarFunction::constant(as,c[i]);
    }
    return RealVectorFunction(res);
}

RealVectorFunction RealVectorFunction::identity(Nat n)
{
    VectorOfScalarFunctionBody* res = new VectorOfScalarFunctionBody(n,n);
    for(uint i=0; i!=n; ++i) {
        res->_vec[i]=RealScalarFunction::coordinate(n,i);
    }
    return RealVectorFunction(res);
}

RealScalarFunction RealVectorFunction::get(Nat i) const
{
    return this->_ptr->operator[](i);
}

void RealVectorFunction::set(Nat i, RealScalarFunction f)
{
    VectorOfScalarFunctionBody* vptr=dynamic_cast<VectorOfScalarFunctionBody*>(this->_ptr.operator->());
    ARIADNE_ASSERT_MSG(vptr,"Can only set component of a vector of scalar functions.");
    vptr->set(i,f);
}

RealScalarFunction RealVectorFunction::operator[](Nat i) const
{
    return this->get(i);
}

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




//------------------------ Function operators -------------------------------//

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

RealScalarFunction compose(const RealScalarFunction& f, const RealVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return RealScalarFunction(new ScalarComposedFunctionBody(f,g));
}

RealVectorFunction compose(const RealVectorFunction& f, const RealVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return RealVectorFunction(new VectorComposedFunctionBody(f,g));
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



//------------------------ Special functions --------------------------------//

ScalarAffineFunction::ScalarAffineFunction(const Vector<Real>& a, const Real& b)
    : RealScalarFunction(new ScalarAffineFunctionBody(a,b))
{
}

void ScalarAffineFunction::_check_type(const RealScalarFunctionInterface* pointer) const {
    ARIADNE_ASSERT(dynamic_cast<const ScalarAffineFunctionBody*>(pointer)); }



VectorConstantFunction::VectorConstantFunction(const Vector<Real>& c, uint as)
    : RealVectorFunction(new VectorConstantFunctionBody(c,as))
{
}

void VectorConstantFunction::_check_type(const RealVectorFunctionInterface* pointer) const { }


VectorAffineFunction::VectorAffineFunction(const Matrix<Real>& A, const Vector<Real>& b)
    : RealVectorFunction(new VectorAffineFunctionBody(A,b))
{
}

void VectorAffineFunction::_check_type(const RealVectorFunctionInterface* pointer) const {
    ARIADNE_ASSERT(dynamic_cast<const VectorAffineFunctionBody*>(pointer)); }

const Matrix<Real> VectorAffineFunction::A() const
{
    return static_cast<const VectorAffineFunctionBody*>(this->pointer())->A();
}

const Vector<Real> VectorAffineFunction::b() const
{
    return static_cast<const VectorAffineFunctionBody*>(this->pointer())->b();
}



IdentityFunction::IdentityFunction(uint n)
    : RealVectorFunction(new IdentityFunctionBody(n))
{
}

void IdentityFunction::_check_type(const RealVectorFunctionInterface* pointer) const {
    ARIADNE_ASSERT(dynamic_cast<const IdentityFunctionBody*>(pointer)); }




ProjectionFunction::ProjectionFunction(uint m, uint n, uint k)
    : RealVectorFunction(new ProjectionFunctionBody(m,n,k))
{
}

ProjectionFunction::ProjectionFunction(const array<uint>& p, uint n)
    : RealVectorFunction(new ProjectionFunctionBody(p,n))
{
}

ProjectionFunction::ProjectionFunction(uint m, uint n, const array<uint>& p)
    : RealVectorFunction(new ProjectionFunctionBody(m,n,p))
{
}

const uint ProjectionFunction::p(uint i) const
{
    return static_cast<const ProjectionFunctionBody*>(this->pointer())->p(i);
}

void ProjectionFunction::_check_type(const RealVectorFunctionInterface* pointer) const {
    ARIADNE_ASSERT(dynamic_cast<const ProjectionFunctionBody*>(pointer)); }


VectorUnscalingFunction::VectorUnscalingFunction(const Vector<Interval>& bx)
    : RealVectorFunction(new VectorUnscalingFunctionBody(bx))
{
}

void VectorUnscalingFunction::_check_type(const RealVectorFunctionInterface* pointer) const {
    ARIADNE_ASSERT(dynamic_cast<const VectorUnscalingFunction*>(pointer)); }




} // namespace Ariadne
