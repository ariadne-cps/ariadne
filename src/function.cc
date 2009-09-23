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

#include "function.h"
#include "real.h"
#include "polynomial.h"

namespace Ariadne {

typedef uint Nat;

ScalarFunction ScalarFunction::constant(Nat n, Real c)
{
    Polynomial<Interval> p(n);
    p[MultiIndex::zero(n)]=c;
    return ScalarFunction(p);
}

ScalarFunction ScalarFunction::variable(Nat n, Nat j)
{
    Polynomial<Interval> p(n);
    p[MultiIndex::unit(n,j)]=1.0;
    return ScalarFunction(p);
}

ScalarFunction::ScalarFunction(Nat n)
    : _ptr(new ScalarPolynomialFunction(Polynomial<Interval>(n)))
{
}

ScalarFunction::ScalarFunction(const Polynomial<Real>& p)
    : _ptr(new ScalarPolynomialFunction(Polynomial<Interval>(p)))
{
}

ScalarFunction::ScalarFunction(const Expression<Real>& e, const Space<Real>& s)
    : _ptr(new ScalarExpressionFunction(e,s))
{
}


Vector<Float> ScalarFunction::gradient(const Vector<Float>& x) const {
    return this->_ptr->evaluate(Differential<Float>::variables(1u,x)).gradient();
}

Vector<Interval> ScalarFunction::gradient(const Vector<Interval>& x) const
{
    return this->_ptr->evaluate(Differential<Interval>::variables(1u,x)).gradient();
}


ScalarFunction ScalarFunction::derivative(Nat j) const
{
    const ScalarPolynomialFunction* p=dynamic_cast<const ScalarPolynomialFunction*>(this->pointer());
    if(p) { return ScalarFunction(Ariadne::derivative(static_cast<const Polynomial<Interval>&>(*p),j)); }
    ARIADNE_THROW(std::runtime_error,"ScalarFunction::derivative()","Function "<<*this<<" is not a polynomial.");
}


Polynomial<Real> ScalarFunction::polynomial() const
{
    const ScalarPolynomialFunction* p=dynamic_cast<const ScalarPolynomialFunction*>(this->pointer());
    if(p) { return Polynomial<Real>(static_cast<const Polynomial<Interval>&>(*p)); }
    ARIADNE_THROW(std::runtime_error,"ScalarFunction::polynomial()","Function "<<*this<<" is not a polynomial.");
}


ScalarFunction operator+(const ScalarFunction& f)
{
    const ScalarPolynomialFunction* p=dynamic_cast<const ScalarPolynomialFunction*>(f.pointer());
    if(p) { return ScalarFunction( +static_cast<const Polynomial<Interval>&>(*p) ); }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator-(const ScalarFunction& f)
{
    const ScalarPolynomialFunction* p=dynamic_cast<const ScalarPolynomialFunction*>(f.pointer());
    if(p) { return ScalarFunction( -static_cast<const Polynomial<Interval>&>(*p) ); }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator+(const ScalarFunction& f1, const ScalarFunction& f2)
{
    const ScalarPolynomialFunction* p1=dynamic_cast<const ScalarPolynomialFunction*>(f1.pointer());
    const ScalarPolynomialFunction* p2=dynamic_cast<const ScalarPolynomialFunction*>(f2.pointer());
    if(p1 && p2) {
        return ScalarFunction(static_cast<const Polynomial<Interval>&>(*p1) + static_cast<const Polynomial<Interval>&>(*p2) );
    }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator-(const ScalarFunction& f1, const ScalarFunction& f2)
{
    const ScalarPolynomialFunction* p1=dynamic_cast<const ScalarPolynomialFunction*>(f1.pointer());
    const ScalarPolynomialFunction* p2=dynamic_cast<const ScalarPolynomialFunction*>(f2.pointer());
    if(p1 && p2) {
        return ScalarFunction(static_cast<const Polynomial<Interval>&>(*p1) - static_cast<const Polynomial<Interval>&>(*p2) );
    }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator*(const ScalarFunction& f1, const ScalarFunction& f2)
{
    const ScalarPolynomialFunction* p1=dynamic_cast<const ScalarPolynomialFunction*>(f1.pointer());
    const ScalarPolynomialFunction* p2=dynamic_cast<const ScalarPolynomialFunction*>(f2.pointer());
    if(p1 && p2) {
        return ScalarFunction(static_cast<const Polynomial<Interval>&>(*p1) * static_cast<const Polynomial<Interval>&>(*p2) );
    }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator/(const ScalarFunction& f1, const ScalarFunction& f2)
{
    ARIADNE_NOT_IMPLEMENTED;
}


ScalarFunction operator+(const ScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunction* p1=dynamic_cast<const ScalarPolynomialFunction*>(f1.pointer());
    if(p1) { return ScalarFunction(static_cast<const Polynomial<Interval>&>(*p1) + static_cast<Interval>(s2) ); }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator-(const ScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunction* p1=dynamic_cast<const ScalarPolynomialFunction*>(f1.pointer());
    if(p1) { return ScalarFunction(static_cast<const Polynomial<Interval>&>(*p1) - static_cast<Interval>(s2) ); }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator*(const ScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunction* p1=dynamic_cast<const ScalarPolynomialFunction*>(f1.pointer());
    if(p1) { return ScalarFunction(static_cast<const Polynomial<Interval>&>(*p1) * static_cast<Interval>(s2) ); }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator/(const ScalarFunction& f1, const Real& s2)
{
    const ScalarPolynomialFunction* p1=dynamic_cast<const ScalarPolynomialFunction*>(f1.pointer());
    if(p1) { return ScalarFunction(static_cast<const Polynomial<Interval>&>(*p1) / static_cast<Interval>(s2) ); }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator+(const Real& s1, const ScalarFunction& f2)
{
    return f2+s1;
}

ScalarFunction operator-(const Real& s1, const ScalarFunction& f2)
{
    const ScalarPolynomialFunction* p2=dynamic_cast<const ScalarPolynomialFunction*>(f2.pointer());
    if(p2) { return ScalarFunction( static_cast<Interval>(s1) - static_cast<const Polynomial<Interval>&>(*p2) ); }
    ARIADNE_FAIL_MSG("");
}

ScalarFunction operator*(const Real& s1, const ScalarFunction& f2)
{
    return f2*s1;
}

ScalarFunction operator/(const Real& s1, const ScalarFunction& f2)
{
    ARIADNE_NOT_IMPLEMENTED;
}





VectorFunction::VectorFunction(Nat rs, Nat as)
    : _rep(Vector<ScalarFunction>(rs,ScalarFunction::constant(as,Real(0.0))))
{
}

VectorFunction VectorFunction::constant(const Vector<Real>& c, Nat as)
{
    const Nat rs=c.size();
    VectorFunction res(rs,as);
    for(uint i=0; i!=rs; ++i) {
        res._rep[i]=ScalarFunction::constant(as,c[i]);
    }
    return res;
}

VectorFunction VectorFunction::identity(Nat n)
{
    VectorFunction res(n,n);
    for(uint i=0; i!=n; ++i) {
        res._rep[i]=ScalarFunction::variable(n,i);
    }
    return res;
}


ScalarFunction VectorFunction::get(Nat i) const
{
    return _rep[i];
}

void VectorFunction::set(Nat i, ScalarFunction f)
{
    _rep[i]=f;
}

ScalarFunction VectorFunction::operator[](Nat i) const
{
    return _rep[i];
}

ScalarFunction& VectorFunction::operator[](Nat i)
{
    return _rep[i];
}

Nat VectorFunction::result_size() const
{
    return this->_rep.size();
}

Nat VectorFunction::argument_size() const
{
    if(this->_rep.size()==0u) { return 0; }
    return this->_rep[0].argument_size();
}

Vector<Float> VectorFunction::operator()(const Vector<Float>& x) const
{
    Vector<Float> r(this->result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=this->_rep[i](x);
    }
    return r;
}

Vector<Interval> VectorFunction::operator()(const Vector<Interval>& x) const
{
    Vector<Interval> r(this->result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=this->_rep[i](x);
    }
    return r;
}

std::ostream& VectorFunction::write(std::ostream& os) const
{
    return os << this->_rep;
}





} // namespace Ariadne
