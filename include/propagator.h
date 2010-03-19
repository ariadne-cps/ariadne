/***************************************************************************
 *            propagator.h
 *
 *  Copyright  2009  Pieter Collins
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

/*! \file propagator.h
 *  \brief Class for forward-backward constraint propagation.
 */

#ifndef ARIADNE_PROPAGATOR_H
#define ARIADNE_PROPAGATOR_H

#include "pointer.h"

#include "operators.h"
#include "numeric.h"
#include "real.h"
#include "vector.h"

namespace Ariadne {

inline
void restrict(Interval& r, const Interval& x) {
    r.set_lower(max(r.lower(),x.lower()));
    r.set_upper(min(r.upper(),x.upper()));
};

template<class X> class PropagatorBody;
typedef int Int;

template<class X>
class Propagator {
    typedef Vector<X> ValuationType;
    typedef uint VariableType;
  public:
    Propagator();
    Propagator(int c);
    Propagator(double c);
    Propagator(const Real& c);
    Propagator(const X& c);
    Propagator(Operator op, const Propagator& a);
    Propagator(Operator op, const Propagator& a1, const Propagator& a2);

    Propagator<X>& operator=(const int& c);
    Propagator<X>& operator=(const Real& c);

    static Vector<Propagator<X> > variables(uint n);
    static Propagator<X> variable(uint i);

    Operator op() const;
    const X& value() const;
    X& value();
    const Propagator<X>& arg1() const;
    Propagator<X>& arg1();
    const Propagator<X>& arg2() const;
    Propagator<X>& arg2();

    ValuationType reduce(const ValuationType& x);

    std::ostream& repr(std::ostream& os) const;
    std::ostream& write(std::ostream& os) const;

  public:
    Propagator(PropagatorBody<X>* ptr);
    Propagator<X>& assign(const ValuationType& x);
    Propagator<X>& propagate(ValuationType& x);
  public:
    shared_ptr< PropagatorBody<X> > _ptr;
};

template<class X>
struct PropagatorBody {
    Operator _op;
    Int _var;
    Propagator<X> _arg1;
    Propagator<X> _arg2;
    Interval _val;
    PropagatorBody(Int v) : _op(VAR), _var(v) { }
    PropagatorBody(Interval c) : _op(CNST), _val(c) { }
    PropagatorBody(Operator o, const Propagator<X>& a) : _op(o), _arg1(a) { }
    PropagatorBody(Operator o, const Propagator<X>& a1, const Propagator<X>& a2) : _op(o), _arg1(a1), _arg2(a2) { }
};

template<class X> inline Propagator<X>::Propagator() : _ptr() { }
template<class X> inline Propagator<X>::Propagator(PropagatorBody<X>* ptr) : _ptr(ptr) { }
template<class X> inline Propagator<X>::Propagator(int c) : _ptr(new PropagatorBody<X>(X(c))) { }
template<class X> inline Propagator<X>::Propagator(double c) : _ptr(new PropagatorBody<X>(X(c))) { }
template<class X> inline Propagator<X>::Propagator(const Real& c) : _ptr(new PropagatorBody<X>(X(c))) { }
template<class X> inline Propagator<X>::Propagator(const X& c) : _ptr(new PropagatorBody<X>(c)) { }
template<class X> inline Propagator<X>::Propagator(Operator op, const Propagator<X>& a) : _ptr(new PropagatorBody<X>(op,a)) { }
template<class X> inline Propagator<X>::Propagator(Operator op, const Propagator<X>& a1, const Propagator<X>& a2)
    : _ptr(new PropagatorBody<X>(op,a1,a2)) { }

template<class X> inline Propagator<X>& Propagator<X>::operator=(const int& c) { *this=Propagator<X>(c); return *this; }
template<class X> inline Propagator<X>& Propagator<X>::operator=(const Real& c) { *this=Propagator<X>(c); return *this; }

template<class X> inline Propagator<X> Propagator<X>::variable(VariableType v) { return Propagator<X>(new PropagatorBody<X>(v)); }
template<class X> inline Vector<Propagator<X> > Propagator<X>::variables(VariableType n) {
    Vector< Propagator<X> > r(n,Propagator<X>()); for(uint i=0; i!=n; ++i) { r[i]=Propagator<X>::variable(i); } return r; }

template<class X> inline Operator Propagator<X>::op() const { return this->_ptr->_op; }
template<class X> inline const X& Propagator<X>::value() const { return this->_ptr->_val; }
template<class X> inline X& Propagator<X>::value() { return this->_ptr->_val; }
template<class X> inline const Propagator<X>& Propagator<X>::arg1() const { return this->_ptr->_arg1; }
template<class X> inline Propagator<X>& Propagator<X>::arg1() { return this->_ptr->_arg1; }
template<class X> inline const Propagator<X>& Propagator<X>::arg2() const { return this->_ptr->_arg2; }
template<class X> inline Propagator<X>& Propagator<X>::arg2() { return this->_ptr->_arg2; }


template<class X> inline Propagator<X> operator+(const Propagator<X>& a) { return Propagator<X>(POS,a); }
template<class X> inline Propagator<X> operator-(const Propagator<X>& a) { return Propagator<X>(NEG,a); }
template<class X> inline Propagator<X> operator+(const Propagator<X>& a1, const Propagator<X>& a2) { return Propagator<X>(ADD,a1,a2); }
template<class X> inline Propagator<X> operator-(const Propagator<X>& a1, const Propagator<X>& a2) { return Propagator<X>(SUB,a1,a2); }
template<class X> inline Propagator<X> operator*(const Propagator<X>& a1, const Propagator<X>& a2) { return Propagator<X>(MUL,a1,a2); }
template<class X> inline Propagator<X> operator/(const Propagator<X>& a1, const Propagator<X>& a2) { return Propagator<X>(DIV,a1,a2); }
template<class X> inline Propagator<X> neg(const Propagator<X>& a) { return Propagator<X>(NEG,a); }
template<class X> inline Propagator<X> rec(const Propagator<X>& a) { return Propagator<X>(REC,a); }
template<class X> inline Propagator<X> sqr(const Propagator<X>& a) { return Propagator<X>(SQR,a); }
template<class X> inline Propagator<X> sqrt(const Propagator<X>& a) { return Propagator<X>(SQRT,a); }
template<class X> inline Propagator<X> exp(const Propagator<X>& a) { return Propagator<X>(EXP,a); }
template<class X> inline Propagator<X> log(const Propagator<X>& a) { return Propagator<X>(LOG,a); }
template<class X> inline Propagator<X> sin(const Propagator<X>& a) { return Propagator<X>(SIN,a); }
template<class X> inline Propagator<X> cos(const Propagator<X>& a) { return Propagator<X>(COS,a); }
template<class X> inline Propagator<X> tan(const Propagator<X>& a) { return Propagator<X>(TAN,a); }
template<class X> inline Propagator<X> atan(const Propagator<X>& a) { return Propagator<X>(ATAN,a); }
template<class X> inline Propagator<X> operator==(const Propagator<X>& a1, const Propagator<X>& a2) { return Propagator<X>(EQ,a1,a2); }
template<class X> inline Propagator<X> operator<=(const Propagator<X>& a1, const Propagator<X>& a2) { return Propagator<X>(EQ,a1,a2); }

template<class X> inline Propagator<X> operator+(const Propagator<X>& a1, const Real& c2) { return Propagator<X>(ADD,a1,Propagator<X>(c2)); }
template<class X> inline Propagator<X> operator+(const Real& c1, const Propagator<X>& a2) { return Propagator<X>(ADD,Propagator<X>(c1),a2); }
template<class X> inline Propagator<X> operator-(const Propagator<X>& a1, const Real& c2) { return Propagator<X>(SUB,a1,Propagator<X>(c2)); }
template<class X> inline Propagator<X> operator-(const Real& c1, const Propagator<X>& a2) { return Propagator<X>(SUB,Propagator<X>(c1),a2); }
template<class X> inline Propagator<X> operator*(const Propagator<X>& a1, const Real& c2) { return Propagator<X>(MUL,a1,Propagator<X>(c2)); }
template<class X> inline Propagator<X> operator*(const Real& c1, const Propagator<X>& a2) { return Propagator<X>(MUL,Propagator<X>(c1),a2); }
template<class X> inline Propagator<X> operator/(const Propagator<X>& a1, const Real& c2) { return Propagator<X>(DIV,a1,Propagator<X>(c2)); }
template<class X> inline Propagator<X> operator/(const Real& c1, const Propagator<X>& a2) { return Propagator<X>(DIV,Propagator<X>(c1),a2); }
template<class X> inline Propagator<X>& operator+=(Propagator<X>& a1, const Propagator<X>& a2) { a1=Propagator<X>(ADD,a1,a2); return a1; }
template<class X> inline Propagator<X>& operator+=(Propagator<X>& a1, const Real& c2) { a1=Propagator<X>(ADD,a1,Propagator<X>(c2)); return a1; }
template<class X> inline Propagator<X>& operator*=(Propagator<X>& a1, const Propagator<X>& a2) { a1=Propagator<X>(MUL,a1,a2); return a1; }
template<class X> inline Propagator<X>& operator*=(Propagator<X>& a1, const Real& c2) { a1=Propagator<X>(MUL,a1,Propagator<X>(c2)); return a1; }
template<class X> inline Propagator<X> operator==(const Propagator<X>& a1, const X& c2) { return Propagator<X>(EQ,a1,Propagator<X>(c2)); }

template<class X> inline std::ostream& Propagator<X>::repr(std::ostream& os) const {
    os << this->_ptr << " " << *this << ":" << this->value() << "\n";
    switch(this->_ptr->_op) {
        case CNST: case VAR: break;
        default: this->arg1().repr(os); os << this->_ptr->_op; this->arg2().repr(os);
    }
    return os;
}


template<class X> inline std::ostream& operator<<(std::ostream& os, const Propagator<X>& p) {
    return p.write(os); }

template<class X>
Propagator<X>& Propagator<X>::assign(const ValuationType& x)
{
    //Interval& v=this->value();
    Propagator<X>& a1=this->arg1();
    Propagator<X>& a2=this->arg2();
    switch(_ptr->_op) {
        case CNST: break;
        case VAR: value()=x[this->_ptr->_var]; break;
        case ADD: value()=a1.assign(x).value()+a2.assign(x).value(); break;
        case SUB: value()=a1.assign(x).value()-a2.assign(x).value(); break;
        case MUL: value()=a1.assign(x).value()*a2.assign(x).value(); break;
        case DIV: value()=a1.assign(x).value()/a2.assign(x).value(); break;
        case POS: value()=a1.assign(x).value(); break;
        case NEG: value()=-a1.assign(x).value(); break;
        case REC: value()=rec(a1.assign(x).value()); break;
        case SQR: value()=sqr(a1.assign(x).value()); break;
        case SQRT: value()=sqrt(a1.assign(x).value()); break;
        case SIN: value()=sin(a1.assign(x).value()); break;
        case COS: value()=cos(a1.assign(x).value()); break;
        case ATAN: value()=atan(a1.assign(x).value()); break;
        case EQ: value()=intersection(a1.assign(x).value(),a2.assign(x).value()); break;
        case LEQ: value()=Interval(a1.assign(x).value().lower(),a2.assign(x).value().upper()); break;
        default: ARIADNE_THROW(std::runtime_error,"Propagator<X>::propagate","Unhandled operator "<<_ptr->_op);
    }
    return *this;
}

template<class X>
Propagator<X>& Propagator<X>::propagate(Vector<X>& x)
{
    static const Float inf=Ariadne::inf<Float>();

    const Interval& v0=this->value();
    Propagator<X>& a1=this->arg1();
    Propagator<X>& a2=this->arg2();
    Interval t1,t2;

    switch(_ptr->_op) {
        case CNST: break;
        case VAR: restrict(x[_ptr->_var],v0); break;
        case ADD: restrict(a1.value(),v0-a2.value()); restrict(a2.value(),v0-a1.value()); break;
        case SUB: restrict(a1.value(),a2.value()+v0); restrict(a2.value(),a1.value()-v0); break;
        case MUL: restrict(a1.value(),v0/a2.value()); restrict(a2.value(),v0/a1.value()); break;
        case DIV: restrict(a1.value(),a2.value()*v0); restrict(a2.value(),a1.value()/v0); break;
        case POS: restrict(a1.value(),v0); break;
        case NEG: restrict(a1.value(),neg(v0)); break;
        case REC: restrict(a1.value(),rec(v0)); break;
        case SQR: restrict(a1.value(),sqrt(v0)); break;
        case SQRT: restrict(a1.value(),sqr(v0)); break;
        case EXP: restrict(a1.value(),log(v0)); break;
        case LOG: restrict(a1.value(),exp(v0)); break;
        case SIN: restrict(a1.value(),asin(v0)); break;
        case COS: restrict(a1.value(),acos(v0)); break;
        case ATAN: restrict(a1.value(),tan(v0)); break;
        case EQ: restrict(a1.value(),v0); restrict(a2.value(),v0); break;
        case LEQ: restrict(a1.value(),Interval(-inf,a2.value().upper())); restrict(a1.value(),Interval(a2.value().lower(),+inf)); break;
        default: ARIADNE_THROW(std::runtime_error,"Propagator<X>::propagate","Unhandled operator "<<_ptr->_op);
    }

    if(this->_ptr->_arg1._ptr.get()) {
        a1.propagate(x);
        if(this->_ptr->_arg2._ptr.get()) {
            a2.propagate(x);
        }
    }
    return *this;
}

template<class X>
Vector<X> Propagator<X>::reduce(const Vector<X>& x)
{
    Vector<X> r(x);
    this->assign(r);
    this->propagate(r);
    return r;
}


template<class X>
std::ostream& Propagator<X>::write(std::ostream& os) const
{
    const Propagator<X>& p=*this;
    if(p.op()==CNST) { return os << p.value(); }
    if(p.op()==VAR) { return os << 'x' << p._ptr->_var; }
    os << p.arg1();
    switch(p.op()) {
        case ADD: os << "+"; break;
        case SUB: os << "-"; break;
        case MUL: os << "*"; break;
        case DIV: os << "/"; break;
        case EQ: os << "=="; break;
        case LEQ: os << "<="; break;
        default: break;
    }
    os << p.arg2();
    return os;
}



} // namespace Ariadne

#endif /* ARIADNE_PROPAGATOR_H */

