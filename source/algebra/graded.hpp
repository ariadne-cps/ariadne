/***************************************************************************
 *            graded.hpp
 *
 *  Copyright  2011  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file graded.hpp
 *  \brief Graded algebras.
 */

#ifndef ARIADNE_GRADED_HPP
#define ARIADNE_GRADED_HPP

#include "../function/procedure.hpp"

namespace Ariadne {

struct AntiDiff { };

template<class Op, class A1, class A2=Void, class A3=Void> struct ClosureExpression;

template<class Op, class A> struct ClosureExpression<Op,A,Void,Void> {
    const A& arg;
    ClosureExpression(const A& a) : arg(a) { }
};

template<class Op, class A1, class A2> struct ClosureExpression<Op,A1,A2,Void> {
    const A1& arg1; const A2& arg2;
    ClosureExpression(const A1& a1, const A2& a2) : arg1(a1), arg2(a2) { }
};

template<class Op, class A1, class A2, class A3> struct ClosureExpression {
    const A1& arg1; const A2& arg2; const A2& arg3;
    ClosureExpression(const A1& a1, const A2& a2, const A3& a3) : arg1(a1), arg2(a2), arg3(a3) { }
};

template<class Op, class A> ClosureExpression<Op,A>
make_expression(Op op, const A& a) {
    return ClosureExpression<Op,A>(a);
}

template<class Op, class A1, class A2> ClosureExpression<Op,A1,A2>
make_expression(Op op, const A1& a1, const A2& a2) {
    return ClosureExpression<Op,A1,A2>(a1,a2);
}

template<class Op, class A1, class A2, class A3> ClosureExpression<Op,A1,A2,A3>
make_expression(Op op, const A1& a1, const A2& a2, const A3& a3) {
    return ClosureExpression<Op,A1,A2,A3>(a1,a2,a3);
}



inline Bool compatible(const FloatDP& x1, const FloatDP& x2) { return true; }
template<class X> inline Bool compatible(const Differential<X>& x1, const Differential<X>& x2) { return x1.argument_size()==x2.argument_size(); }

inline FloatDP create(const FloatDP& x) { return FloatDP(0); }
inline ExactIntervalType create(const ExactIntervalType& x) { return ExactIntervalType(0); }
template<class X> inline Differential<X> create(const Differential<X>& x) { return Differential<X>(x.argument_size(),x.degree()); }

template<class X> inline X create(const X& x) { return nul(x); }

template<class A> class Graded : public List<A>
{
  public:
    typedef typename A::NumericType NumericType;
    typedef Graded<A> SelfType;
    Graded() : List<A>() { }
    Graded(const A& a) : List<A>(1u,a) { }
    Graded(const List<A>& l) : List<A>(l) { }
    template<class Op> Void operator=(const ClosureExpression<Op,SelfType,SelfType>& expr);
    template<class Op> Void operator=(const ClosureExpression<Op,SelfType>& expr);
    Void operator=(const ClosureExpression<AntiDiff,SelfType>& ad);
    Graded<A>& operator=(const Graded<A>& a) { this->List<A>::operator=(a); return *this; }
    Graded<A> create_zero() const { return Graded<A>(List<A>(this->degree()+1, Ariadne::create_zero((*this)[0]))); }
    Nat degree() const { return this->size()-1u; }
    Void extend(const A& a) { this->List<A>::append(a); }
};
template<class A> OutputStream& operator<<(OutputStream& os, const Graded<A>& g) {
    if(g.size()==0) { return os << "G[-]{}"; }
    os << "G[" << g.degree() << "]{";
    os << "(" << g[0] << ")";
    for(Nat i=1; i<=g.degree(); ++i) {
        os << " + (" << g[i] << ")*t";
        if(i>1) { os << "^"<<i; }
    }
    os << "}";
    return os;
}
template<> inline OutputStream& operator<<(OutputStream& os, const Graded<FloatDP>& g) {
    if(g.size()==0) { return os << "G[-]{}"; }
    os << "G[" << g.degree() << "]{";
    Bool nonzero=false;
    if(g[0]!=0) { os << g[0]; nonzero=true; }
    for(Nat i=1; i<=g.degree(); ++i) {
        if(g[i]>0) {
            if(nonzero) { os << "+"; }
            if(g[i]==1) { os << "t"; } else { os << g[i]<<"*t"; }
            if(i>1) { os << "^"<<i; }
            nonzero=true;
        } else if(g[i]<0) {
            if(g[i]==-1) { os << "-t"; } else { os << g[i]<<"*t"; }
            if(i>1) { os << "^"<<i; }
            nonzero=true;
         }
    }
    if(!nonzero) { os << "0"; }
    os << "}";
    return os;
}
template<> inline OutputStream& operator<<(OutputStream& os, const Graded<ExactIntervalType>& g) {
    if(g.size()==0) { return os << "0"; }
    os << g[0];
    for(Nat i=1; i<=g.degree(); ++i) {
        os << "+"<<g[i]<<"*t";
        if(i>1) { os << "^"<<i; }
    }
    return os;
}

template<class X> Void compute(X& r, const Neg&, const X& a) { return neg(r,a); }
template<class X> Void compute(X& r, const Add&, const X& a1, const X& a2) { return add(r,a1,a2); }
template<class X> Void compute(X& r, const Sub&, const X& a1, const X& a2) { return sub(r,a1,a2); }
template<class X> Void compute(X& r, const Mul&, const X& a1, const X& a2) { return mul(r,a1,a2); }
template<class X> Void compute(X& r, const Div&, const X& a1, const X& a2) { return div(r,a1,a2); }
template<class X> Void compute(X& r, const Sqr&, const X& a) { return sqr(r,a); }
template<class X> Void compute(X& r, const Sqrt&, const X& a) { return sqrt(r,a); }
template<class X> Void compute(X& r, const Exp&, const X& a) { return exp(r,a); }
template<class X> Void compute(X& r, const Log&, const X& a) { return log(r,a); }
template<class X> Void compute(X& r, const Rec&, const X& a) { return rec(r,a); }
template<class X> Void compute(X& r, const Sin&, const X& a) { return sin(r,a); }
template<class X> Void compute(X& r, const Cos&, const X& a) { return cos(r,a); }

template<class A, class B> Graded<A>& operator+=(Graded<A>& a, const B& c) {
    ARIADNE_ASSERT(a.size()>0); if(a.degree()==0) { a[0]+=c; } return a; }
template<class A, class B> Graded<A>& operator*=(Graded<A>& a, const B& c) {
    a.back()*=c; return a; }

template<class A, class B> Graded<A> operator+(const Graded<A>& a, const B& c) {
    Graded<A> r(a); r[0]+=c; return r;  }
template<class A, class B> Graded<A> operator*(const Graded<A>& a, const B& c) {
    Graded<A> r(a); r*=c; return r;  }
template<class A> Graded<A> operator*(const Graded<A>& a, Int c) {
    Graded<A> r(a); if(c==0) { for(Nat i=0; i!=r.size(); ++i) { r[i]*=c; } } else { r*=c; } return r;  }

template<class A> Graded<A> operator+(const GenericNumericType<A>& c, const Graded<A>& a) {
    Graded<A> r(a); r[0]+=c; return r;  }
template<class A> Graded<A> operator-(const GenericNumericType<A>& c, const Graded<A>& a) {
    ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> operator*(const GenericNumericType<A>& c, const Graded<A>& a) {
    Graded<A> r(a); r[0]*=c; return r;  }
template<class A> Graded<A> operator/(const GenericNumericType<A>& c, const Graded<A>& a) {
    ARIADNE_NOT_IMPLEMENTED; }

template<class A> ClosureExpression<Neg,Graded<A>> operator-(const Graded<A>& a) {
    return make_expression(Neg(),a); }
template<class A> ClosureExpression<Add,Graded<A>,Graded<A>> operator+(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Add(),a1,a2); }
template<class A> ClosureExpression<Sub,Graded<A>,Graded<A>> operator-(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Sub(),a1,a2); }
template<class A> ClosureExpression<Mul,Graded<A>,Graded<A>> operator*(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Mul(),a1,a2); }
template<class A> ClosureExpression<Div,Graded<A>,Graded<A>> operator/(const Graded<A>& a1, const Graded<A>& a2) {
    return make_expression(Div(),a1,a2); }
template<class A> ClosureExpression<Neg,Graded<A>> neg(const Graded<A>& a) {
    return make_expression(Neg(),a); }
template<class A> ClosureExpression<Sqr,Graded<A>> sqr(const Graded<A>& a) {
    return make_expression(Sqr(),a); }
template<class A> ClosureExpression<Sqrt,Graded<A>> sqrt(const Graded<A>& a) {
    return make_expression(Sqrt(),a); }
template<class A> ClosureExpression<Exp,Graded<A>> exp(const Graded<A>& a) {
    return make_expression(Exp(),a); }
template<class A> ClosureExpression<Log,Graded<A>> log(const Graded<A>& a) {
    return make_expression(Log(),a); }
template<class A> ClosureExpression<Rec,Graded<A>> rec(const Graded<A>& a) {
    return make_expression(Rec(),a); }
template<class A> ClosureExpression<Sin,Graded<A>> sin(const Graded<A>& a) {
    return make_expression(Sin(),a); }
template<class A> ClosureExpression<Cos,Graded<A>> cos(const Graded<A>& a) {
    return make_expression(Cos(),a); }

template<class A> Graded<A> pos(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> abs(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> pow(const Graded<A>& a, Int n) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> tan(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }
template<class A> Graded<A> atan(const Graded<A>& a) { ARIADNE_NOT_IMPLEMENTED; }

template<class A> Void neg(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree()+1u == a.degree());
    r.append(-a.back());
}

template<class A> Void add(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT(r.degree()+1u == a1.degree());
    ARIADNE_ASSERT(r.degree()+1u == a2.degree());
    r.append(a1.back()+a2.back());
}

template<class A> Void sub(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT(r.degree()+1u == a1.degree());
    ARIADNE_ASSERT(r.degree()+1u == a2.degree());
    r.append(a1.back()-a2.back());
}

template<class A> Void mul(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT_MSG(r.degree()+1u == a1.degree(),"r="<<r<<", a1="<<a1<<", a2="<<a2<<"\n");
    ARIADNE_ASSERT_MSG(r.degree()+1u == a2.degree(),"r="<<r<<", a1="<<a1<<", a2="<<a2<<"\n");
    ARIADNE_ASSERT(compatible(a1[0],a2[0]));
    r.append(create(a1[0]));
    Nat d = r.degree();
    for(Nat i=0; i<=d; ++i) {
        r[d] += a1[i]*a2[d-i];
    }
}

template<class A> Void div(Graded<A>& r, const Graded<A>& a1, const Graded<A>& a2) {
    ARIADNE_ASSERT(r.degree()+1u == a1.degree());
    ARIADNE_ASSERT(r.degree()+1u == a2.degree());
    ARIADNE_ASSERT(compatible(a1[0],a2[0]));
    r.append(create(a1[0]));
    Nat d = r.degree();
    r[d]+=a1[d];
    for(Nat i=0; i!=d; ++i) {
        r[d] -= a2[d-i]*r[i];
    }
    r[d]=r[d]/a2[0];
}

template<class A> Void sqr(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.size() < a.size());
    r.append(create(a[0]));
    Nat d = r.degree();
    for(Nat i=0; i<=d; ++i) {
        r[d] += a[i]*a[d-i];
    }
}

template<class A> Void rec(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree()+1u == a.degree());
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]=rec(a[0]); return; }
    for(Nat i=0; i!=d; ++i) {
        r[d] -= a[d-i]*r[i];
    }
    r[d]=r[d]*r[0];
}

template<class A> Void sqrt(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree()+1u == a.degree());
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]=sqrt(a[0]); return; }
    r[d]=a[d];
    for(Nat i=1; i!=d; ++i) {
        r[d] -= r[d-i]*r[i];
    }
    r[d]=r[d]/(2*r[0]);
}

template<class A> Void exp(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree()+1u == a.degree());
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]+=exp(a[0]); return; }
    for(Nat i=0; i!=d; ++i) {
        r[d] += (d-i)*a[d-i]*r[i];
    }
    r[d]/=d;
}

template<class A> Void log(Graded<A>& r, const Graded<A>& a) {
    // y=log x; r=1/x; s=r^2;
    r.append(create(a[0]));
    Nat d = r.degree();
    if(d==0) { r[d]=log(a[0]); return; }
    if(d==1) { r[d]=a[d]/a[0];  return; }
    r[d]+=d*a[d];
    for(Nat i=1; i!=d; ++i) {
        r[d] -= a[d-i]*r[i]*i;
    }
    r[d]=r[d]/a[0];
    r[d]/=d;
}

/*
template<class A> Void sin(Graded<A>& r, const Graded<A>& a) {
    // Let f[0](t)=f(t), s[0](t)=sin(f(t)), c[0](t)=cos(f(t))
    // f[n+1](t)=df[n](t)/dt, s[n+1]=ds[n](t)/dt, c[n+1]=dc[n](t)/dt
    // Then s[1]=c[0]*f[1]; c[1]=-s[0]*f[1]
    const Graded<A>& f=a;
    Graded<A>& s=r;
    A z=create(f[0])*0;
    s.append(z);
    Nat d = s.degree();
    if(d==0) { s[d]=sin(f[0]); return; }
    if(d==1) { s[d]=cos(f[0])*f[1]; return; }
    Graded<A> c(cos(f[0]));
    c.append(-s[0]*f[1]);
    for(Nat i=2; i!=d; ++i) {
        c.append(z);
        for(Nat j=0; j!=i; ++j) {
            c.back() -= (i-j)*f[i-j]*s[j];
        }
        c.back()/=i;
    }
    for(Nat i=0; i!=d; ++i) {
        s[d] += (d-i)*f[d-i]*c[i];
    }
    s[d]/=d;
}
*/

template<class A> Void sincos(Graded<A>& s, Graded<A>& c, const Graded<A>& f) {
    // Let f[0](t)=f(t), s[0](t)=sin(f(t)), c[0](t)=cos(f(t))
    // f[n+1](t)=df[n](t)/dt, s[n+1]=ds[n](t)/dt, c[n+1]=dc[n](t)/dt
    // Then s[n]=+Sum_{m=1}^{n} (m*f[m]*c[n-m])/n
    // Then c[n]=-Sum_{m=1}^{n} (m*f[m]*s[n-m])/n
    Nat d = f.degree();
    A z=create(f[0]);
    ARIADNE_ASSERT(s.size()==0 && c.size()==0);
    for(Nat i=0; i<=d; ++i) { s.append(z); c.append(z); }
    s[0]=sin(f[0]);
    c[0]=cos(f[0]);
    for(Nat i=1; i<=d; ++i) {
        for(Nat j=1; j<=i; ++j) {
            s[i] += j*f[j]*c[i-j];
            c[i] += j*f[j]*s[i-j];
        }
        s[i]/=(+i);
        c[i]/=(-i);
    }
}

template<class A> Void sin(Graded<A>& r, const Graded<A>& a) {
    Graded<A> s;
    Graded<A> c;
    sincos(s,c,a);
    r=s;
}

template<class A> Void cos(Graded<A>& r, const Graded<A>& a) {
    Graded<A> s;
    Graded<A> c;
    sincos(s,c,a);
    r=c;
}


template<class A> template<class Op> Void Graded<A>::operator=(const ClosureExpression<Op,Graded<A>,Graded<A>>& expr) {
    compute(*this,Op(),expr.arg1,expr.arg2);
}

template<class A> template<class Op> Void Graded<A>::operator=(const ClosureExpression<Op,Graded<A>>& expr) {
    compute(*this,Op(),expr.arg);
}

template<class A> Void Graded<A>::operator=(const ClosureExpression<AntiDiff,Graded<A>>& expr) {
    antidifferentiate(*this,expr.arg);
}

template<class A> Void antidifferentiate(Graded<A>& r, const Graded<A>& a) {
    ARIADNE_ASSERT(r.degree() == a.degree());
    r.append(a.back()/(r.degree()+1u));
}

template<class A> ClosureExpression<AntiDiff,Graded<A>> antidifferential(const Graded<A>& a) {
    return ClosureExpression<AntiDiff,Graded<A>>(a);
}


Pair<List<FloatDP>,FloatDP> inline midpoint_error(const Graded<ExactIntervalType>& x) {
    List<FloatDP> m(x.degree()+1);
    FloatDP e;
    for(Nat i=0; i<=x.degree(); ++i) {
        m[i]=midpoint(x[i]).raw();
        e=add(up,e,max(sub(up,m[i],x[i].lower().raw()),sub(up,x[i].upper().raw(),m[i])));
    }
    return Pair<List<FloatDP>,FloatDP>(m,e);
}


template<class X> Differential<X> make_differential_variable(Nat as, Nat deg, X val, Nat ind) {
    return Differential<X>::variable(as,deg,val,ind); }
template<class X> Graded<X> make_graded(const X& val) {
    return Graded<X>(val); }
template<class X> Graded<X> create_graded(const X&) {
    return Graded<X>(); }



template<class X, class A> Void compute(const Vector<Procedure<X>>& p, Vector<Graded<A>>& r, List<Graded<A>>& t, const Vector<Graded<A>>& a) {
    execute(t,p,a);
    for(Nat i=0; i!=p._results.size(); ++i) { r[i]=t[p._results[i]]; }
}

} // namespace Ariadne

#endif /* ARIADNE_GRADED_HPP */
