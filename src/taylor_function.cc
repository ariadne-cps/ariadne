/***************************************************************************
 *            taylor_function.cc
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

#include <iostream>
#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "polynomial.h"
#include "differential.h"

#include "function_interface.h"
#include "taylor_function.h"

namespace Ariadne {

typedef unsigned int uint;

typedef Vector<Float> Point;
typedef Vector<Interval> Box;


void _set_scaling(TaylorExpression& x, const Interval& ivl, uint j)
{
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    const double& l=ivl.lower();
    const double& u=ivl.upper();
    volatile double pc=u; pc+=l;
    volatile double nc=-u; nc-=l;
    volatile double pg=u; pg-=l;
    volatile double ng=l; ng-=u;
    x.error()=(pc+nc+pg+ng)/4;
    set_rounding_mode(to_nearest);
    MultiIndex a(x.argument_size());
    x.expansion().append(a,(l+u)/2);
    ++a[j];
    x.expansion().append(a,(l+u)/2);
    set_rounding_mode(rounding_mode);
}




TaylorExpression::TaylorExpression()
    : _domain(0), _model(0)
{ }

TaylorExpression::TaylorExpression(const DomainType& d)
    : _domain(d), _model(d.size())
{
}

TaylorExpression::TaylorExpression(const DomainType& d, const TaylorModel& m)
    : _domain(d), _model(m)
{
}

TaylorExpression::TaylorExpression(const DomainType& d, const ExpansionType& f, const ErrorType& e)
    : _domain(d), _model(f,e)
{
    ARIADNE_ASSERT(d.size()==f.argument_size());
}

TaylorExpression::TaylorExpression(const DomainType& d, const ScalarFunctionInterface& f)
    : _domain(d), _model(f.argument_size())
{
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_model=f.evaluate(x);
}

TaylorExpression::TaylorExpression(const DomainType& d, const Polynomial<Float>& p)
    : _domain(d), _model(p.argument_size())
{
    ARIADNE_ASSERT(d.size()==p.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_model=Ariadne::evaluate(p,x);
}

TaylorExpression::TaylorExpression(const DomainType& d, const Polynomial<Interval>& p)
    : _domain(d), _model(p.argument_size())
{
    ARIADNE_ASSERT(d.size()==p.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_model=Ariadne::evaluate(p,x);
}


TaylorExpression TaylorExpression::constant(const Vector<Interval>& d, const Float& c)
{
    return TaylorExpression(d,TaylorModel::constant(d.size(),c));
}

TaylorExpression TaylorExpression::constant(const Vector<Interval>& d, const Interval& c)
{
    return TaylorExpression(d,TaylorModel::constant(d.size(),c));
}

TaylorExpression TaylorExpression::variable(const Vector<Interval>& d, uint j)
{
    return TaylorExpression(d,TaylorModel::scaling(d.size(),j,d[j]));
}

Vector<TaylorExpression> TaylorExpression::constants(const Vector<Interval>& d, const Vector<Interval>& c)
{
    Vector<TaylorExpression> x(c.size(),TaylorExpression(d));
    for(uint i=0; i!=c.size(); ++i) {
        x[i]=c[i];
    }
    return x;
}

Vector<TaylorExpression> TaylorExpression::variables(const Vector<Interval>& d)
{
    Vector<TaylorExpression> x(d.size());
    for(uint i=0; i!=d.size(); ++i) {
        x[i]=TaylorExpression::variable(d,i);
    }
    return x;
}


Polynomial<Interval> polynomial(const TaylorModel& tm);

Polynomial<Interval>
TaylorExpression::polynomial() const
{
    Polynomial<Interval> p(this->argument_size());
    p=Ariadne::polynomial(this->model());

    Vector<Polynomial<Interval> > s(this->argument_size());
    for(uint j=0; j!=this->argument_size(); ++j) {
        if(this->domain()[j].radius()<=0) { std::cerr<<"Warning: zero radius in domain of TaylorFunction"<<std::endl; }
        else { s[j]=Ariadne::polynomial(TaylorModel::unscaling(this->argument_size(),j,this->domain()[j])); }
    }

    return compose(p,s);
}


bool TaylorExpression::operator==(const TaylorExpression& tv) const
{
    return this->_domain==tv._domain && this->_model==tv._model;
}


TaylorExpression& operator+=(TaylorExpression& x, const TaylorExpression& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model+=y._model; }
    else { x._model+=restrict(y,x.domain())._model; }
    return x;
}

TaylorExpression& operator-=(TaylorExpression& x, const TaylorExpression& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model-=y._model; }
    else { x._model-=restrict(y,x.domain())._model; }
    return x;
}


TaylorExpression operator+(const TaylorExpression& x1, const TaylorExpression& x2) {
    if(x1._domain==x2._domain) {
        return TaylorExpression(x1._domain,x1._model+x2._model); }
    else {
        TaylorExpression::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorExpression(domain,restrict(x1,domain)._model+restrict(x2,domain)._model);}
}

TaylorExpression operator-(const TaylorExpression& x1, const TaylorExpression& x2) {
    if(x1._domain==x2._domain) {
        return TaylorExpression(x1._domain,x1._model-x2._model); }
    else {
        TaylorExpression::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorExpression(domain,restrict(x1,domain)._model-restrict(x2,domain)._model);}
}

TaylorExpression operator*(const TaylorExpression& x1, const TaylorExpression& x2) {
    if(x1._domain==x2._domain) {
        return TaylorExpression(x1._domain,x1._model*x2._model); }
    else {
        TaylorExpression::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorExpression(domain,restrict(x1,domain)._model*restrict(x2,domain)._model);}
}

TaylorExpression operator/(const TaylorExpression& x1, const TaylorExpression& x2) {
    if(x1._domain==x2._domain) {
        return TaylorExpression(x1._domain,x1._model/x2._model); }
    else {
        TaylorExpression::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorExpression(domain,restrict(x1,domain)._model/restrict(x2,domain)._model);}
}





TaylorExpression max(const TaylorExpression& x1, const TaylorExpression& x2) {
    if(x1._domain==x2._domain) {
        return TaylorExpression(x1._domain,max(x1._model,x2._model)); }
    else {
        TaylorExpression::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorExpression(domain,max(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}

TaylorExpression min(const TaylorExpression& x1, const TaylorExpression& x2) {
    if(x1._domain==x2._domain) {
        return TaylorExpression(x1._domain,min(x1._model,x2._model)); }
    else {
        TaylorExpression::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorExpression(domain,min(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}





Interval
TaylorExpression::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}

Interval
TaylorExpression::evaluate(const Vector<Interval>& x) const
{
    Vector<Interval> sx=Ariadne::evaluate(TaylorModel::unscalings(this->_domain),x);
    return Ariadne::evaluate(this->_model,sx);
}

TaylorExpression
partial_evaluate(const TaylorExpression& te, uint k, const Interval& c)
{
    // Scale c to domain
    const uint as=te.argument_size();
    ARIADNE_ASSERT(k<as);
    const Vector<Interval>& domain=te.domain();
    const Interval& dk=domain[k];
    Interval sc=(c-med_ivl(dk))/rad_ivl(dk);

    Vector<Interval> new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    TaylorModel new_model=partial_evaluate(te.model(),k,sc);

    return TaylorExpression(new_domain,new_model);
}


// To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
// and translation t=((c+d)-(a+b))/(b-a)
// Because we are scaling the model on [-1,+1], this is not the same as
// the mapping taking [a,b] to [c,d]
TaylorExpression restrict(const TaylorExpression& tv, uint k, const Interval& new_ivl) {
    ARIADNE_ASSERT(k<tv.argument_size())
    const Interval& old_ivl=tv.domain()[k];
    ARIADNE_ASSERT(subset(new_ivl,old_ivl));
    if(new_ivl==old_ivl) { return tv; }
    Float a=old_ivl.lower(); Float b=old_ivl.upper();
    Float c=new_ivl.lower(); Float d=new_ivl.upper();
    if(a==b) { ARIADNE_ASSERT( a<b || (a==b && c==d) ); return tv; }
    Interval s=sub_ivl(d,c)/sub_ivl(b,a);
    // Interval t=(mul_ivl(b,c)-mul_ivl(a,d))/sub_ivl(b,a);  // WRONG!!
    Interval t=(add_ivl(c,d)-add_ivl(a,b))/sub_ivl(b,a);
    Vector<Interval> new_dom=tv.domain();
    new_dom[k]=new_ivl;
    return TaylorExpression(new_dom,preaffine(tv.model(),k,s,t));
}

TaylorExpression restrict(const TaylorExpression& tv, const Vector<Interval>& d) {
    ARIADNE_ASSERT(subset(d,tv.domain()));
    const Vector<Interval>& od=tv.domain();
    TaylorExpression r=tv;
    for(uint j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { r=restrict(r,j,d[j]); }
    }
    return r;
}

TaylorExpression extend(const TaylorExpression& tv, const Vector<Interval>& d) {
    const Vector<Interval>& domain=tv.domain();
    ARIADNE_ASSERT(subset(domain,d));
    for(uint i=0; i!=d.size(); ++i) {
        ARIADNE_ASSERT(domain[i]==d[i] || domain[i].lower()==domain[i].upper());
    }
    return TaylorExpression(d,tv._model);
}

pair<TaylorExpression,TaylorExpression>
split(const TaylorExpression& tv, uint j)
{
    ARIADNE_NOT_IMPLEMENTED;
}

bool refines(const TaylorExpression& tv1, const TaylorExpression& tv2)
{
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); }
    if(subset(tv2.domain(),tv1.domain())) { return refines(restrict(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}

bool disjoint(const TaylorExpression& tv1, const TaylorExpression& tv2)
{
    if(tv1.domain()==tv2.domain()) {
        return disjoint(tv1.model(),tv2.model());
    } else {
        Vector<Interval> domain=intersection(tv1.domain(),tv2.domain());
        return disjoint(restrict(tv1,domain).model(),restrict(tv2,domain).model());
    }
}

TaylorExpression intersection(const TaylorExpression& tv1, const TaylorExpression& tv2)
{
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return TaylorExpression(tv1.domain(),intersection(tv1.model(),tv2.model()));
}

Vector<TaylorExpression> compose(const Vector<TaylorExpression>& x, const Vector<TaylorExpression>& y) {
    ARIADNE_NOT_IMPLEMENTED; }
TaylorExpression compose(const TaylorExpression& x, const Vector<TaylorExpression>& y) {
    ARIADNE_NOT_IMPLEMENTED; }

Vector<TaylorExpression>
prod(const Matrix<Interval>& A,
     const Vector<TaylorExpression>& x)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(A.column_size()==x.size());
    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }

    Vector<TaylorExpression> r(A.row_size(),TaylorExpression(x[0].domain()));
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.column_size(); ++j) {
            r[i]+=A[i][j]*x[j];
        }
    }
    return r;
}

Matrix<Interval>
jacobian(const Vector<TaylorExpression>& tv, const Vector<Interval>& x);

TaylorExpression
midpoint(const TaylorExpression& f)
{
    TaylorModel tm=f.model();
    tm.set_error(0.0);
    return TaylorExpression(f.domain(),tm);
}


TaylorExpression
compose(const ScalarFunctionInterface& f, const Vector<TaylorExpression>& g)
{
    ARIADNE_ASSERT(f.argument_size()==g.size());
    for(uint i=0; i!=g.size(); ++i) {
        ARIADNE_ASSERT(g[0].domain()==g[i].domain());
    }

    Vector<Interval> gdomain=g[0].domain();
    Vector<TaylorModel> gmodels(g.size());
    for(uint i=0; i!=g.size(); ++i) { gmodels[i]=g[i].model(); }

    return TaylorExpression(gdomain,f.evaluate(gmodels));
}


Vector<TaylorExpression>
implicit(const Vector<TaylorExpression>& f)
{
    ARIADNE_ASSERT(f.size()>0);
    ARIADNE_ASSERT(f.size()<=f[0].argument_size());
    // Solve the equation f(x,h(x))=0
    // Use D1f + D2f Dh = 0, so Dh=-D2f^-1 D1f
    uint rs=f.size();
    uint fas=f[0].argument_size();
    uint has=fas-rs;

    Vector<Interval> h_domain=project(f[0].domain(),range(0u,has));
    Vector<Interval> h_range=project(f[0].domain(),range(has,fas));
    Vector<TaylorExpression> id=TaylorExpression::variables(h_domain);
    Vector<TaylorExpression> h=TaylorExpression::constants(h_domain,h_range);
    //std::cerr<<"\nid="<<id<<"\nh0="<<h0<<"\n";

    for(uint k=0; k!=10; ++k) {
        Vector<Interval> ih_range=join(h_domain,h_range);
        Matrix<Interval> Df=jacobian(f,ih_range);
        //std::cerr<<"  Df="<<Df<<std::endl;
        Matrix<Interval> D2f=project(Df,range(0,rs),range(has,fas));
        //std::cerr<<"  D2f="<<J<<std::endl;
        Matrix<Interval> D2finv=inverse(D2f);

        for(uint i=0; i!=rs; ++i) {
            h[i].set_error(0);
        }
        Vector<TaylorExpression> idh=join(id,h);
        Vector<TaylorExpression> fidxhx=compose(f,idh);
        //std::cerr<<"  f(x,h(x))="<<fh<<std::endl;
        Vector<TaylorExpression> dh=prod(D2finv,fidxhx);
        //std::cerr<<"  dh="<<dh<<std::endl;
        h=h-dh;
    }
    //std::cerr<<"\n  f="<<f<<"\n  h[0]="<<h0<<"\n  h[1]="<<h1<<"\n\n";
    ARIADNE_ASSERT(h.size()==f.size());
    ARIADNE_ASSERT(h[0].argument_size()+h.size()==f[0].argument_size());
    return h;

}

TaylorExpression implicit(const TaylorExpression& f) {
    Vector<Interval> h_domain=project(f.domain(),range(0u,f.argument_size()-1u));
    Interval h_codomain=f.domain()[f.argument_size()-1u];
    TaylorModel h_model=implicit(f.model());
    ARIADNE_ASSERT(h_model.argument_size()+1==f.model().argument_size());
    TaylorModel hrs_model=h_model.rescale(Interval(-1,+1),h_codomain);
    ARIADNE_ASSERT(hrs_model.argument_size()+1==f.model().argument_size());
    return TaylorExpression(h_domain,hrs_model);
}


TaylorExpression
implicit(const ScalarFunctionInterface& f, const Vector<TaylorExpression>& g)
{
    ARIADNE_ASSERT(f.argument_size()>g.size());
    ARIADNE_ASSERT(g.size()>0);
    for(uint i=1; i!=g.size(); ++i) {
        ARIADNE_ASSERT(g[i].domain()==g[0].domain());
    }

    Vector<Interval> domain=g[0].domain();
    Vector<TaylorModel> models(g.size());
    for(uint i=0; i!=g.size(); ++i) {
        models[i]=g[i].model();
    }

    return TaylorExpression(domain,implicit(f,models));
}

TaylorExpression
implicit(const ScalarFunctionInterface& f, const Vector<Interval>& d)
{
    ARIADNE_ASSERT(f.argument_size()>=1u);
    ARIADNE_ASSERT(d.size()+1u==f.argument_size());

    Vector<TaylorExpression> id=TaylorExpression::variables(d);
    //for(uint i=0; i!=d.size(); ++i) { id[i].model().set_accuracy(a); }
    return implicit(f,id);
}


TaylorExpression
crossing_time(const ScalarFunctionInterface& g, const TaylorFunction& phi)
{
    Vector<Interval> d=project(phi.domain(),range(0,phi.result_size()));
    Interval h=phi.domain()[phi.result_size()];

    TaylorExpression gphi(phi.domain(),g.evaluate(phi.models()));

    return implicit(gphi);
}


std::ostream&
operator<<(std::ostream& os, const TaylorExpression& tv) {
    return os << "TaylorExpression( domain=" << tv.domain() << ", polynomial=" << midpoint(tv.polynomial()) << "+/-" << tv.error() << ", model=" << tv.model() << ")";
}





bool
check(const Vector<TaylorExpression>& tv)
{
    for(uint i=0; i!=tv.size(); ++i) {
        if(tv[0].domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

Vector< Expansion<Float> >
expansion(const Vector<TaylorExpression>& x)
{
    Vector< Expansion<Float> > r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

Vector<Float>
error(const Vector<TaylorExpression>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

Vector<Float>
value(const Vector<TaylorExpression>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

Vector<Interval>
ranges(const Vector<TaylorExpression>& x)
{
    Vector<Interval> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}


Vector<Interval>
evaluate(const Vector<TaylorExpression>& tv, const Vector<Interval>& x)
{
    Vector<Interval> r(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}

Matrix<Interval>
jacobian(const Vector<TaylorExpression>& tv, const Vector<Interval>& x)
{
    ARIADNE_ASSERT(check(tv));
    const Vector<Interval>& dom=tv[0].domain();
    const uint n=dom.size();
    Vector< Differential<Interval> > s(n,n,1u);
    for(uint j=0; j!=n; ++j) {
        Interval dj=dom[j];
        s[j].set_value((x[j]-dj.midpoint())/dj.radius());
        s[j].set_gradient(j,1/dj.radius());
    }
    Vector< Expansion<Float> > p=expansion(tv);
    Vector< Differential<Interval> > d=evaluate(p,s);
    //std::cerr<<"  x="<<x<<"\n  p="<<p<<"\n"<<"  s="<<s<<"\n  p.s="<<d<<"\n  J="<<d.jacobian()<<"\n"<<std::endl;
    return d.jacobian();
}





TaylorFunction::TaylorFunction()
    : _domain(), _models()
{
}




TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector<TaylorModel>& f)
    : _domain(d), _models(f)
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
    }
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Expansion<Float> >& f)
    : _domain(d), _models(f.size())
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=TaylorModel(f[i],0.0);
    }
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Expansion<Float> >& f,
                               const Vector<Float>& e)
    : _domain(d), _models(f.size())
{
    ARIADNE_ASSERT(f.size()==e.size());
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=TaylorModel(f[i],e[i]);
    }
}



TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const FunctionInterface& f)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=f.evaluate(x);
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const FunctionInterface& f,
                               const shared_ptr<TaylorModel::Accuracy> accuracy_ptr)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    for(uint i=0; i!=x.size(); ++i) { x[i].accuracy_ptr()=accuracy_ptr; }
    this->_models=f.evaluate(x);
}


TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Float> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=Ariadne::evaluate(p,x);
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Interval> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=Ariadne::evaluate(p,x);
}

TaylorFunction::TaylorFunction(const Vector<TaylorExpression>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(uint i=1; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v[0].domain()); }
    this->_domain=v[0].domain();
    for(uint i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}










TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Float>& c)
{
    return TaylorFunction(d,TaylorModel::constants(d.size(),c));
}

TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Interval>& c)
{
    return TaylorFunction(d,TaylorModel::constants(d.size(),c));
}

TaylorFunction
TaylorFunction::identity(const Vector<Interval>& d)
{
    return TaylorFunction(d,TaylorModel::scalings(d));
}


Polynomial<Interval> polynomial(const TaylorModel& tm) {
    return Polynomial<Interval>(tm.expansion())+Interval(-tm.error(),+tm.error());
}

Vector< Polynomial<Interval> >
TaylorFunction::polynomial() const
{
    Vector<Polynomial<Interval> > p(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        p[i]=Ariadne::polynomial(this->models()[i]);
    }

    Vector<Polynomial<Interval> > s(this->argument_size());
    for(uint j=0; j!=this->argument_size(); ++j) {
        if(this->domain()[j].radius()<=0) { std::cerr<<"Warning: zero radius in domain of TaylorFunction"<<std::endl; }
        else { s[j]=Ariadne::polynomial(TaylorModel::unscaling(this->argument_size(),j,this->domain()[j])); }
    }

    return compose(p,s);
}

bool
TaylorFunction::operator==(const TaylorFunction& tm) const
{
    return this->_models==tm._models;
}



bool
TaylorFunction::operator!=(const TaylorFunction& p2) const
{
    return !(*this==p2);
}



shared_ptr<TaylorModel::Accuracy>
TaylorFunction::accuracy_ptr() const
{
    return this->_models[0].accuracy_ptr();
}


void
TaylorFunction::set_accuracy(shared_ptr<TaylorModel::Accuracy> acc_ptr)
{
    for(uint i=0; i!=this->result_size(); ++i) {
        this->_models[i].set_accuracy(acc_ptr);
    }
}

const Vector<Interval>&
TaylorFunction::domain() const
{
    return this->_domain;
}

const Vector<Interval>
TaylorFunction::codomain() const
{
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return result;
}


const Vector<Float>
TaylorFunction::centre() const
{
    Vector<Float> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].value();
    }
    return result;
}

const Vector<Interval>
TaylorFunction::range() const
{
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return result;
}



const Vector<TaylorModel>&
TaylorFunction::models() const
{
    return this->_models;
}




uint
TaylorFunction::argument_size() const
{
    return this->_domain.size();
}


uint
TaylorFunction::result_size() const
{
    return this->_models.size();
}


TaylorExpression
TaylorFunction::operator[](uint i) const
{
    return TaylorExpression(this->_domain,this->_models[i]);
}

TaylorExpression
TaylorFunction::get(uint i) const
{
    return TaylorExpression(this->_domain,this->_models[i]);
}

void
TaylorFunction::set(uint i, const TaylorExpression& e)
{
    ARIADNE_ASSERT_MSG(this->size()>i,"Cannot set "<<i<<"th element of TaylorFunction "<<(*this));
    if(this->domain().size()!=0) {
        ARIADNE_ASSERT_MSG(e.domain()==this->domain(),"Domain of "<<e<<" conflicts with existing domain "<<this->domain());
    } else {
        this->_domain=e.domain();
    }
    this->_models[i]=e.model();
}



















TaylorFunction&
TaylorFunction::truncate(ushort degree)
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].truncate(degree);
    }
    return *this;
}

TaylorFunction&
TaylorFunction::clobber()
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].clobber();
    }
    return *this;
}



Vector<Interval>
TaylorFunction::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}


Vector<Interval>
TaylorFunction::evaluate(const Vector<Interval>& x) const
{
    Vector<Interval> sx=Ariadne::evaluate(TaylorModel::unscalings(this->_domain),x);
    return Ariadne::evaluate(this->_models,sx);
}


Matrix<Interval>
TaylorFunction::jacobian(const Vector<Interval>& x) const
{
    Matrix<Interval> J=Ariadne::jacobian(this->_models,unscale(x,this->_domain));
    for(uint j=0; j!=J.column_size(); ++j) {
        Interval rad=rad_ivl(this->_domain[j]);
        for(uint i=0; i!=J.row_size(); ++i) {
            J[i][j]/=rad;
        }
    }
    return J;
}


TaylorFunction
join(const TaylorFunction& f1, const TaylorExpression& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return TaylorFunction(f1.domain(),join(f1.models(),f2.model()));
}

TaylorFunction
join(const TaylorFunction& f, const TaylorFunction& g)
{
    ARIADNE_ASSERT(f.domain()==g.domain());
    return TaylorFunction(f.domain(),join(f.models(),g.models()));
}

TaylorFunction
join(const TaylorExpression& f1, const TaylorExpression& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return TaylorFunction(f1.domain(),join(f1.model(),f2.model()));
}

TaylorFunction
combine(const TaylorExpression& f1, const TaylorExpression& f2)
{
    return TaylorFunction(join(f1.domain(),f2.domain()),combine(f1.model(),f2.model()));
}

TaylorFunction
combine(const TaylorExpression& f1, const TaylorFunction& f2)
{
    return TaylorFunction(join(f1.domain(),f2.domain()),combine(f1.model(),f2.models()));
}

TaylorFunction
combine(const TaylorFunction& f1, const TaylorExpression& f2)
{
    return TaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.model()));
}

TaylorFunction
combine(const TaylorFunction& f1, const TaylorFunction& f2)
{
    return TaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
}


TaylorFunction
embed(const TaylorFunction& f, const Interval& d)
{
    return TaylorFunction(join(f.domain(),d),embed(f.models(),1u));
}

TaylorFunction
embed(const TaylorFunction& f, const Vector<Interval>& d)
{
    return TaylorFunction(join(f.domain(),d),embed(f.models(),d.size()));
}

TaylorFunction
embed(const Vector<Interval>& d, const TaylorFunction& f)
{
    return TaylorFunction(join(d,f.domain()),embed(d.size(),f.models()));
}

TaylorFunction
restrict(const TaylorFunction& f, const Vector<Interval>& d)
{
    ARIADNE_ASSERT(subset(d,f.domain()));
    if(d==f.domain()) { return f; }
    Vector<TaylorModel> s=TaylorModel::rescalings(f.domain(),d);
    TaylorFunction r(d,compose(f._models,s));
    r.set_accuracy(f.accuracy_ptr());
    return r;
}

bool
refines(const TaylorFunction& f1, const TaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(uint i=0; i!=f1.result_size(); ++i) {
        if(!refines(f1[i],f2[i])) { return false; }
    }
    return true;
}

bool
disjoint(const TaylorFunction& f1, const TaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(uint i=0; i!=f1.result_size(); ++i) {
        if(disjoint(f1[i],f2[i])) { return true; }
    }
    return false;
}

TaylorFunction
intersection(const TaylorFunction& f1, const TaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    TaylorFunction r(f1.result_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r[i]=intersection(f1[i],f2[i]);
    }
    return r;
}

TaylorFunction&
operator+=(TaylorFunction& f, const TaylorFunction& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f._models+=g._models;
    return f;
}

TaylorFunction&
operator-=(TaylorFunction& f, const TaylorFunction& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f._models+=g._models;
    return f;
}

TaylorFunction&
operator+=(TaylorFunction& f, const Vector<Interval>& e)
{
    ARIADNE_ASSERT(f.result_size()==e.size());
    f._models+=e;
    return f;
}

TaylorFunction&
operator-=(TaylorFunction& f, const Vector<Interval>& e)
{
    ARIADNE_ASSERT(f.result_size()==e.size());
    f._models-=e;
    return f;
}

TaylorFunction&
operator*=(TaylorFunction& f, const Float& c)
{
    f._models*=c;
    return f;
}

TaylorFunction&
operator/=(TaylorFunction& f, const Float& c)
{
    f._models/=c;
    return f;
}


TaylorFunction
operator+(const TaylorFunction& f1, const TaylorFunction& f2)
{
    ARIADNE_ASSERT_MSG(!intersection(f1.domain(),f2.domain()).empty(),
                       "operator+(TaylorFunction f1, TaylorFunction f2) with f1="<<f1<<" f2="<<f2<<
                       ": domains are disjoint");
    if(f1.domain()==f2.domain()) {
        return TaylorFunction(f1.domain(),Vector<TaylorModel>(f1.models()+f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator+(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}


TaylorFunction
operator-(const TaylorFunction& f1, const TaylorFunction& f2)
{
    ARIADNE_ASSERT(!intersection(f1.domain(),f2.domain()).empty());
    if(f1.domain()==f2.domain()) {
        return TaylorFunction(f1.domain(),Vector<TaylorModel>(f1.models()+f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}



TaylorFunction
operator-(const TaylorFunction& f)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(-f.models()));
}

TaylorFunction
operator*(const TaylorFunction& f, const Float& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()*c));
}

TaylorFunction
operator/(const TaylorFunction& f, const Float& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()/c));
}

TaylorFunction
operator+(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()+c));
}

TaylorFunction
operator+(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()+c));
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()-c));
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()-c));
}

TaylorFunction
operator*(const Matrix<Float>& A, const TaylorFunction& f)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(prod(A,f.models())));
}

TaylorFunction
operator*(const Matrix<Interval>& A, const TaylorFunction& f)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(boost::numeric::ublas::prod(A,f.models())));
}





TaylorFunction
partial_evaluate(const TaylorFunction& tf, uint k, const Interval& c)
{
    // Scale c to domain
    const uint as=tf.argument_size();
    ARIADNE_ASSERT(k<as);
    const Vector<Interval>& domain=tf.domain();
    const Interval& dk=domain[k];
    Interval sc=(c-med_ivl(dk))/rad_ivl(dk);

    Vector<Interval> new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    Vector<TaylorModel> new_models=partial_evaluate(tf.models(),k,sc);

    return TaylorFunction(new_domain,new_models);
}


TaylorExpression
compose(const ScalarFunctionInterface& g, const TaylorFunction& f)
{
    return TaylorExpression(f.domain(),g.evaluate(f.models()));
}

TaylorFunction
compose(const FunctionInterface& g, const TaylorFunction& f)
{
    return TaylorFunction(f.domain(),g.evaluate(f.models()));
}

TaylorFunction
compose(const TaylorFunction& g, const TaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(subset(f.range(),g.domain()),"f.range()="<<f.range()<<" is not a subset of g.domain()="<<g.domain());
    return TaylorFunction(f.domain(),Ariadne::compose(g.models(),unscale(f.models(),g.domain())));
}


TaylorExpression
compose(const TaylorExpression& g, const TaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(subset(f.range(),g.domain()),"f.range()="<<f.range()<<" is not a subset of g.domain()="<<g.domain());
    return TaylorExpression(f.domain(),Ariadne::unchecked_compose(g.model(),unscale(f.models(),g.domain())));
}

TaylorExpression
unchecked_compose(const TaylorExpression& g, const TaylorFunction& f)
{
        return TaylorExpression(f.domain(),Ariadne::unchecked_compose(g.model(),unscale(f.models(),g.domain())));
}

TaylorFunction
unchecked_compose(const TaylorFunction& g, const TaylorFunction& f)
{
        return TaylorFunction(f.domain(),Ariadne::unchecked_compose(g.models(),unscale(f.models(),g.domain())));
}



TaylorFunction
antiderivative(const TaylorFunction& f, uint k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    Interval fdomkrad=rad_ivl(f.domain()[k].lower(),f.domain()[k].upper());
    TaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g._models[i].antidifferentiate(k);
        g._models[i]*=fdomkrad;
    }
    return g;
}

TaylorExpression
implicit(const ScalarFunctionInterface& f, const TaylorFunction& g)
{
    return TaylorExpression(g.domain(),implicit(f,g.models()));
}

TaylorFunction
implicit(const TaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(f.argument_size()>f.result_size(),"f.argument_size()<=f.result_size() in implicit(f): f="<<f);
    uint fas=f.argument_size();
    uint has=f.argument_size()-f.result_size();
    Vector<Interval> hdom=project(f.domain(),range(0,has));
    Vector<Interval> hcodom=project(f.domain(),range(has,fas));
    return TaylorFunction(hdom,scale(implicit(f.models()),hcodom));
}

TaylorFunction
flow(const FunctionInterface& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    TaylorFunction phi0=embed(TaylorFunction::identity(d),Vector<Interval>(1u,Interval(-h,+h)));
    TaylorFunction phi=phi0;

    for(uint i=0; i!=10; ++i) {
        phi=antiderivative(compose(vf,phi),vf.result_size());
    }

    return phi;
}

TaylorFunction
flow(const TaylorFunction& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    return flow(vf,d,Interval(-h,+h),o);
}

TaylorFunction
unchecked_flow(const TaylorFunction& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    return unchecked_flow(vf,d,Interval(-h,+h),o);
}

TaylorFunction
flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& h, const uint o)
{
    ARIADNE_ASSERT(subset(d,vf.domain()));
    Vector<Interval> vf_range=vf.range();
    Vector<Interval> euler_step=d+h*vf_range;
    if(!subset(euler_step,vf.domain())) {
        ARIADNE_THROW(FlowBoundsException,"flow(TaylorFunction,Box,Interval,Nat)",
                      std::setprecision(20)<<"Euler step "<<euler_step<<" of vector field "<<vf<<
                      " with range "<<vf_range<<" starting in domain "<<d<<
                      " over time interval "<<h<<" does not remain in domain of vector field.");
    }

    return unchecked_flow(vf,d,h,o);
}


template<class X> inline Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const X& s3) {
    Vector<X> r(v1.size()+v2.size()+1u);
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    r[v1.size()+v2.size()]=s3;
    return r;
}


TaylorFunction
parameterised_flow(const TaylorFunction& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    ARIADNE_ASSERT(vf.result_size()<=vf.argument_size());
    ARIADNE_ASSERT(d.size()==vf.result_size());
    ARIADNE_ASSERT(h>0.0);

    uint np=vf.argument_size()-vf.result_size();
    uint nx=vf.result_size();

    Vector<Interval> const bp=project(vf.domain(),range(0,np)); // The bounding box of the parameters
    Vector<Interval> const bx=project(vf.domain(),range(np,np+nx)); // The bounding box of the variables
    Vector<Interval> const& dp=bp;  // The domain of the parameters
    Vector<Interval> const& dx=d;  // The domain of the variables

    std::cerr<<"dx="<<dx<<"\nbx="<<bx<<"\ndp="<<dp<<"\n";
    ARIADNE_ASSERT(subset(dx,bx));
    Vector<Interval> vf_range=vf.range();
    Vector<Interval> euler_step=dx+Interval(0,+h)*vf_range;
    if(!subset(euler_step,bx)) {
        ARIADNE_THROW(FlowBoundsException,"parameterised_flow(TaylorFunction,Box,Float,Nat)",
                      std::setprecision(20)<<"Euler step "<<euler_step<<" of vector field "<<vf<<
                      " with range "<<vf_range<<" starting in domain "<<d<<
                      " over time interval "<<h<<" does not remain in domain of vector field.");
    }

    // Sanity check that vector field domain has nonempty interior
    for(uint i=0; i!=nx; ++i) { ARIADNE_ASSERT_MSG(bx[i].radius()>0.0,"Domain of vector field "<<bx<<" has non-empty interior."); }
    for(uint i=0; i!=nx; ++i) { if(dx[i].radius()<=0) { std::cerr<<"Initial set "<<dx<<" has non-empty interior.\n";  } }
    for(uint i=0; i!=np; ++i) { if(dp[i].radius()<=0) { std::cerr<<"Parameter set "<<dp<<" has non-empty interior.\n";  } }

    // Scale multiply models by inverse reciprocal of radius
    Vector<TaylorModel> unit_scaled_vf(nx);
    for(uint i=0; i!=nx; ++i) { unit_scaled_vf[i]=vf.models()[i]*(h/rad_ivl(bx[i])); }
    std::cerr<<"unit_scaled_vf="<<unit_scaled_vf<<std::endl;

    Vector<TaylorModel> y0=embed(np,TaylorModel::scalings(d));
    std::cerr<<"y0="<<y0<<std::endl;
    Vector<TaylorModel> unit_scaled_y0=unscale(y0,bx);
    std::cerr<<"unit_scaled_y0"<<unit_scaled_y0<<std::endl;

    Vector<TaylorModel> unit_scaled_flow=parameterised_flow(unit_scaled_vf,unit_scaled_y0,o);
    std::cerr<<"unit_scaled_flow="<<unit_scaled_flow<<std::endl;

    Vector<TaylorModel> model_flow=scale(unit_scaled_flow,bx);
    std::cerr<<"model_flow="<<model_flow<<std::endl;

    // Check if flow is only positive
    if(0.0==0) { model_flow=split(model_flow,np+nx,true); }
    std::cerr<<"split_model_flow="<<model_flow<<std::endl;

    TaylorFunction flow(join(dp,dx,Interval(0,+h)),model_flow);
    std::cerr<<"\nflow="<<flow<<"\n"<<std::endl;

    return flow;
}

// This method should be used if we know already that the flow over time h remains in
// the domain of the vector field approximation, for example, if this has been
// checked for the original flow
TaylorFunction
unchecked_flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& h, const uint o)
{
    uint n=vf.size();
    Float hmag=mag(h);
    const Vector<Interval>& b=vf.domain();

    assert(h.lower()==-h.upper() || h.lower()==0);

    // Sanity check that vector field domain has nonempty interior
    for(uint i=0; i!=n; ++i) { ARIADNE_ASSERT_MSG(b[i].radius()>0.0,"Domain of vector field "<<b<<" has non-empty interior."); }
    for(uint i=0; i!=n; ++i) { if(d[i].radius()<=0) { std::cerr<<"Initial set "<<d<<" has non-empty interior.\n";  } }


    // Scale multiply models by inverse reciprocal of radius
    Vector<TaylorModel> unit_scaled_vf(vf.size());
    for(uint i=0; i!=vf.size(); ++i) { unit_scaled_vf[i]=vf.models()[i]*(hmag/rad_ivl(b[i])); }
    //std::cerr<<"unit_scaled_vf="<<unit_scaled_vf<<std::endl;

    Vector<TaylorModel> y0=TaylorModel::scalings(d);
    //std::cerr<<"y0="<<y0<<std::endl;
    Vector<TaylorModel> unit_scaled_y0=embed(unscale(y0,b),1u);
    //std::cerr<<"unit_scaled_y0"<<unit_scaled_y0<<std::endl;

    Vector<TaylorModel> unit_scaled_flow=unchecked_flow(unit_scaled_vf,unit_scaled_y0,o);
    //std::cerr<<"unit_scaled_flow="<<unit_scaled_flow<<std::endl;

    Vector<TaylorModel> model_flow=scale(unit_scaled_flow,b);
    //std::cerr<<"model_flow="<<model_flow<<std::endl;

    if(h.lower()==0) { model_flow=split(model_flow,d.size(),true); }

    TaylorFunction flow(join(d,h),model_flow);
    //std::cerr<<"\nflow="<<flow<<"\n"<<std::endl;

    return flow;
}






std::ostream&
TaylorFunction::write(std::ostream& os) const
{
    //os << "TaylorFunction( "<<this->domain()<<" , ";
    //for(uint i=0; i!=this->result_size(); ++i) {
    //    os << (i==0?'[':',')<<this->_models[i].expansion()<<","<<this->_models[i].error();
    //}
    //return os << "] )";
    return os << "TaylorFunction(d=" << this->domain() << ", p=" << this->polynomial() << " m=" << this->models() << ")";
}



std::ostream&
operator<<(std::ostream& os, const TaylorFunction& p)
{
    return p.write(os);
}



/*
latexstream&
operator<<(Output::latexstream& texs, const TaylorFunction& p)
{
    using namespace Function;
    texs << "%TaylorFunction\n";
    texs << "\\ensuremath{\n";
    texs << "\\left( \\begin{array}{c}\n";
    char var='x';
    for(uint i=0; i!=p.result_size(); ++i) {
        bool first = true;
        if(i!=0) { texs << "\\\\"; }
        for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
            const Interval& a=p.centre_derivatives()[i][j];
            if(a!=0) {
                if(first) { first=false; }
                else { if(a>0) { texs << '+'; } }
                if(a==1) { if(j.degree()==0) { texs << a; } }
                else if(a==-1) { if(j.degree()==0) { texs << a; } else { texs << '-'; } }
                else { texs << a << ' '; }
                for(uint k=0; k!=p.argument_size(); ++k) {
                    if(j[k]!=0) {
                        texs << var << "_{ " << k+1 << "}";
                        if(j[k]!=1) {
                            texs << "^{" << j[k] << "}";
                        }
                    }
                    texs << " ";
                }
            }
        }
        texs << "\n";
    }
    texs << "\\end{array}\\right)\n}\n";
    return texs;
}
*/


} // namespace Ariadne
