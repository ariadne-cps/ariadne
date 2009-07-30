/***************************************************************************
 *            taylor_expression.cc
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include <iomanip>

#include "config.h"
#include "rounding.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "expansion.h"
#include "polynomial.h"
#include "differential.h"
#include "taylor_model.h"
#include "taylor_expression.h"
#include "taylor_function.h"
#include "expression_interface.h"
#include "exceptions.h"

namespace Ariadne {


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

TaylorExpression::TaylorExpression(const DomainType& d, const ExpressionInterface& f)
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



TaylorExpression restrict(const TaylorExpression& tv, const Vector<Interval>& d) {
    ARIADNE_ASSERT(subset(d,tv.domain()));
    if(d==tv.domain()) { return tv; }
    Vector<TaylorModel> s=TaylorModel::rescalings(tv.domain(),d);
    return TaylorExpression(d,compose(tv._model,s));
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
compose(const ExpressionInterface& f, const Vector<TaylorExpression>& g)
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
implicit(const ExpressionInterface& f, const Vector<TaylorExpression>& g)
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
implicit(const ExpressionInterface& f, const Vector<Interval>& d)
{
    ARIADNE_ASSERT(f.argument_size()>=1u);
    ARIADNE_ASSERT(d.size()+1u==f.argument_size());

    Vector<TaylorExpression> id=TaylorExpression::variables(d);
    //for(uint i=0; i!=d.size(); ++i) { id[i].model().set_accuracy(a); }
    return implicit(f,id);
}


TaylorExpression
crossing_time(const ExpressionInterface& g, const TaylorFunction& phi)
{
    Vector<Interval> d=project(phi.domain(),range(0,phi.result_size()));
    Interval h=phi.domain()[phi.result_size()];

    TaylorExpression gphi(phi.domain(),g.evaluate(phi.models()));

    return implicit(gphi);
}


std::ostream&
operator<<(std::ostream& os, const TaylorExpression& tv) {
    return os << "TaylorExpression(" << tv.domain() << "," << tv.expansion() << "," << tv.error() << ")";
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





} //namespace Ariadne


