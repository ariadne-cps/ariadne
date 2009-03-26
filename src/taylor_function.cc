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
#include "taylor_variable.h"
#include "taylor_function.h"

namespace Ariadne {

typedef unsigned int uint;

typedef Vector<Float> Point;
typedef Vector<Interval> Box;

TaylorFunction::TaylorFunction()
    : _models()
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
        _models[i]=TaylorModel(f[i]);
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
    Vector<TaylorModel> x=TaylorModel::variables(d);
    this->_models=f.evaluate(x);
}


TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Float> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<TaylorModel> x=TaylorModel::variables(d);
    this->_models=Ariadne::evaluate(p,x);
}










TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Float>& c)
{
    return TaylorFunction(d,TaylorModel::constants(d.size(),c));
}

TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Interval>& r)
{
    Vector<TaylorModel> e(r.size(),TaylorModel(d.size()));
    for(uint i=0; i!=r.size(); ++i) {
        e[i]=r[i];
    }
    return TaylorFunction(d,e);
}

TaylorFunction
TaylorFunction::identity(const Vector<Interval>& d)
{
    return TaylorFunction(TaylorModel::variables(d));
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



const Vector<Interval>&
TaylorFunction::domain() const
{
    return this->_domain;
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
    return this->_models[0].argument_size();
}


uint
TaylorFunction::result_size() const
{
    return this->_models.size();
}






















TaylorFunction
TaylorFunction::truncate(ushort degree) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



Vector<Interval>
TaylorFunction::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}


Vector<Interval>
TaylorFunction::evaluate(const Vector<Interval>& x) const
{
    return Ariadne::evaluate(this->_models,x);

    if(this->argument_size()!=x.size()) {
        ARIADNE_THROW(std::runtime_error,"TaylorFunction::evaluate(Vector)","Incompatible argument size");
    }

    //TODO: Make this MUCH more efficient!

    // Scale x to domain
    Vector<Interval> scaled_x(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        const Interval& d=this->domain()[i];
        Interval dm=add_ivl(d.l/2,d.u/2);
        Interval dr=sub_ivl(d.u/2,d.l/2);
        scaled_x[i]=(x[i]-dm)/dr;
    }

    // Compute result
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=Ariadne::evaluate(this->_models[i],scaled_x);
    }
    return result;
}


Matrix<Interval>
TaylorFunction::jacobian(const Vector<Interval>& x) const
{
    return Ariadne::jacobian(this->_models,x);
}


TaylorFunction
join(const TaylorFunction& f, const TaylorFunction& g)
{
    ARIADNE_ASSERT(f.domain()==g.domain());
    return TaylorFunction(join(f.models(),g.models()));
}

TaylorFunction
restrict(const TaylorFunction& f, const Vector<Interval>& d)
{
    ARIADNE_ASSERT(subset(d,f.domain()));
    if(d==f.domain()) { return f; }
    Vector<TaylorModel> s=TaylorModel::variables(d);
    return TaylorFunction(d,compose(f._models,f._domain,s));
}



TaylorFunction&
operator+=(TaylorFunction& f, const Vector<Interval>& e)
{
    ARIADNE_ASSERT(f.result_size()==e.size());
    f._models+=e;
    return f;
}


TaylorFunction
operator+(const TaylorFunction& f1, const TaylorFunction& f2)
{
    ARIADNE_ASSERT(!intersection(f1.domain(),f2.domain()).empty());
    if(f1.domain()==f2.domain()) {
        return TaylorFunction(Vector<TaylorModel>(f1.models()+f2.models()));
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
        return TaylorFunction(Vector<TaylorModel>(f1.models()+f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}



TaylorFunction
operator+(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f.models()+c);
}

TaylorFunction
operator+(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(f.models()+c);
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f.models()-c);
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(f.models()-c);
}

TaylorFunction
operator*(const Matrix<Interval>& A, const TaylorFunction& f)
{
    return TaylorFunction(prod(A,f.models()));
}





TaylorFunction
combine(const TaylorFunction& f1, const TaylorFunction& f2)
{
    return TaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
}


TaylorFunction
compose(const TaylorFunction& g, const TaylorFunction& f)
{
    if(!subset(f.range(),g.domain())) {
        std::cerr<<"f.range()="<<f.range()<<" is not a subset of g.domain()="<<g.domain()<<std::endl;
        ARIADNE_ASSERT(subset(f.range(),g.domain()));
    }
    return TaylorFunction(f.domain(),Ariadne::compose(g.models(),g.domain(),f.models()));
}



TaylorFunction
antiderivative(const TaylorFunction& f, uint k)
{
    return TaylorFunction(f.domain(),antiderivative(f.models(),f.domain()[k],k));
}

TaylorFunction
implicit(const TaylorFunction& f)
{
    uint fas=f.argument_size();
    uint has=f.argument_size()-f.result_size();
    Vector<Interval> hdom=project(f.domain(),range(0,has));
    Vector<Interval> hcodom=project(f.domain(),range(has,fas));
    return TaylorFunction(hdom,implicit(f.models(),hcodom));
}

TaylorFunction
flow(const TaylorFunction& f, const Vector<Interval>& d, const Interval& h, const uint o)
{
    return TaylorFunction(join(d,h),flow(f.models(),d,h,o));
}



bool
refines(const TaylorFunction& f, const TaylorFunction& g)
{
    return refines(f.models(),g.models());
}



std::ostream&
TaylorFunction::write(std::ostream& os) const
{
    os << "TaylorFunction( "<<this->domain()<<" , ";
    for(uint i=0; i!=this->result_size(); ++i) {
        os << (i==0?'[':',')<<this->_models[i].expansion()<<","<<this->_models[i].error();
    }
    return os << "] )";
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



TaylorExpression::TaylorExpression()
    : _domain(), _model() { }

TaylorExpression::TaylorExpression(const Vector<Interval>& domain)
    : _domain(domain), _model(domain.size()) { }

TaylorExpression::TaylorExpression(const Vector<Interval>& domain, const TaylorModel& model)
    : _domain(domain), _model(model) { }

TaylorExpression TaylorExpression::variable(const Vector<Interval>& domain, Nat i) {
    TaylorExpression t(domain); t._model=TaylorModel::variables(domain)[i]; return t; }

Vector<TaylorExpression> TaylorExpression::variables(const Vector<Interval>& domain) {
    Vector<TaylorExpression> r(domain.size(),TaylorExpression(domain));
    for(uint i=0; i!=r.size(); ++i) { r[i]=variable(domain,i); }
    return r; }

const Vector<Interval>& TaylorExpression::domain() const {
    return  this->_domain; }
const TaylorModel& TaylorExpression::model() const {
    return  this->_model; }

TaylorExpression operator+(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c+t._model); }
TaylorExpression operator-(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c-t._model); }
TaylorExpression operator*(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c*t._model); }
TaylorExpression operator/(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c/t._model); }

TaylorExpression operator+(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c+t._model); }
TaylorExpression operator-(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c-t._model); }
TaylorExpression operator*(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c*t._model); }
TaylorExpression operator/(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c/t._model); }

TaylorExpression operator+(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model+c); }
TaylorExpression operator-(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model-c); }
TaylorExpression operator*(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model*c); }
TaylorExpression operator/(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model/c); }

TaylorExpression operator+(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model+c); }
TaylorExpression operator-(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model-c); }
TaylorExpression operator*(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model*c); }
TaylorExpression operator/(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model/c); }

TaylorExpression operator+(const TaylorExpression& t) {
    return TaylorExpression(t._domain,+t._model);  }
TaylorExpression operator-(const TaylorExpression& t) {
    return TaylorExpression(t._domain,-t._model);  }
TaylorExpression operator+(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model+t2._model); }
TaylorExpression operator-(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model-t2._model); }
TaylorExpression operator*(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model*t2._model); }
TaylorExpression operator/(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model/t2._model); }

TaylorFunction operator*(const TaylorExpression& t, const Vector<Float>& v) {
    return TaylorFunction(Vector<TaylorModel>(t._model*v)); }
TaylorFunction operator*(const Vector<Float>& v, const TaylorExpression& t) {
    return t*v; };

std::ostream& operator<<(std::ostream& os, const TaylorExpression& t) {
    return os<<"("<<t._domain<<","<<t._model<<")"; }

} // namespace Ariadne
