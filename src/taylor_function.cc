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
 
#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "sparse_differential.h"
#include "differential_vector.h"
#include "function_interface.h"
#include "taylor_variable.h"
#include "taylor_function.h"

namespace Ariadne {

typedef unsigned int uint;
    
typedef Vector<Float> Point;
typedef Vector<Interval> Box;

TaylorFunction::TaylorFunction() 
    : _domain(), 
      _expansion()
{
}


TaylorFunction::TaylorFunction(uint rs, uint as) 
    : _domain(as,I(-1,1)),
      _expansion(rs,as)
{
}


TaylorFunction::TaylorFunction(uint rs, uint as, ushort o, ushort s) 
    : _domain(as,I(-1,1)),
      _expansion(rs,as)
{
}


TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector<TaylorVariable>& e)
    : _domain(d),
      _expansion(e)
{
    for(uint i=0; i!=e.size(); ++i) {
        ARIADNE_ASSERT(d.size()==e[i].argument_size());
    }
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const FunctionInterface& f)
    : _domain(d),
      _expansion(f.result_size())
{
    ARIADNE_ASSERT(d.size()==f.argument_size());
    
    
    Vector<TaylorVariable> x=TaylorVariable::variables(Vector<Float>(f.argument_size(),0.0));
    for(uint i=0; i!=x.size(); ++i) {
        const Interval& di=d[i];
        Interval dm=add_ivl(di.l/2,di.u/2);
        Interval dr=sub_ivl(di.u/2,di.l/2);
        x[i]*=dr;
        x[i]+=dm;
    }

    this->_expansion=f.evaluate(x);
}





TaylorFunction::TaylorFunction(const Vector<Interval>& d, 
                         ushort o, ushort s,
                         const FunctionInterface& f)
    : _domain(d),
      _expansion(f.result_size())
{
    ARIADNE_NOT_IMPLEMENTED;
}




/*
TaylorFunction
TaylorFunction::zero(uint rs, uint as)  
{
    return TaylorFunction(Vector<Interval>(as,Interval(-inf(),+inf()),TaylorVariable::constants(rs,as,Vector<Float>(rs,0.0))));
}
*/




TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Float>& c)  
{
    return TaylorFunction(d,TaylorVariable::constants(d.size(),c));
}

TaylorFunction
TaylorFunction::identity(const Vector<Interval>& d)  
{
    return TaylorFunction(d,TaylorVariable::variables(midpoint(d)));
}



bool
TaylorFunction::operator==(const TaylorFunction& tm) const
{
    return this->_domain==tm._domain
        && this->_expansion==tm._expansion;
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
        result[i]=this->_expansion[i].range();
    }
    return result;
}



const Vector<TaylorVariable>&
TaylorFunction::variables() const
{ 
    return this->_expansion;
}


const Vector<TaylorVariable>&
TaylorFunction::expansion() const
{ 
    return this->_expansion;
}




uint 
TaylorFunction::argument_size() const
{ 
    return this->_expansion[0].argument_size(); 
}


uint 
TaylorFunction::result_size() const 
{ 
    return this->_expansion.size();
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
    if(this->argument_size()!=x.size()) {
        ARIADNE_THROW(std::runtime_error,"TaylorFunction::evaluate(Vector)","Incompatible argument size");
    }

    //TODO: Make this MUCH more efficient!

    // Scale x to domain
    Vector<Interval> scaled_x(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        const Interval& d=this->_domain[i];
        Interval dm=add_ivl(d.l/2,d.u/2);
        Interval dr=sub_ivl(d.u/2,d.l/2);
        scaled_x[i]=(x[i]-dm)/dr;
    }

    // Compute result
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_expansion[i].evaluate(scaled_x);
    }
    return result;
}


Matrix<Float> 
TaylorFunction::jacobian(const Vector<Float>& x) const
{
    DifferentialVector< SparseDifferential<Float> > y(this->argument_size(),SparseDifferential<Float>(this->argument_size(),1u));
    for(uint j=0; j!=this->argument_size(); ++j) {
        y[j].set_value(this->domain()[j].midpoint());
        y[j].set_gradient(j,this->domain()[j].radius());
    }
    DifferentialVector< SparseDifferential<Float> > t(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        t[i]=SparseDifferential<Float>(this->_expansion[i].expansion());
    }
    return get_jacobian(compose(t-x,y));
}


TaylorFunction
restrict(const TaylorFunction& tm, const Box& nd)
{
    const Vector<Interval>& od=tm.domain();
    ARIADNE_ASSERT(subset(nd,od));
    if(nd==od) { return tm; }
    Vector<TaylorVariable> sm=Vector<TaylorVariable>(nd);
    return TaylorFunction(nd,compose(tm._expansion,tm._domain,sm));
}



TaylorFunction
operator+(const TaylorFunction& p1, const TaylorFunction& p2)
{
    ARIADNE_ASSERT(!intersection(p1.domain(),p2.domain()).empty());
    if(p1.domain()==p2.domain()) {
        return TaylorFunction(p1._domain,Vector<TaylorVariable>(p1._expansion+p2._expansion));
    } else {
        Box new_domain=intersection(p1.domain(),p2.domain());
        return operator+(restrict(p1,new_domain),restrict(p2,new_domain));
    }
}


TaylorFunction
operator-(const TaylorFunction& p1, const TaylorFunction& p2)
{
    ARIADNE_ASSERT(!intersection(p1.domain(),p2.domain()).empty());
    if(p1.domain()==p2.domain()) {
        return TaylorFunction(p1._domain,Vector<TaylorVariable>(p1._expansion-p2._expansion));
    } else {
        ARIADNE_NOT_IMPLEMENTED;
        Box new_domain=intersection(p1.domain(),p2.domain());
        return operator-(restrict(p1,new_domain),restrict(p2,new_domain));
    }
}



TaylorFunction
operator+(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f._domain,f._expansion+c);
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f._domain,f._expansion-c);
}


/*

TaylorFunction
mul(const TaylorFunction& tm, const Float& x)
{
    return TaylorFunction(tm.domain(),tm.centre(),tm.smoothness(),tm.centre_derivatives()*x);
}

*/


/*

void
mul(TaylorFunction& p0, const TaylorFunction& p1, const TaylorFunction& p2)
{
    if(p1.result_size()!=1u) {
        ARIADNE_THROW(IncompatibleSizes,"mul(TaylorFunction,TaylorFunction)","p1.result_size()="<<p1.result_size());
    }
    if(p2.result_size()!=1u) {
        ARIADNE_THROW(IncompatibleSizes,"mul(TaylorFunction,TaylorFunction)","p2.result_size()="<<p2.result_size());
    }
    if(p1.argument_size()!=p2.argument_size()) {
        ARIADNE_THROW(IncompatibleSizes,"add(TaylorFunction p1,TaylorFunction p2)","p1.result_size()="<<p1.result_size()<<", p2.result_size()="<<p2.result_size());
    }

    uint rs=1u;
    uint as=p1.argument_size();
    uint d1=p1.order();
    uint d2=p2.order();
    uint d=d1+d2;
    uint s1=p1.smoothness();
    uint s2=p2.smoothness();
    uint s=std::max(s1,s2);

    p0.resize(rs,as,d,s);
    MultiIndex j0(as);
    for(MultiIndex j1(as); j1.degree()<=d1; ++j1) {
        const Float& x1=p1.get(0u,j1);
        for(MultiIndex j2(as); j2.degree()<=d2; ++j2) {
            const R& x2=p2.get(0u,j2);
            j0=j1+j2;
            p0.at(0u,j0)+=x1*x2;
        }
    }
}


*/





/*
TaylorFunction&
operator*=(TaylorFunction& p0, const Interval& x1)
{
    for(uint i=0; i!=p0._data.size(); ++i) {
        p0._data[i]*=x1;
    }
    return p0;
}
*/

TaylorFunction
compose(const TaylorFunction& p1, const TaylorFunction& p2)
{
    ARIADNE_ASSERT(subset(p2.range(),p1.domain()));
    return TaylorFunction(p2.domain(),Ariadne::compose(p1.expansion(),p1.domain(),p2.expansion()));
}



TaylorFunction
antiderivative(const TaylorFunction& tm, uint k) 
{
    return TaylorFunction(tm.domain(),antiderivative(tm.expansion(),tm.domain()[k],k));
}



/*
void
compose(TaylorFunction& p0, const TaylorFunction& p1, const TaylorFunction& p2)
{
    // TODO: Improve this algorithm as it's critical!!
    if(p1.argument_size()!=p2.result_size()) {
        ARIADNE_THROW(IncompatibleSizes,"compose(TaylorFunction p1,TaylorFunction p2)","p1.argument_size()="<<p1.argument_size()<<", p2.result_size()="<<p2.result_size());
    }

    if(p1.order()==0) { 
        p0=static_cast< TaylorFunction >(p1); 
        return;
    }

    p0.resize(p1.result_size(),p2.argument_size(),p1.order()*p2.order(),std::max(p1.smoothness(),p2.smoothness()));
  
    TaylorFunction* all_powers=new TaylorFunction[p2.result_size()*(p1.order()+1)];
    TaylorFunction* powers[p2.result_size()];
    for(uint i=0; i!=p2.result_size(); ++i) {
        powers[i]=all_powers+i*(p1.order()+1);
    }

    for(uint i=0; i!=p2.result_size(); ++i) {
        powers[i][0]=TaylorFunction::one(p2.argument_size());
        powers[i][1]=p2.component(i);
        if(p1.order()>=2) {
            powers[i][2]=pow(powers[i][1],2);
        }
        for(uint j=3; j<=p1.order(); ++j) {
            powers[i][j]=powers[i][2]*powers[i][j-2];
        }
    }
  
    TaylorFunction* results=new TaylorFunction[p1.result_size()];
    for(uint i=0; i!=p1.result_size(); ++i) {
        results[i]=TaylorFunction::zero(1u,p2.argument_size());
    }

    for(uint i=0; i!=p1.result_size(); ++i) {
        for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
            TaylorFunction t=TaylorFunction::constant(p2.argument_size(),p1.get(i,j));
            for(uint k=0; k!=p1.argument_size(); ++k) {
                t=t*powers[k][j[k]];
            }
            results[i]=results[i]+t;
        }
    }
  
    for(uint i=0; i!=p0.result_size(); ++i) {
        for(uint j=0; j!=p0.data().size()/p0.result_size(); ++j) {
            p0._data[i+j*p0.result_size()]=results[i].data()[j];
        }
    }
  
    delete[] results;
    delete[] all_powers;
}

*/

/*

void
derivative(TaylorFunction& p0, const TaylorFunction& p1, uint k)
{
    if(p1.smoothness()==0) {
        ARIADNE_THROW(std::runtime_error,"derivative(TaylorFunction,uint)"," model has smoothness 0");
    }

    p0.resize(p1.result_size(),p1.argument_size(),p1.order()-1,p1.smoothness()-1);
  
    MultiIndex dj(p1.argument_size());

    for(uint i=0; i!=p0.result_size(); ++i) {
        for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
            if(j[k]!=0) {
                dj=j;
                dj.decrement_index(k);
                p0.at(i,dj)+=static_cast<int>(j[k])*p1.get(i,j);
            }
        }
    }
}

*/

 /*
Matrix<Interval> 
TaylorFunction::jacobian(const Vector<Float>& s) const
{
    return this->jacobian(Vector<Interval>(s));
}


Matrix<Interval> 
TaylorFunction::jacobian(const Vector<Interval>& x) const
{
    Matrix<Interval> J(this->result_size(),this->argument_size());
    Box w=x-this->centre();
    array< array<I> > powers=this->_powers(w);

    for(uint j=0; j!=this->argument_size(); ++j) {
        for(MultiIndex m(this->argument_size()); m.degree()<this->order(); ++m) {
            MultiIndex n=m;
            int c=n[j];
            if(c!=0) {
                n.decrement_index(j);
                I a=c;
                for(uint k=0; k!=this->argument_size(); ++k) {
                    a*=powers[k][n[k]];
                }
                for(uint i=0; i!=this->result_size(); ++i) {
                    I amim=a*this->_centre_derivatives[i][m];
                    J[i][j]+=amim;
                }
            }
        }
        for(MultiIndex m=MultiIndex::first(this->argument_size(),this->order()); m.degree()<=this->order(); ++m) {
            MultiIndex n=m;
            int c=n[j];
            if(c!=0) {
                n.decrement_index(j);
                I a=c;
                for(uint k=0; k!=this->argument_size(); ++k) {
                    a*=powers[k][n[k]];
                }
                for(uint i=0; i!=this->result_size(); ++i) {
                    J[i][j]+=a*this->_domain_derivatives[i][m];
                }
            }
        }
    }
  
    return J;
}
 */



 
TaylorFunction 
inverse(const TaylorFunction& p, const Point& v)
{
    assert(p.result_size()==p.argument_size());
    assert(p.argument_size()==v.size());

    // The following are only to simplfy testing.
    assert(v==Point(v.size(),0));
    assert(p.evaluate(v)==Point(p.result_size(),0));
    typedef Interval I;

    Point c=midpoint(Box(p.evaluate(v)));
    Matrix<Float> J=p.jacobian(v);
  
    Matrix<Float> invJ=inverse(J);

    // FIXME: Need to re-solve for image of centre. What should initial set be? Different code needed for Rational?
    Box invf=v;

    // FIXME: Give correct initial conditions
    TaylorFunction result;

    for(MultiIndex m(p.result_size()); m.degree()<=2; ++m) {
        if(m.degree()==0) {
            for(uint i=0; i!=p.argument_size(); ++i) {
                result._expansion[i][m]=v[i];
            }
        } else if(m.degree()==1) {
            for(uint i=0; i!=p.argument_size(); ++i) {
                result._expansion[i][m]=invJ[i][m.position()-1];
            }
        } else {
            // FIXME: Add code for higher indices
        }
    }
    return result;
}


 

TaylorFunction 
implicit(const TaylorFunction& p)
{
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorFunction 
flow(const TaylorFunction& p, const Vector<Interval>& domain, const Interval& time)
{
    return TaylorFunction(join(domain,time),flow(p.expansion(),domain,time,p.domain()));
}


array< array<Interval> >
TaylorFunction::_powers(const Vector<Interval>& v) const
{
    uint order=21;
    array< array<I> > powers(this->argument_size(), array<I>(order));
    for(uint i=0; i!=this->argument_size(); ++i) {
        powers[i][0]=1;
        if(order>=1) {
            powers[i][1]=v(i);
            if(order>=2) {
                powers[i][2]=pow(v(i),2);
                for(uint j=3; j<=order; ++j) {
                    powers[i][j]=powers[i][2]*powers[i][j-2];
                }
            }
        }
    }
    return powers;
}



std::ostream&
TaylorFunction::write(std::ostream& os) const 
{
    os << "TaylorFunction( "<<this->_domain<<" , ";
    for(uint i=0; i!=this->result_size(); ++i) { 
        os << (i==0?'[':',')<<this->_expansion[i].expansion()<<","<<this->_expansion[i].error();
    }
    return os << "] )";
    /*
    os << "TaylorFunction(\n";
    os << "  domain=" << this->domain() << ",\n" << std::flush;
    os << "  range=" << this->range() << ",\n" << std::flush;
    os << "  model=" << this->_expansion << "\n" << std::flush;
    //os << "  series=" << Ariadne::expansion(this->_expansion,this->_domain) << "\n" << std::flush;
    os << ")\n";
    return os;
    */
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
