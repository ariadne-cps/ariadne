/***************************************************************************
 *            solver.cc
 *
 *  Copyright  2006-9  Pieter Collins
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

#include "solver.h"

#include "logging.h"
#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "taylor_model.h"
#include "function.h"
#include "taylor_function.h"

namespace {

using namespace Ariadne;

typedef unsigned int uint;

void
solve_all(Set< Vector<Interval> >& r,
          const SolverInterface& s,
          const RealVectorFunction& f,
          const Vector<Interval>& ix)
{
    uint verbosity=s.verbosity;
    ARIADNE_LOG(5,"solve_all(f,ix): f="<<f<<", ix="<<ix<<"\n");

    // Test for no solution
    const Vector<Interval> z(ix.size());
    if(disjoint(f.evaluate(ix),z)) {
        return;
    }

    bool invertible_jacobian=true;
    //Vector<Interval> nx=2*ix-midpoint(ix);
    Vector<Interval> nx=ix;
    try {
        Matrix<Interval> Jinv=inverse(f.jacobian(nx));
    }
    catch(const SingularMatrixException& e) {
        invertible_jacobian=false;
    }

    bool need_to_split=true;

    if(invertible_jacobian) {
        //std::cerr<<"Nonsingular matrix -- applying contractor\n";
        try {
            Vector<Interval> y=s.zero(f,nx);
            bool is_new=true;
            for(Set<Vector<Interval> >::const_iterator iter=r.begin(); iter!=r.end(); ++iter) {
                if(!disjoint(y,*iter)) {
                    is_new=false;
                    break;
                }
            }
            if(is_new) {
                r.insert(y);
            }
            need_to_split=false;
        }
        catch(const NoSolutionException& e) {
            //ARIADNE_WARN("NoSolutionException exception: "<<e.what()<<"\n");
            need_to_split=false;
        }
        catch(const SolverException& e) {
            ARIADNE_WARN("SolverException exception: "<<e.what()<<"\n");
            need_to_split=true;
            // No solution found, try splitting
        }
    } else {
        need_to_split=true;
    }

    if(need_to_split) {
        // If radius is too small, assume solution is not verified
        if(radius(ix)<s.maximum_error()) {
            if(!invertible_jacobian) {
                ARIADNE_WARN("Cannot verify solution in "<<ix<<" with f="<<f(ix)<<"); "
                             <<"Jacobian "<<f.jacobian(nx)<<" is not invertible; "
                             <<"approximate inverse="<<inverse(midpoint(f.jacobian(nx)))<<"\n");
            } else {
                ARIADNE_WARN("Cannot verify or falsify solution in "<<ix<<"; f("<<ix<<")="<<f(ix)<<".\n");
            }
            return;
        }

        //std::cerr<<"  Splitting "<<ix<<"\n";
        std::pair< Vector<Interval>, Vector<Interval> > splt=split(ix);
        solve_all(r,s,f,splt.first);
        solve_all(r,s,f,splt.second);
    }

}

Vector<Interval> ranges(const Vector<IntervalTaylorModel>& f) {
    Vector<Interval> r(f.size()); for(uint i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

Vector<IntervalTaylorModel>& clobber(Vector<IntervalTaylorModel>& h) {
    for(uint i=0; i!=h.size(); ++i) { h[i].set_error(0.0); } return h; }

// Compute the Jacobian over an arbitrary domain
Matrix<Interval>
jacobian2(const Vector<IntervalTaylorModel>& f, const Vector<Interval>& x)
{
    Vector< Differential<Interval> > dx(x.size());
    for(uint i=0; i!=x.size()-f.size(); ++i) {
        dx[i]=Differential<Interval>::constant(f.size(),1u,x[i]); }
    for(uint i=0; i!=f.size(); ++i) {
        uint j=i+(x.size()-f.size());
        dx[j]=Differential<Interval>::variable(f.size(),1u,x[j],i); }
    Vector< Differential<Interval> > df(f.size());
    for(uint i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<Interval> J=jacobian(df);
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<Float>
jacobian2_value(const Vector<IntervalTaylorModel>& f)
{
    const uint rs=f.size();
    const uint fas=f[0].argument_size();
    const uint has=fas-rs;
    Matrix<Float> J(rs,rs);
    MultiIndex a(fas);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=rs; ++j) {
            a[has+j]=1; const Float x=f[i][a]; J[i][j]=x; a[has+j]=0;
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<Interval>
jacobian2_range(const Vector<IntervalTaylorModel>& f)
{
    uint rs=f.size();
    uint fas=f[0].argument_size();
    uint has=fas-rs;
    Matrix<Interval> J(rs,rs);
    for(uint i=0; i!=rs; ++i) {
        for(IntervalTaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=rs; ++k) {
                const uint c=iter->key()[has+k];
                if(c>0) {
                    const Float& x=iter->data();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=Interval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}

//Compute the implicit function by preconditioning f by the inverse
// of the Jacobian value matrix and using a Gauss-Seidel iteration scheme
// on the system y=g(y)
Vector<IntervalTaylorModel> _implicit5(const Vector<IntervalTaylorModel>& f, uint n)
{
    //std::cerr<<__FUNCTION__<<std::endl;
    uint rs=f.size(); uint fas=f[0].argument_size(); uint has=fas-rs;

    Vector<Interval> domain_h(rs,Interval(-1,+1));
    Vector<IntervalTaylorModel> id=IntervalTaylorModel::variables(has);
    Vector<IntervalTaylorModel> h=IntervalTaylorModel::constants(has,domain_h);
    Vector<IntervalTaylorModel> idh=join(id,h);

    // Compute the Jacobian of f with respect to the second arguments at the centre of the domain
    Matrix<Float> D2f=jacobian2_value(f);

    // Compute g=-D2finv*(f-D2f*y) = -D2finv*f-y
    Vector<IntervalTaylorModel> g=-f;
    for(uint i=0; i!=rs; ++i) {
        for(uint j=has; j!=fas; ++j) {
            g[i][MultiIndex::unit(fas,j)]=0;
        }
    }
    Matrix<Float> J=inverse(D2f);
    g=prod(J,g);
    for(uint i=0; i!=rs; ++i) {
        g[i].clean();
    }

    // Iterate h'=h(g(x,h(x)))
    Vector<IntervalTaylorModel> h_new;
    uint number_non_contracting=rs;
    array<bool> contracting(rs,false);
    for(uint k=0; k!=n; ++k) {
        h_new=compose(g,join(id,h));
        for(uint i=0; i!=rs; ++i) {
            if(!contracting[i]) {
                if(refines(h_new,h)) {
                    contracting[i]=true;
                    --number_non_contracting;
                }
                if(disjoint(h_new,h)) {
                    ARIADNE_THROW(ImplicitFunctionException,"implicit(Vector<IntervalTaylorModel> f)",
                                "(with f="<<f<<"): Application of Newton solver to "<<h<<" yields "<<h_new<<
                                " which is disjoint. No solution");
                }
            }
        }
        for(uint i=0; i!=rs; ++i) {
            h[i]=intersection(h[i],h_new[i]);
        }
    }

    if(number_non_contracting==0) {
        return h;
    } else {
        ARIADNE_THROW(ImplicitFunctionException,"implicit(Vector<IntervalTaylorModel>)",
                      "Could not verify solution of implicit function to "<<h<<"\n");
    }
}

Vector<IntervalTaylorModel>
newton_implicit(const Vector<IntervalTaylorModel>& f)
{
    // Check that the arguments are suitable
    ARIADNE_ASSERT(f.size()>0);
    for(uint i=1; i!=f.size(); ++i) { ARIADNE_ASSERT(f[i].argument_size()==f[0].argument_size()); }

    // Set some useful size constants
    const uint rs=f.size();
    const uint fas=f[0].argument_size();
    const uint has=fas-rs;

    // Check to see if a solution exists
    Matrix<Interval> D2f=jacobian2_range(f);
    Matrix<Interval> D2finv;
    try {
        D2finv=inverse(D2f);
    }
    catch(...) {
        ARIADNE_THROW(ImplicitFunctionException,
                      "implicit(Vector<IntervalTaylorModel>)",
                      "Jacobian "<<D2f<<" is not invertible");
    }

    uint number_of_steps=6;
    Vector<IntervalTaylorModel> id=IntervalTaylorModel::variables(has);
    Vector<IntervalTaylorModel> h=_implicit5(f,number_of_steps);

    // Perform proper Newton step improvements
    Vector<Interval> domain_h(h[0].argument_size(),Interval(-1,+1));
    for(uint i=0; i!=3; ++i) {
        D2finv=inverse(jacobian2(f,join(domain_h,ranges(h))));
        clobber(h);
        Vector<IntervalTaylorModel> fidh=compose(f,join(id,h));
        Vector<IntervalTaylorModel> dh=prod(D2finv,fidh);
        h-=dh;
    }

    // Check that the result has the correct sizes.
    ARIADNE_ASSERT(h.size()==f.size());
    for(uint i=0; i!=h.size(); ++i) {
        ARIADNE_ASSERT(h[0].argument_size()+f.size()==f[i].argument_size());
    }

    return h;
}

VectorTaylorFunction
newton_implicit(const VectorTaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(f.argument_size()>f.result_size(),"f.argument_size()<=f.result_size() in implicit(f): f="<<f);
    uint fas=f.argument_size();
    uint has=f.argument_size()-f.result_size();
    Vector<Interval> hdom=project(f.domain(),range(0,has));
    Vector<Interval> hcodom=project(f.domain(),range(has,fas));
    return VectorTaylorFunction(hdom,scale(newton_implicit(f.models()),hcodom));
}



} // namespace



namespace Ariadne {

VectorTaylorFunction evaluate(const RealVectorFunction& f,const VectorTaylorFunction& x) {
    for(uint i=0; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].domain()==x[0].domain()); }
    Vector<IntervalTaylorModel> m(x.size());
    for(uint i=0; i!=m.size(); ++i) { m[i]=x[i].model(); }
    m=f.evaluate(m);
    VectorTaylorFunction r(m.size());
    for(uint i=0; i!=r.size(); ++i) { r[i]=ScalarTaylorFunction(x[0].domain(),m[i]); }
    return r;
}

VectorTaylorFunction product(const Matrix<Float>& A,const VectorTaylorFunction& v) {
    VectorTaylorFunction r(A.row_size());
    for(uint i=0; i!=r.size(); ++i) {
        ScalarTaylorFunction t(v[0].domain());
        for(uint j=0; j!=v.size(); ++j) {
            t+=A[i][j]*v[j];
        }
        r[i]=t;
    }
    return r;
}

VectorTaylorFunction product(const Matrix<Interval>& A,const VectorTaylorFunction& v) {
    VectorTaylorFunction r(A.row_size());
    for(uint i=0; i!=r.size(); ++i) {
        ScalarTaylorFunction t(v[0].domain());
        for(uint j=0; j!=v.size(); ++j) {
            t+=A[i][j]*v[j];
        }
        r[i]=t;
    }
    return r;
}

Float radius(const VectorTaylorFunction& x) {
    Float r=0.0;
    for(uint i=0; i!=x.size(); ++i) { r=std::max(r,x[i].error()); }
    return r;
}




/*
class DifferenceFunction
    : public VectorFunctionTemplate<DifferenceFunction>
{
  public:
    DifferenceFunction(const RealVectorFunction& f) : fptr(f.clone()) { }
    virtual DifferenceFunction* clone() const { return new DifferenceFunction(*this); }
    virtual uint result_size() const { return fptr->result_size(); }
    virtual uint argument_size() const { return fptr->argument_size(); }
    virtual ushort smoothness() const { return fptr->smoothness(); }
    template<class Res, class Args> void _compute(Res& r, const Args& a) const { r=fptr->evaluate(a)-a; }
    template<class Res, class Args> void _compute_approx(Res& r, const Args& a) const { _compute(r,a); }
  private:
    boost::shared_ptr<RealVectorFunction> fptr;
};
*/


SolverBase::SolverBase(double max_error, uint max_steps)
  : _max_error(max_error), _max_steps(max_steps)
{
}


Vector<Interval>
SolverBase::solve(const RealVectorFunction& f,
                  const Vector<Interval>& ix) const
{
    Set< Vector<Interval> > r;
    ::solve_all(r,*this,f,ix);
    if(r.size()!=1u) { ARIADNE_THROW(SolverException,"SolverBase::solve","non-unique solution in solve("<<f<<","<<ix<<")"); }
    return *r.begin();
}

Set< Vector<Interval> >
SolverBase::solve_all(const RealVectorFunction& f,
                      const Vector<Interval>& ix) const
{
    Set< Vector<Interval> > r;
    ::solve_all(r,*this,f,ix);
    return r;
}



Vector<Interval>
SolverBase::zero(const RealVectorFunction& f,
                 const Vector<Interval>& ix) const
{
    const double& e=this->maximum_error();
    uint n=this->maximum_number_of_steps();
    ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
    Vector<Interval> r(ix);
    Vector<Interval> nr(r.size());
    bool has_solution=false;
    while(n>0) {
        nr=this->step(f,r);
        ARIADNE_LOG(5,"  nr="<<nr<<"\n");

        if(!has_solution && subset(nr,r)) {
            has_solution=true;
        }

        if(has_solution && radius(nr) < e) {
            return nr;
        }

        if(disjoint(nr,r)) {
            ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<ix<<"; "<<nr<<" is disjoint from "<<r);
        }
        r=intersection(nr,r);
        n=n-1;
    }
    if(disjoint(f.evaluate(r),Vector<Interval>(f.result_size()))) {
        ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<ix<<"; f("<<r<<") is disjoint from zero");
    } else {
        r+=(eps()+radius(r))*Vector<Interval>(r.size(),Interval(-1,+1));
        nr=this->step(f,r);
        if(subset(nr,r)) {
            return nr;
        }
        ARIADNE_THROW(SolverException,"SolverBase::zero","No result verified in "<<ix<<"; maximum number of steps reached with approximation "<<r<<" which cannot be robustly checked");
    }
}



Vector<Interval>
SolverBase::fixed_point(const RealVectorFunction& f, const Vector<Interval>& pt) const
{
    ARIADNE_NOT_IMPLEMENTED;
    //return Vector<Interval>(this->zero(DifferenceFunction(f),pt));
}


VectorTaylorFunction
SolverBase::implicit(const RealVectorFunction& f,
                      const Vector<Interval>& ip,
                      const Vector<Interval>& ix) const
{
    ARIADNE_ASSERT(f.result_size()==ix.size());
    ARIADNE_ASSERT(f.argument_size()==ip.size()+ix.size());

    VectorTaylorFunction result;

    const uint np=ip.size();
    const uint nx=ix.size();
    const double err=this->maximum_error();

    VectorTaylorFunction p(VectorTaylorFunction::identity(ip));
    VectorTaylorFunction x(VectorTaylorFunction::constant(ip,ix));
    VectorTaylorFunction nwx(x.size());
    VectorTaylorFunction fnwx(f.result_size());

    uint steps_remaining=this->maximum_number_of_steps();
    uint number_unrefined=nx;
    array<bool> refinement(nx,false);

    while(steps_remaining>0) {
        nwx=this->implicit_step(f,p,x);
        fnwx=compose(f,join(p,nwx));
        ARIADNE_LOG(5,"\n  step="<<this->maximum_number_of_steps()-steps_remaining<<"\n");
        ARIADNE_LOG(5,"  nwx="<<nwx<<"\n");
        ARIADNE_LOG(5,"  fnwx="<<fnwx<<"\n");
        ARIADNE_LOG(6,"\n");

        for(uint i=0; i!=nx; ++i) {
            if(!refinement[i]) {
                if(refines(nwx[i],x[i])) {
                    refinement[i]=true;
                    --number_unrefined;
                } else if(disjoint(nwx[i],x[i])) {
                    ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","No result found in "<<ix<<"; "<<nwx<<" is disjoint from "<<x);
                }
            }
        }

        if( (number_unrefined==0) && (radius(nwx)<err) && (radius(fnwx)<err) ) {
            return VectorTaylorFunction(nwx);
        }

        x=intersection(nwx,x);
        steps_remaining=steps_remaining-1;
    }

    ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","Could not prove existence of a solution in "<<ix<<".");


}

ScalarTaylorFunction
SolverBase::implicit(const RealScalarFunction& f,
                      const Vector<Interval>& ip,
                      const Interval& ix) const
{
    VectorTaylorFunction res=this->implicit(RealVectorFunction(1u,f),ip,Vector<Interval>(1u,ix));
    return res[0];
}


ScalarTaylorFunction implicit(const ScalarTaylorFunction& f);
IntervalTaylorModel implicit(const IntervalTaylorModel& f);


ScalarTaylorFunction
SolverBase::implicit(const IntervalScalarFunction& f,
                      const Vector<Interval>& ip,
                      const Interval& ix) const
{
    ScalarTaylorFunction tf(join(ip,ix),f);
    return Ariadne::implicit(tf);
}

ScalarTaylorFunction
SolverBase::implicit(const ScalarTaylorFunction& tf,
                      const Vector<Interval>& ip,
                      const Interval& ix) const
{
    ARIADNE_ASSERT(subset(join(ip,ix),tf.domain()));
    return Ariadne::implicit(tf);
}



ScalarTaylorFunction implicit(const ScalarTaylorFunction& f) {
    Vector<Interval> h_domain=project(f.domain(),range(0u,f.argument_size()-1u));
    Interval h_codomain=f.domain()[f.argument_size()-1u];
    IntervalTaylorModel h_model=implicit(f.model());
    ARIADNE_ASSERT(h_model.argument_size()+1==f.model().argument_size());
    IntervalTaylorModel hrs_model=h_model.rescale(Interval(-1,+1),h_codomain);
    ARIADNE_ASSERT(hrs_model.argument_size()+1==f.model().argument_size());
    return ScalarTaylorFunction(h_domain,hrs_model);
}

IntervalTaylorModel
implicit(const IntervalTaylorModel& f) {

    // Check that the arguments are suitable
    ARIADNE_ASSERT(f.argument_size()>1);

    // Set some useful size constants
    const uint fas=f.argument_size();
    const uint has=fas-1u;

    // Check to see if a solution exists
    IntervalTaylorModel g=derivative(f,has);
    Interval g_range=g.range();
    if(g_range.lower()<=0 && g_range.upper()>=0) {
        ARIADNE_THROW(ImplicitFunctionException,
                      "implicit(IntervalTaylorModel)",
                      "derivative "<<g_range<<" is not invertible");
    }

    uint number_of_steps=10;
    Vector<IntervalTaylorModel> id=IntervalTaylorModel::variables(has);
    //IntervalTaylorModel h=IntervalTaylorModel::constant(has,Interval(-1,+1));
    IntervalTaylorModel h=IntervalTaylorModel::constant(has,Interval(0));
    Vector<IntervalTaylorModel> idh=join(id,h);
    // Perform proper Newton step improvements
    //std::cerr<<"\nf="<<f<<"\n";
    //std::cerr<<"g="<<g<<"\n";
    //std::cerr<<"id="<<id<<"\n";
    //std::cerr<<"h="<<h<<"\n";
    for(uint i=0; i!=number_of_steps; ++i) {
        IntervalTaylorModel& h=idh[has];
        IntervalTaylorModel gidh=compose(g,idh);
        h.clobber();
        IntervalTaylorModel fidh=compose(f,idh);
        IntervalTaylorModel dh=fidh/gidh;
        IntervalTaylorModel nh=h-dh;
        h=nh;
        //std::cerr<<"fh["<<i<<"]="<<fidh<<"\n";
        //std::cerr<<"gh["<<i<<"]="<<gidh<<"\n";
        //std::cerr<<"dh["<<i<<"]="<<dh<<"\n";
        //std::cerr<<"nh["<<i<<"]="<<idh[has]<<"\n";
    }
    //std::cerr<<"err="<<compose(f,idh)<<"\n\n";
    //std::cerr<<"res="<<idh[has]<<"\n\n";
    return idh[has];
    //std::cerr<<"IMPLICIT(IntervalTaylorModel f="<<f<<")\n";
    //return implicit(Vector<IntervalTaylorModel>(1u,f))[0];
}


Vector<Interval>
IntervalNewtonSolver::step(const RealVectorFunction& f,
                           const Vector<Interval>& x) const
{
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    Vector<Float> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<Interval> im(m);
    Vector<Interval> w=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<w<<"\n");
    Matrix<Interval> A=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<A<<"\n");
    Matrix<Interval> Ainv=inverse(A);
    ARIADNE_LOG(5,"  inverse(Df(r))="<<Ainv<<"\n");
    Vector<Interval> dx=prod(Ainv, w);
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<Interval> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    return nx;
}

Vector<Interval>
KrawczykSolver::step(const RealVectorFunction& f,
                     const Vector<Interval>& x) const
{
    Matrix<Interval> I=Matrix<Interval>::identity(x.size());
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    Vector<Float> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<Interval> im(m);
    Vector<Interval> fm=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<fm<<"\n");
    Matrix<Interval> J=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<J<<"\n");
    Matrix<Interval> M=inverse(midpoint(J));
    ARIADNE_LOG(5,"  inverse(Df(m))="<<M<<"\n");
    Vector<Interval> dx=prod(M,fm)-prod(Matrix<Interval>(I-prod(M,J)),Vector<Interval>(x-m));
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<Interval> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<Interval> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");
    return nr;
}


Vector<Interval>
FactoredKrawczykSolver::step(const RealVectorFunction& f,
                             const Vector<Interval>& x) const
{
    Matrix<Interval> I=Matrix<Interval>::identity(x.size());
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    Vector<Float> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<Interval> im(m);
    Vector<Interval> fm=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<fm<<"\n");
    Matrix<Interval> J=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<J<<"\n");
    Matrix<Float> mJ=midpoint(J);
    Matrix<Interval> M=inverse(mJ);
    ARIADNE_LOG(5,"  inverse(Df(m))="<<M<<"\n");
    Vector<Interval> dx=prod(M,fm+prod(Matrix<Interval>(J-mJ),Vector<Interval>(x-m)));
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<Interval> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<Interval> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");
    return nr;
}


VectorTaylorFunction
IntervalNewtonSolver::implicit_step(const RealVectorFunction& f,
                            const VectorTaylorFunction& p,
                            const VectorTaylorFunction& x) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

VectorTaylorFunction
KrawczykSolver::implicit_step(const RealVectorFunction& f,
                              const VectorTaylorFunction& p,
                              const VectorTaylorFunction& x) const
{
    const uint np=p.size();
    const uint nx=x.size();
    Matrix<Interval> I=Matrix<Interval>::identity(nx);
    ARIADNE_LOG(4,"  Contracting "<<x<<"\n");
    //ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    VectorTaylorFunction mx(x);
    for(uint i=0; i!=mx.size(); ++i) { mx[i].set_error(0.0); }
    ARIADNE_LOG(5,"    mx="<<mx<<"\n");
    Vector<Float> ex(nx);
    for(uint i=0; i!=nx; ++i) { ex[i]=+x[i].error(); }
    ARIADNE_LOG(5,"    ex="<<ex<<"\n");
    VectorTaylorFunction fm=evaluate(f,join(p,mx));
    ARIADNE_LOG(5,"    f(p,mx)="<<fm<<"\n");
    Vector<Interval> rp(np);
    for(uint i=0; i!=np; ++i) { rp[i]=p[i].range(); }
    Vector<Interval> rx(nx);
    for(uint i=0; i!=nx; ++i) { rx[i]=x[i].range(); }
    Matrix<Interval> J=project(f.jacobian(join(rp,rx)),range(0,nx),range(np,np+nx));
    ARIADNE_LOG(5,"    D2f(r)="<<J<<"\n");
    Matrix<Interval> M=inverse(midpoint(J));
    ARIADNE_LOG(5,"    inverse(D2f(m))="<<M<<"\n");
    VectorTaylorFunction dx=product(M,fm) - Vector<Interval>( prod(Matrix<Interval>(I-prod(M,J)),ex*Interval(-1,+1)) );
    ARIADNE_LOG(5,"    dx="<<dx<<"\n");
    VectorTaylorFunction nwx= mx - dx;
    ARIADNE_LOG(5,"    nwx="<<nwx<<"\n");
    return nwx;
}


ScalarTaylorFunction
IntervalNewtonSolver::implicit(const RealScalarFunction& f,
                               const Vector<Interval>& ip,
                               const Interval& ix) const
{
    ARIADNE_ASSERT_MSG(f.argument_size()==ip.size()+1u,"f="<<f<<", ip="<<ip<<", ix="<<ix<<"\n");
    const uint n=ip.size();
    RealScalarFunction df=f.derivative(n);

    Vector<Interval> mp=midpoint(ip);
    Interval dfmpix=df(join(mp,ix));

    if(dfmpix.lower()<=0.0 && dfmpix.upper()>=0.0) {
        ARIADNE_THROW(ImplicitFunctionException,"IntervalNewtonSolver","df(x,[y])="<<dfmpix<<" which contains zero.");
    }

    Interval x=ix;
    bool validated=false;
    for(uint i=0; i!=this->maximum_number_of_steps(); ++i) {
        Interval mx=Interval(midpoint(x));
        Interval nx=mx-f(join(mp,mx))/df(join(mp,x));
        if(disjoint(nx,x)) {
            ARIADNE_THROW(ImplicitFunctionException,"IntervalNewtonSolver","No solution of f(x,h(x))=0 with f="<<f<<", x="<<ip<<", h in "<<ix);
        }
        if(subset(nx,x)) {
            validated=true;
            x=nx;
            if(x.radius()<this->maximum_error()) { break; }
        } else {
            x=intersection(x,nx);
        }
    }

    ScalarTaylorFunction h=ScalarTaylorFunction::constant(ip,midpoint(x));
    VectorTaylorFunction id=VectorTaylorFunction::identity(ip);
    ScalarTaylorFunction dh;

    for(uint i=0; i!=this->maximum_number_of_steps(); ++i) {
        dh=compose(f,join(id,h))/df(join(mp,Interval(h.value())));
        h=h-dh;
        h.clobber();
        Float herr=norm(dh);
        ARIADNE_LOG(6,"  dh="<<dh.range()<<"\n h="<<h<<"\n");
        if(herr<this->maximum_error()/2) { break; }
    }
    ARIADNE_LOG(4," happrox="<<h<<"\n");

    h.error()+=mag(dh.range())*4;
    Interval dfipix=df(join(ip,h.range()));
    if(dfipix.lower()<=0.0 && dfipix.upper()>=0.0) {
        ARIADNE_THROW(ImplicitFunctionException,"IntervalNewtonSolver",
                      "Cannot verify approximate solution solution of f(x,h(x))=0 with f="<<f<<", x="<<ip<<", h in "<<ix<<" since range of derivative contains zero.");
    }

    // Try to verify implicit function
    Interval dfmpmx=df(join(mp,Interval(h.value())));
    Interval herr=Interval(-h.error(),+h.error());
    Float sf=mag(1-dfipix/dfmpmx);
    Interval dherr=(1-dfipix/dfmpmx)*herr;
    ScalarTaylorFunction hnew=h;
    hnew.clobber();
    dh=compose(f,join(id,hnew));
    hnew=hnew-dh+dherr;
    if(refines(hnew,h)) { return hnew; }

    // Try one more step using larger error
    hnew.error()+=mag(dh.range())*4;
    h=hnew;
    dfmpmx=df(join(mp,Interval(h.value())));
    herr=Interval(-h.error(),+h.error());
    dfipix=df(join(ip,h.range()));
    sf=mag(1-dfipix/dfmpmx);
    dherr=(1-dfipix/dfmpmx)*herr;
    hnew=h;
    hnew.clobber();
    dh=compose(f,join(id,hnew));
    hnew=hnew-dh+dherr;

    if(!refines(hnew,h)) {
        ARIADNE_THROW(ImplicitFunctionException,"IntervalNewtonSolver",
                      "Cannot verify approximate solution solution of f(x,h(x))=0 with f="<<f<<", x="<<ip<<", h in "<<ix<<" since range of derivative contains zero.");
    }
    return hnew;
}

} // namespace Ariadne
