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

#include "function/functional.h"
#include "config.h"

#include "geometry/interval.h"
#include "function/function_model.h"

#include "solvers/solver.h"

#include "utility/logging.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/differential.h"
#include "function/taylor_model.h"
#include "function/function.h"
#include "function/function_model.h"

namespace Ariadne {

namespace {

Vector<UpperInterval> ranges(const Vector<ValidatedTaylorModel>& f) {
    Vector<UpperInterval> r(f.size()); for(Nat i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

Vector<ValidatedTaylorModel>& clobber(Vector<ValidatedTaylorModel>& h) {
    for(Nat i=0; i!=h.size(); ++i) { h[i].set_error(0u); } return h; }

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumber>
jacobian2(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumber>& x)
{
    Vector< Differential<ValidatedNumber> > dx(x.size());
    for(Nat i=0; i!=x.size()-f.size(); ++i) {
        dx[i]=Differential<ValidatedNumber>::constant(f.size(),1u,x[i]); }
    for(Nat i=0; i!=f.size(); ++i) {
        Nat j=i+(x.size()-f.size());
        dx[j]=Differential<ValidatedNumber>::variable(f.size(),1u,x[j],i); }
    Vector< Differential<ValidatedNumber> > df(f.size());
    for(Nat i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumber> J=jacobian(df);
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<ExactFloat>
jacobian2_value(const Vector<ValidatedTaylorModel>& f)
{
    const Nat rs=f.size();
    const Nat fas=f.zero_element().argument_size();
    const Nat has=fas-rs;
    Matrix<ExactFloat> J(rs,rs);
    MultiIndex a(fas);
    for(Nat i=0; i!=rs; ++i) {
        for(Nat j=0; j!=rs; ++j) {
            a[has+j]=1; const ExactFloat x=f[i][a]; J[i][j]=x; a[has+j]=0;
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<ValidatedNumber>
jacobian2_range(const Vector<ValidatedTaylorModel>& f)
{
    Nat rs=f.size();
    Nat fas=f.zero_element().argument_size();
    Nat has=fas-rs;
    Matrix<ValidatedNumber> J(rs,rs);
    for(Nat i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(Nat k=0; k!=rs; ++k) {
                const Nat c=iter->key()[has+k];
                if(c>0) {
                    const ExactFloat& x=iter->data();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=ValidatedNumber(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}



} // namespace


static const Bool ALLOW_PARTIAL_FUNCTION = true;

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory();

ValidatedVectorFunctionModel operator*(const Matrix<ExactFloat>& A,const ValidatedVectorFunctionModel& v) {
    ARIADNE_ASSERT(v.size()!=0);
    ValidatedVectorFunctionModel r(A.row_size(),v[0].create_zero());
    for(Nat i=0; i!=r.size(); ++i) {
        ValidatedScalarFunctionModel t=r[i];
        for(Nat j=0; j!=v.size(); ++j) {
            t+=ExactFloat(A[i][j])*v[j];
        }
        r[i]=t;
    }
    return r;
}

ValidatedVectorFunctionModel operator*(const Matrix<ValidatedNumber>& A,const ValidatedVectorFunctionModel& v) {
    ARIADNE_ASSERT(v.size()!=0);
    ValidatedVectorFunctionModel r(A.row_size(),v[0].create_zero());
    for(Nat i=0; i!=r.size(); ++i) {
        ValidatedScalarFunctionModel t=r[i];
        for(Nat j=0; j!=v.size(); ++j) {
            t+=A[i][j]*v[j];
        }
        r[i]=t;
    }
    return r;
}

ErrorFloat sup_error(const ValidatedVectorFunctionModel& x) {
    ErrorFloat r=0u;
    for(Nat i=0; i!=x.size(); ++i) { r=std::max(r,x[i].error()); }
    return r;
}







SolverBase::SolverBase(double max_error, Nat max_steps)
  : _max_error(max_error), _max_steps(max_steps), _function_factory_ptr(make_taylor_function_factory())
{
}

Void
SolverBase::set_function_factory(const FunctionModelFactoryInterface<ValidatedTag>& factory)
{
    this->_function_factory_ptr=std::shared_ptr< FunctionModelFactoryInterface<ValidatedTag> >(factory.clone());
}

const FunctionModelFactoryInterface<ValidatedTag>&
SolverBase::function_factory() const
{
    return *this->_function_factory_ptr;
}


Vector<ValidatedNumber>
SolverBase::solve(const ValidatedVectorFunction& f,
                  const Vector<ValidatedNumber>& ix) const
{
    ExactBox bx=make_exact_box(ix);
    return this->solve(f,bx);
}

Vector<ValidatedNumber>
SolverBase::solve(const ValidatedVectorFunction& f,
                  const ExactBox& bx) const
{
    Set< Vector<ValidatedNumber> > r = this->solve_all(f,bx);
    if(r.size()==0u) { ARIADNE_THROW(NoSolutionException,"SolverBase::solve","no solution in solve("<<f<<","<<bx<<")"); }
    if(r.size()!=1u) { ARIADNE_THROW(SolverException,"SolverBase::solve","non-unique solution in solve("<<f<<","<<bx<<")"); }
    return *r.begin();
}

Void
solve_all(Set< Vector<ValidatedNumber> >& r,
          const SolverInterface& s,
          const ValidatedVectorFunction& f,
          const ExactBox& ix);


template<class X1, class X2>
Bool operator<(const Vector<X1>& v1, const Vector<X2>& v2)
{
    if(v1.size()!=v2.size()) { return v1.size()<v2.size(); }
    for(SizeType i=0; i!=v1.size(); ++i) {
        if(decide(v1[i]<v2[i])) { return true; }
        else if(decide(v1[i]>v2[i])) { return false; }
    }
    return true;
}


Set< Vector<ValidatedNumber> >
SolverBase::solve_all(const ValidatedVectorFunction& f,
                      const ExactBox& bx) const
{
    ARIADNE_LOG(5,"SolverBase::solve_all(f,bx): f="<<f<<", ix="<<bx<<"\n");

    // Create result set
    Set< Vector<ValidatedNumber> > r;

    Vector<ValidatedFloat> x=make_singleton(bx);

    // Test for no solution
    const Vector<ValidatedNumber> z(bx.size());
    if(inconsistent(f.evaluate(x),z)) {
        return r;
    }

    Bool invertible_jacobian=true;
    //Vector<ValidatedNumber> nx=2*ix-make_exact(bx);
    Vector<ValidatedNumber> nx=x;
    try {
        Matrix<ValidatedNumber> Jinv=inverse(f.jacobian(nx));
    }
    catch(const SingularMatrixException& e) {
        invertible_jacobian=false;
    }

    Bool need_to_split=true;

    if(invertible_jacobian) {
        //std::cerr<<"Nonsingular matrix -- applying contractor\n";
        try {
            Vector<ValidatedNumber> y=this->zero(f,bx);
            Bool is_new=true;
            for(Set<Vector<ValidatedNumber> >::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
                if(consistent(y,*iter)) {
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
        // If sup_error is too small, assume solution is not verified
        if(sup_error(make_singleton(bx))<this->maximum_error()) {
            if(!invertible_jacobian) {
                ARIADNE_WARN("Cannot verify solution in "<<bx<<" with f="<<f(make_singleton(bx))<<"); "
                             <<"Jacobian "<<f.jacobian(nx)<<" is not invertible; "
                             <<"approximate inverse="<<inverse(midpoint(f.jacobian(nx)))<<"\n");
            } else {
                ARIADNE_WARN("Cannot verify or falsify solution in "<<bx<<"; f("<<bx<<")="<<f(make_singleton(bx))<<".\n");
            }
            return r;
        }

        //std::cerr<<"  Splitting "<<bx<<"\n";
        Pair< Vector<ExactInterval>, Vector<ExactInterval> > splt=split(bx);
        r.adjoin(this->solve_all(f,splt.first));
        r.adjoin(this->solve_all(f,splt.second));
    }

    return r;
}





Vector<ValidatedNumber>
SolverBase::zero(const ValidatedVectorFunction& f,
                 const ExactBox& bx) const
{
    const double& e=this->maximum_error();
    Nat n=this->maximum_number_of_steps();
    ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
    Vector<ValidatedNumber> r=make_singleton(bx);
    Vector<ValidatedNumber> nr(r.size());
    Bool has_solution=false;
    while(n>0) {
        nr=this->step(f,r);
        ARIADNE_LOG(5,"  nr="<<nr<<"\n");

        if(!has_solution && refines(nr,r)) {
            has_solution=true;
        }

        if(has_solution && sup_error(nr) < e) {
            return nr;
        }

        if(!consistent(nr,r)) {
            ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<bx<<"; "<<nr<<" is inconsistent with "<<r);
        }
        r=refinement(nr,r);
        n=n-1;
    }
    if(!consistent(f.evaluate(r),Vector<ValidatedNumber>(f.result_size()))) {
        ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<bx<<"; f("<<r<<") is inconsistent with zero");
    } else {
        UpperFloat widen=ExactFloat(Float::eps())*sup_error(r);
        r+=Vector<ValidatedNumber>(r.size(),ValidatedNumber(-widen,+widen));
        nr=this->step(f,r);
        if(refines(nr,r)) {
            return nr;
        }
        ARIADNE_THROW(SolverException,"SolverBase::zero","No result verified in "<<bx<<"; maximum number of steps reached with approximation "<<r<<" which cannot be robustly checked");
    }
}



Vector<ValidatedNumber>
SolverBase::fixed_point(const ValidatedVectorFunction& f, const ExactBox& bx) const
{
    ValidatedVectorFunction id=ValidatedVectorFunction::identity(f.argument_size());
    return this->solve(f-id,bx);
}


ValidatedVectorFunctionModel
SolverBase::implicit(const ValidatedVectorFunction& f,
                      const ExactBox& ip,
                      const ExactBox& ix) const
{
    ARIADNE_LOG(4,"SolverBase::implicit(ValidatedVectorFunction f, ExactIntervalVector ip, ExactIntervalVector ix)\n");
    ARIADNE_LOG(5,"f="<<f<<"\n");
    ARIADNE_ASSERT(f.result_size()==ix.size());
    ARIADNE_ASSERT(f.argument_size()==ip.size()+ix.size());

    const Nat n=ix.size();
    const double err=this->maximum_error();

    ValidatedVectorFunctionModel id(this->function_factory().create_identity(ip));
    ValidatedVectorFunctionModel h(this->function_factory().create_constants(ip,make_singleton(ix)));
    ValidatedVectorFunctionModel nh(this->function_factory().create_zeros(n,ip));
    ValidatedVectorFunctionModel fnh(this->function_factory().create_zeros(n,ip));

    Nat steps_remaining=this->maximum_number_of_steps();
    Nat number_unrefined=n;
    Array<Bool> refinement(n,false);

    while(steps_remaining>0) {
        ARIADNE_LOG(5,"\n");
        ARIADNE_LOG(5,"step="<<this->maximum_number_of_steps()-steps_remaining<<"\n");
        nh=this->implicit_step(f,id,h);
        fnh=compose(f,join(id,nh));
        ARIADNE_LOG(5,"  nh="<<nh<<"\n");
        ARIADNE_LOG(5,"  fnh="<<fnh<<"\n");

        if(ALLOW_PARTIAL_FUNCTION) {
            for(Nat i=0; i!=n; ++i) {
                if(!refinement[i]) {
                    ARIADNE_LOG(6,"refines(nh["<<i<<"],x["<<i<<"])="<<refines(nh[i],h[i])<<"\n");
                    if(refines(nh[i],h[i])) {
                        refinement[i]=true;
                        --number_unrefined;
                    } else if(definitely(disjoint(ValidatedScalarFunctionModel(nh[i]).range(),ix[i]))) {
                        ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","No result found in "<<ix<<"; function "<<nh<<" is disjoint from "<<h<<" for at least one point.");
                    }
                }
            }

            ARIADNE_LOG(6,"refinement="<<refinement<<"\n");
            ARIADNE_LOG(6,"nh.range()="<<nh.range()<<"\n");
            ARIADNE_LOG(6,"sup_error(nh)="<<sup_error(nh)<<"\n");
            ARIADNE_LOG(6,"sup_error(fnh)="<<sup_error(fnh)<<"\n");

            h=nh;

        } else { // !ALLOW_PARTIAL_FUNCTION
            for(Nat i=0; i!=n; ++i) {
                if(!refinement[i]) {
                    if(refines(nh[i],h[i])) {
                        refinement[i]=true;
                        --number_unrefined;
                    } else if(disjoint(nh[i],h[i])) {
                        ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","No result found in "<<ix<<"; "<<nh<<" is disjoint from "<<h);
                    }
                }
            }

            h=intersection(nh,h);
        }

        steps_remaining=steps_remaining-1;
        if( (number_unrefined==0) && ( (steps_remaining==0) || ((sup_error(nh)<err) && (sup_error(fnh)<err)) ) ) {
            return h;
        }
    }

    ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","Could not prove existence of a solution in "<<ix<<".");


}

ValidatedScalarFunctionModel
SolverBase::implicit(const ValidatedScalarFunction& f,
                     const ExactBox& ip,
                     const ExactInterval& ix) const
{
    ARIADNE_LOG(4,"SolverBase::implicit(ValidatedScalarFunction f, ExactIntervalVector ip, ExactInterval ix)\n");
    ARIADNE_LOG(5,"f="<<f<<"\n");
    ValidatedVectorFunctionModel res=this->implicit(ValidatedVectorFunction(List<ValidatedScalarFunction>(1u,f)),ip,ExactBox(1u,ix));
    return res[0];
}

ValidatedVectorFunctionModel
SolverBase::continuation(const ValidatedVectorFunction& f,
                         const Vector<ApproximateNumber>& p,
                         const ExactBox& ix,
                         const ExactBox& ip) const
{
    ARIADNE_NOT_IMPLEMENTED;
}




Vector<ValidatedNumber>
IntervalNewtonSolver::step(const ValidatedVectorFunction& f,
                           const Vector<ValidatedNumber>& x) const
{
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    Vector<ExactFloat> m(make_exact(x));
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<ValidatedNumber> im(m);
    Vector<ValidatedNumber> w=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<w<<"\n");
    Matrix<ValidatedNumber> A=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<A<<"\n");
    Matrix<ValidatedNumber> Ainv=inverse(A);
    ARIADNE_LOG(5,"  inverse(Df(r))="<<Ainv<<"\n");
    Vector<ValidatedNumber> dx=Ainv*w;
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<ValidatedNumber> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    return nx;
}

Vector<ValidatedNumber>
KrawczykSolver::step(const ValidatedVectorFunction& f,
                     const Vector<ValidatedNumber>& x) const
{
    Matrix<ValidatedNumber> I=Matrix<ValidatedNumber>::identity(x.size());
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    Vector<ExactFloat> m(make_exact(x));
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<ValidatedNumber> im(m);
    Vector<ValidatedNumber> fm=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<fm<<"\n");
    Matrix<ValidatedNumber> J=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<J<<"\n");
    Matrix<ValidatedNumber> M=inverse(midpoint(J));
    ARIADNE_LOG(5,"  inverse(Df(m))="<<M<<"\n");
    Vector<ValidatedNumber> dx=M*fm-(I-M*J)*(x-m);
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<ValidatedNumber> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<ValidatedNumber> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");
    return nr;
}


Vector<ValidatedNumber>
FactoredKrawczykSolver::step(const ValidatedVectorFunction& f,
                             const Vector<ValidatedNumber>& x) const
{
    Matrix<ValidatedNumber> I=Matrix<ValidatedNumber>::identity(x.size());
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    Vector<ExactFloat> m(make_exact(x));
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<ValidatedNumber> im(m);
    Vector<ValidatedNumber> fm=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<fm<<"\n");
    Matrix<ValidatedNumber> J=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<J<<"\n");
    Matrix<ValidatedNumber> mJ(midpoint(J));
    Matrix<ValidatedNumber> M=inverse(mJ);
    ARIADNE_LOG(5,"  inverse(Df(m))="<<M<<"\n");
    Vector<ValidatedNumber> dx=M*(fm+(J-mJ)*(x-m));
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<ValidatedNumber> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<ValidatedNumber> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");
    return nr;
}


ValidatedVectorFunctionModel
IntervalNewtonSolver::implicit_step(const ValidatedVectorFunction& f,
                                    const ValidatedVectorFunctionModel& id,
                                    const ValidatedVectorFunctionModel& h) const
{
    const Nat m=id.size();
    const Nat n=h.size();

    ARIADNE_LOG(6,"IntervalNewtonSolver::implicit_step(ValidatedVectorFunction f, ValidatedVectorFunctionModel id, ValidatedVectorFunctionModel h)\n");
    ARIADNE_LOG(7,"f="<<f<<"\n");
    ARIADNE_LOG(7,"h="<<h<<"\n");

    ValidatedVectorFunctionModel mh=h; mh.clobber();
    ARIADNE_LOG(7,"midpoint(h)="<<mh<<"\n");

    ValidatedScalarFunction zero_function(f.argument_size());
    Matrix<ValidatedScalarFunction> D2f(n,n,zero_function);
    for(Nat i=0; i!=n; ++i) {
        for(Nat j=0; j!=n; ++j) {
            D2f[i][j]=f[i].derivative(m+j);
        }
    }
    ARIADNE_LOG(7,"D2f="<<D2f<<"\n");

    ValidatedNumber zero(0);
    ValidatedScalarFunctionModel z=h[0]*zero;
    ValidatedVectorFunctionModel idh=join(id,h);

    Matrix<ValidatedScalarFunctionModel> J(n,n,z);
    for(Nat i=0; i!=n; ++i) {
        for(Nat j=0; j!=n; ++j) {
            J[i][j]=compose(D2f[i][j],idh);
        }
    }
    ARIADNE_LOG(7,"J="<<J<<"\n");

    Matrix<UpperInterval> rngJ(n,n);
    for(Nat i=0; i!=n; ++i) {
        for(Nat j=0; j!=n; ++j) {
            UpperInterval D2fij=UpperInterval(unchecked_evaluate(D2f[i][j],make_singleton(product(id.range(),h.range()))));
            rngJ[i][j]=intersection(J[i][j].range(),D2fij);
        }
    }
    ARIADNE_LOG(7,"rngJ="<<rngJ<<"\n");

    ValidatedVectorFunctionModel fidmh=compose(f,join(id,mh));
    ARIADNE_LOG(7,"compose(f,join(id,midpoint(h)))="<<fidmh<<"\n");

    ValidatedVectorFunctionModel dh(n,z);
    if(n==1) {
        if(contains(rngJ[0][0],ExactFloat(0.0))) {
            ARIADNE_THROW(SingularJacobianException,"IntervalNewtonSolver","D2f(P,X)="<<rngJ[0][0]<<" which contains zero.");
        }
        if(contains(J[0][0].range(),ExactFloat(0.0))) {
            dh[0]=fidmh[0]/make_singleton(rngJ[0][0]);
        } else {
            dh[0]=fidmh[0]/J[0][0];
        }
    } else {
        dh=inverse(make_singleton(rngJ))*fidmh;
    }
    ARIADNE_LOG(7,"dh="<<dh<<"\n");

    ValidatedVectorFunctionModel nh=mh-dh;
    ARIADNE_LOG(7,"nh="<<nh<<"\n");
    return nh;
}


ValidatedVectorFunctionModel
KrawczykSolver::implicit_step(const ValidatedVectorFunction& f,
                              const ValidatedVectorFunctionModel& p,
                              const ValidatedVectorFunctionModel& x) const
{
    const Nat np=p.size();
    const Nat nx=x.size();
    Matrix<ValidatedNumber> I=Matrix<ValidatedNumber>::identity(nx);
    ARIADNE_LOG(4,"  Contracting x="<<x<<"\n");
    ARIADNE_LOG(4,"    p="<<p<<"\n");
    ARIADNE_LOG(4,"    f="<<f<<"\n");
    //ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    ValidatedVectorFunctionModel mx(x);
    for(Nat i=0; i!=mx.size(); ++i) { mx[i].set_error(0u); }
    ARIADNE_LOG(5,"    mx="<<mx<<"\n");
    Vector<ErrorFloat> ex(nx);
    for(Nat i=0; i!=nx; ++i) { ex[i]=x[i].error(); }
    Vector<ValidatedNumber> eix=make_bounds(ex);
    ARIADNE_LOG(5,"    ex="<<ex<<"\n");
    ValidatedVectorFunctionModel fm=compose(f,join(p,mx));
    ARIADNE_LOG(5,"    f(p,mx)="<<fm<<"\n");
    Vector<ValidatedNumber> rp(np);
    for(Nat i=0; i!=np; ++i) { rp[i]=make_singleton(p[i].range()); }
    Vector<ValidatedNumber> rx(nx);
    for(Nat i=0; i!=nx; ++i) { rx[i]=make_singleton(x[i].range()); }
    Matrix<ValidatedNumber> J=project(f.jacobian(join(rp,rx)),range(0,nx),range(np,np+nx));
    ARIADNE_LOG(5,"    D2f(r)=J="<<J<<"\n");
    Matrix<ValidatedNumber> M=inverse(midpoint(J));
    ARIADNE_LOG(5,"    inverse(D2f(m))=M="<<M<<"\n");
    ARIADNE_LOG(5,"    M*f(p,mx)="<<M*fm<<"\n");
    ARIADNE_LOG(5,"    (I-M*J)="<<(I-M*J)<<"\n");
    ARIADNE_LOG(5,"    (I-M*J) * (ex*ValidatedNumber(-1,+1))="<<(I-M*J)<<"*"<<eix<<"="<<(I-M*J) * eix<<"\n");
    ValidatedVectorFunctionModel dx= M*fm - (I-M*J) * eix;
    ARIADNE_LOG(5,"    dx="<<dx<<"\n");
    ValidatedVectorFunctionModel nwx= mx - dx;
    ARIADNE_LOG(5,"    nwx="<<nwx<<"\n");
    return nwx;
}


/*
ValidatedScalarFunctionModel
IntervalNewtonSolver::implicit(const ValidatedScalarFunction& f,
                               const ExactBox& ip,
                               const ExactInterval& ix) const
{
    return this->SolverBase::implicit(f,ip,ix);
    ARIADNE_LOG(4,"IntervalNewtonSolver::implicit(ValidatedScalarFunction f, ExactIntervalVector P, ExactInterval X)\n");
    ARIADNE_LOG(5,"f="<<f<<"\n");
    ARIADNE_LOG(5,"P="<<ip<<"\n");
    ARIADNE_LOG(5,"X="<<std::setprecision(17)<<ix<<"\n");
    ARIADNE_ASSERT_MSG(f.argument_size()==ip.size()+1u,"f="<<f<<", P="<<ip<<", X="<<ix<<"\n");
    const Nat n=ip.size();
    ValidatedScalarFunction df=f.derivative(n);
    ARIADNE_LOG(5,"df="<<df<<"\n");

    // Simple check to see if there is definitely no solution
    ExactInterval fipix=evaluate(f,join(ip,ix));
    ARIADNE_LOG(5,"f(P,X)="<<fipix<<"\n");
    ARIADNE_LOG(5,"f(p,X)="<<f(join(ExactIntervalVector(midpoint(ip)),ix))<<"\n");
    if(fipix.lower()>0.0 || fipix.upper()<0.0) {
        ARIADNE_THROW(NoSolutionException,"IntervalNewtonSolver","f(P,X)="<<fipix<<" which does not contain zero.");
    }

    // Simple test of nonsingularity of Jacobian
    ExactInterval dfipix=evaluate(df,join(ip,ix));
    ARIADNE_LOG(5,"df(P,X)="<<dfipix<<"\n");
    if(dfipix.lower()<=0.0 && dfipix.upper()>=0.0) {
        ARIADNE_THROW(SingularJacobianException,"IntervalNewtonSolver","df(P,X)="<<dfipix<<" which contains zero.");
    }


    // Set up auxiliary functions
    ValidatedScalarFunctionModel h=this->function_factory().create_constant(ip,ix);
    ValidatedVectorFunctionModel id=this->function_factory().create_identity(ip);
    ValidatedScalarFunctionModel dh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel mh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel nh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel dfidh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel fidmh=this->function_factory().create_zero(ip);

    ARIADNE_LOG(5,"f="<<f<<"\n");
    ARIADNE_LOG(5,"df="<<df<<"\n");
    ARIADNE_LOG(5,"P="<<ip<<"\n");
    ARIADNE_LOG(5,"X="<<ix<<"\n");
    ARIADNE_LOG(5,"h="<<h<<"\n");
    ARIADNE_LOG(5,"id="<<id<<"\n");

    // Compute solution
    Bool refinement=false;
    for(Nat i=0; i!=this->maximum_number_of_steps(); ++i) {
        mh=h; mh.clobber();
        ARIADNE_LOG(7,"mh="<<mh<<"\n");
        ValidatedNumber dfiph=evaluate(df,join(ip,intersection(ix,h.range())));
        ARIADNE_LOG(7,"df(P,h(P))="<<dfiph<<"\n");
        dfidh=compose(df,join(id,h));
        ARIADNE_LOG(7,"df(id,h)="<<dfidh<<"\n");
        fidmh=compose(f,join(id,mh));
        ARIADNE_LOG(7,"f(p,mh(p))="<<fidmh<<"\n");
        // Prefer to divide by df(id,h), but may divide by constant df(P,h(P)) if latter is tighter
        if(refines(dfidh.range(),dfiph)) {
            dh=fidmh/dfidh;
        } else {
            dh=fidmh/dfiph;
        }
        ARIADNE_LOG(6,"dh="<<dh<<"\n");
        ARIADNE_LOG(6,"dh.range()="<<dh.range()<<"\n");
        ARIADNE_LOG(8,"norm(dh)="<<norm(dh)<<", h.error()="<<h.error()<<"\n");
        ErrorFloat herr=h.error();
        ErrorFloat dhnrm=norm(dh);
        // Prefer to check refinement using norm of dh, since this is faster
        // if(refines(nh,h)) { refinement=true; }
        if(dhnrm<=herr) { refinement=true; }
        nh=mh-dh;
        ARIADNE_LOG(6,"nh="<<nh<<"\n");
        ARIADNE_LOG(6,"nh.range()="<<nh.range()<<"\n");
        if(disjoint(nh.range(),ix)) {
            ARIADNE_THROW(NoSolutionException,"IntervalNewtonSolver",
                          "No solution of f(a,x)=0 with f="<<f<<", a in "<<ip<<", x in "<<ix<<".");
        }
        h=nh;
        if(refinement && dhnrm<this->maximum_error()/2) { break; }
    }
    return h;
    if(!refinement) {
        ARIADNE_LOG(6,"refinement="<<refinement<<"\n");
        ARIADNE_THROW(UnknownSolutionException,"IntervalNewtonSolver",
                      "Cannot find solution of f(x,h(x))=0 with f="<<f<<", x="<<ip<<", h in "<<ix<<".");
    }
    if(!refines(h.range(),ix)) {
        ARIADNE_LOG(6,"(h.range(),ix)="<<h.range()<<" "<<ix<<"\n");
        ARIADNE_THROW(UnknownSolutionException,"IntervalNewtonSolver",
                      "Range "<<h.range()<<" of solution h of f(x,h(x))=0 with f="<<f<<", x="<<ip<<", h in "<<ix<<" lies outside permitted interval.");
    }

}
*/


Void IntervalNewtonSolver::write(OutputStream& os) const
{
    os << "IntervalNewtonSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}

Void KrawczykSolver::write(OutputStream& os) const
{
    os << "KrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}

Void FactoredKrawczykSolver::write(OutputStream& os) const
{
    os << "FactoredKrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}


} // namespace Ariadne
