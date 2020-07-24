/***************************************************************************
 *            solvers/solver.cpp
 *
 *  Copyright  2006-20  Pieter Collins
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../geometry/interval.hpp"
#include "../function/function_model.hpp"

#include "../solvers/solver.hpp"

#include "../output/logging.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/differential.hpp"
#include "../algebra/algebra.hpp"
#include "../function/taylor_model.hpp"
#include "../function/formula.hpp"
#include "../function/function.hpp"
#include "../function/function_model.hpp"

#include "../algebra/expansion.inl.hpp"

#include "../algebra/evaluate.hpp"
#include "../algebra/evaluate.tpl.hpp"

namespace Ariadne {

namespace {

template<class PR> auto
operator*(const Matrix<FloatBounds<PR>>& A,const ValidatedVectorMultivariateFunctionModel<PR>& v)
    -> ValidatedVectorMultivariateFunctionModel<PR>
{
    ARIADNE_ASSERT(v.size()!=0);
    ValidatedVectorMultivariateFunctionModel<PR> r(A.row_size(),factory(v).create_zero());
    for(Nat i=0; i!=r.size(); ++i) {
        ValidatedScalarMultivariateFunctionModel<PR> t=r[i];
        for(Nat j=0; j!=v.size(); ++j) {
            t+=A[i][j]*v[j];
        }
        r[i]=t;
    }
    return r;
}

template<class PR> FloatError<PR> sup_error(const ValidatedVectorMultivariateFunctionModel<PR>& x) {
    assert(x.size()>0);
    FloatError<PR> r(0u,x[0].error().precision());
    for(Nat i=0; i!=x.size(); ++i) { r=max(r,x[i].error()); }
    return r;
}


}


static const Bool ALLOW_PARTIAL_FUNCTION = true;

FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory();


SolverBase::SolverBase(ApproximateDouble max_error, Nat max_steps)
  : _max_error(cast_exact(max_error)), _max_steps(max_steps), _function_factory_ptr(make_taylor_function_factory())
{
}

Void
SolverBase::set_function_factory(const FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>& factory)
{
    this->_function_factory_ptr=std::shared_ptr< FunctionModelFactoryInterface<ValidatedTag,DoublePrecision> >(factory.clone());
}

const FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>&
SolverBase::function_factory() const
{
    return *this->_function_factory_ptr;
}


auto SolverBase::solve(const ValidatedVectorMultivariateFunction& f,
                       const Vector<ValidatedNumericType>& ix) const
    -> Vector<ValidatedNumericType>
{
    ExactBoxType bx=cast_exact_box(ix);
    return this->solve(f,bx);
}

auto SolverBase::solve(const ValidatedVectorMultivariateFunction& f,
                       const ExactBoxType& bx) const
    -> Vector<ValidatedNumericType>
{
    Set< Vector<ValidatedNumericType> > r = this->solve_all(f,bx);
    if(r.size()==0u) { ARIADNE_THROW(NoSolutionException,"SolverBase::solve","no solution in solve("<<f<<","<<bx<<")"); }
    if(r.size()!=1u) { ARIADNE_THROW(SolverException,"SolverBase::solve","non-unique solution in solve("<<f<<","<<bx<<")"); }
    return *r.begin();
}

Void
solve_all(Set< Vector<SolverInterface::ValidatedNumericType> >& r,
          const SolverInterface& s,
          const ValidatedVectorMultivariateFunction& f,
          const ExactBoxType& ix);

auto SolverBase::solve_all(const ValidatedVectorMultivariateFunction& f,
                           const ExactBoxType& bx) const
    -> Set< Vector<ValidatedNumericType> >
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("f="<<f<<", ix="<<bx);

    // Create result set
    Set< Vector<ValidatedNumericType> > r;

    Vector<FloatDPBounds> x=cast_singleton(bx);

    // Test for no solution
    const Vector<ValidatedNumericType> z(bx.size(),dp);
    if(inconsistent(f.evaluate(x),z)) {
        return r;
    }

    Bool invertible_jacobian=true;
    //Vector<ValidatedNumericType> nx=2*ix-cast_exact(bx);
    Vector<ValidatedNumericType> nx=x;
    try {
        Matrix<ValidatedNumericType> Jinv=inverse(f.jacobian(nx));
    }
    catch(const SingularMatrixException& e) {
        invertible_jacobian=false;
    }

    Bool need_to_split=true;

    if(invertible_jacobian) {
        //std::cerr<<"Nonsingular matrix -- applying contractor\n";
        try {
            Vector<ValidatedNumericType> y=this->zero(f,bx);
            Bool is_new=true;
            for(Set<Vector<ValidatedNumericType> >::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
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
            //ARIADNE_WARN("NoSolutionException exception: "<<e.what());
            need_to_split=false;
        }
        catch(const SolverException& e) {
            ARIADNE_WARN("SolverException exception: "<<e.what());
            need_to_split=true;
            // No solution found, try splitting
        }
    } else {
        need_to_split=true;
    }

    if(need_to_split) {
        // If sup_error is too small, assume solution is not verified
        if(definitely(sup_error(cast_singleton(bx))<this->maximum_error())) {
            if(!invertible_jacobian) {
                ARIADNE_WARN("Cannot verify solution in "<<bx<<" with f="<<f(cast_singleton(bx))<<"); "
                             <<"Jacobian "<<f.jacobian(nx)<<" is not invertible; "
                             <<"approximate inverse="<<inverse(midpoint(f.jacobian(nx))));
            } else {
                ARIADNE_WARN("Cannot verify or falsify solution in "<<bx<<"; f("<<bx<<")="<<f(cast_singleton(bx))<<".");
            }
            return r;
        }

        //std::cerr<<"  Splitting "<<bx<<"\n";
        Pair< Box<ExactIntervalType>, Box<ExactIntervalType> > splt=split(bx);
        r.adjoin(this->solve_all(f,splt.first));
        r.adjoin(this->solve_all(f,splt.second));
    }

    return r;
}





auto SolverBase::zero(const ValidatedVectorMultivariateFunction& f,
                      const ExactBoxType& bx) const
    -> Vector<ValidatedNumericType>
{
    ARIADNE_LOG_SCOPE_CREATE;

    const ExactDouble e=this->maximum_error();
    Nat n=this->maximum_number_of_steps();
    Vector<ValidatedNumericType> r=cast_singleton(bx);
    Vector<ValidatedNumericType> nr(r.size(),dp);
    Bool has_solution=false;
    while(n>0) {
        nr=this->step(f,r);
        ARIADNE_LOG_PRINTLN_AT(1,"nr="<<nr);

        if(!has_solution && refines(nr,r)) {
            has_solution=true;
        }

        if(has_solution && definitely(sup_error(nr) < e)) {
            return nr;
        }

        if(!consistent(nr,r)) {
            ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<bx<<"; "<<nr<<" is inconsistent with "<<r);
        }
        r=refinement(nr,r);
        n=n-1;
    }
    if(!consistent(f.evaluate(r),Vector<ValidatedNumericType>(f.result_size(),dp))) {
        ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<bx<<"; f("<<r<<") is inconsistent with zero");
    } else {
        FloatDPError widen=FloatDPError(FloatDP::eps(dp))*sup_error(r);
        r+=Vector<ValidatedNumericType>(r.size(),ValidatedNumericType(-widen,+widen));
        nr=this->step(f,r);
        if(refines(nr,r)) {
            return nr;
        }
        ARIADNE_THROW(SolverException,"SolverBase::zero","No result verified in "<<bx<<"; maximum number of steps reached with approximation "<<r<<" which cannot be robustly checked");
    }
}



auto SolverBase::fixed_point(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const
    -> Vector<ValidatedNumericType>
{
    ValidatedVectorMultivariateFunction id=ValidatedVectorMultivariateFunction::identity(f.argument_size());
    return this->solve(f-id,bx);
}



auto SolverBase::implicit(const ValidatedVectorMultivariateFunction& f,
                          const ExactBoxType& ip,
                          const ExactBoxType& ix) const
    -> ValidatedVectorMultivariateFunctionModelType
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"f="<<f);
    ARIADNE_ASSERT(f.result_size()==ix.size());
    ARIADNE_ASSERT(f.argument_size()==ip.size()+ix.size());

    const Nat n=ix.size();
    const ExactDouble err=this->maximum_error();

    ValidatedVectorMultivariateFunctionModelDP id(this->function_factory().create_identity(ip));
    ValidatedVectorMultivariateFunctionModelDP h(this->function_factory().create_constants(ip,cast_singleton(ix)));
    ValidatedVectorMultivariateFunctionModelDP nh(this->function_factory().create_zeros(n,ip));
    ValidatedVectorMultivariateFunctionModelDP fnh(this->function_factory().create_zeros(n,ip));

    Nat steps_remaining=this->maximum_number_of_steps();
    Nat number_unrefined=n;
    Array<Bool> is_refinement(n,false);

    while(steps_remaining>0) {
        ARIADNE_LOG_PRINTLN_AT(1,"step="<<this->maximum_number_of_steps()-steps_remaining);
        nh=this->implicit_step(f,id,h);
        fnh=compose(f,join(id,nh));
        ARIADNE_LOG_PRINTLN_AT(2,"nh="<<nh);
        ARIADNE_LOG_PRINTLN_AT(2,"fnh="<<fnh);

        if(ALLOW_PARTIAL_FUNCTION) {
            for(Nat i=0; i!=n; ++i) {
                if(!is_refinement[i]) {
                    ARIADNE_LOG_PRINTLN_AT(3,"refines(nh["<<i<<"],x["<<i<<"])="<<refines(nh[i],h[i]));
                    if(refines(nh[i],h[i])) {
                        is_refinement[i]=true;
                        --number_unrefined;
                    } else if(definitely(disjoint(ValidatedScalarMultivariateFunctionModelDP(nh[i]).range(),ix[i]))) {
                        ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","No result found in "<<ix<<"; function "<<nh<<" is disjoint from "<<h<<" for at least one point.");
                    }
                }
            }

            ARIADNE_LOG_PRINTLN_AT(2,"is_refinement="<<is_refinement);
            ARIADNE_LOG_PRINTLN_AT(2,"nh.range()="<<nh.range());
            ARIADNE_LOG_PRINTLN_AT(2,"sup_error(nh)="<<sup_error(nh));
            ARIADNE_LOG_PRINTLN_AT(2,"sup_error(fnh)="<<sup_error(fnh));

            h=nh;

        } else { // !ALLOW_PARTIAL_FUNCTION
            for(Nat i=0; i!=n; ++i) {
                if(!is_refinement[i]) {
                    if(refines(nh[i],h[i])) {
                        is_refinement[i]=true;
                        --number_unrefined;
                    } else if(inconsistent(nh[i],h[i])) {
                        ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","No result found in "<<ix<<"; "<<nh<<" is disjoint from "<<h);
                    }
                }
            }

            h=refinement(nh,h);
        }

        steps_remaining=steps_remaining-1;
        if( (number_unrefined==0) && ( (steps_remaining==0) || (definitely(sup_error(nh)<err) && definitely(sup_error(fnh)<err)) ) ) {
            return h;
        }
    }

    ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","Could not prove existence of a solution in "<<ix<<".");
}

ValidatedScalarMultivariateFunctionModelDP
SolverBase::implicit(const ValidatedScalarMultivariateFunction& f,
                     const ExactBoxType& ip,
                     const ExactIntervalType& ix) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"f="<<f);
    ValidatedVectorMultivariateFunctionModelDP res=this->implicit(ValidatedVectorMultivariateFunction(List<ValidatedScalarMultivariateFunction>(1u,f)),ip,ExactBoxType(1u,ix));
    return res[0];
}

auto SolverBase::continuation(const ValidatedVectorMultivariateFunction& f,
                              const Vector<ApproximateNumericType>& p,
                              const ExactBoxType& ix,
                              const ExactBoxType& ip) const
    -> ValidatedVectorMultivariateFunctionModelType
{
    ARIADNE_NOT_IMPLEMENTED;
}




auto IntervalNewtonSolver::step(const ValidatedVectorMultivariateFunction& f,
                                const Vector<ValidatedNumericType>& x) const
    -> Vector<ValidatedNumericType>
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("Testing for root in "<<x);
    ARIADNE_LOG_PRINTLN_AT(1,"e="<<sup_error(x)<<", x="<<x);
    Vector<FloatDPValue> m(cast_exact(x));
    ARIADNE_LOG_PRINTLN_AT(1,"m="<<m);
    Vector<ValidatedNumericType> im(m);
    Vector<ValidatedNumericType> w=f.evaluate(im);
    ARIADNE_LOG_PRINTLN_AT(1,"f(m)="<<w);
    Matrix<ValidatedNumericType> A=f.jacobian(x);
    ARIADNE_LOG_PRINTLN_AT(1,"Df(r)="<<A);
    Matrix<ValidatedNumericType> Ainv=inverse(A);
    ARIADNE_LOG_PRINTLN_AT(1,"inverse(Df(r))="<<Ainv);
    Vector<ValidatedNumericType> dx=Ainv*w;
    ARIADNE_LOG_PRINTLN_AT(1,"dx="<<dx);
    Vector<ValidatedNumericType> nx= m - dx;
    ARIADNE_LOG_PRINTLN_AT(1,"nx="<<nx);
    return nx;
}

auto KrawczykSolver::step(const ValidatedVectorMultivariateFunction& f,
                          const Vector<ValidatedNumericType>& x) const
    -> Vector<ValidatedNumericType>
{
    ARIADNE_LOG_SCOPE_CREATE
    Matrix<ValidatedNumericType> I=Matrix<ValidatedNumericType>::identity(x.size(),x.zero_element());
    ARIADNE_LOG_PRINTLN("Testing for root in "<<x);
    ARIADNE_LOG_PRINTLN_AT(1,"e="<<sup_error(x)<<", x="<<x);
    Vector<FloatDPValue> m(cast_exact(x));
    ARIADNE_LOG_PRINTLN_AT(1,"m="<<m);
    Vector<ValidatedNumericType> im(m);
    Vector<ValidatedNumericType> fm=f.evaluate(im);
    ARIADNE_LOG_PRINTLN_AT(1,"f(m)="<<fm);
    Matrix<ValidatedNumericType> J=f.jacobian(x);
    ARIADNE_LOG_PRINTLN_AT(1,"Df(r)="<<J);
    Matrix<ValidatedNumericType> M=inverse(midpoint(J));
    ARIADNE_LOG_PRINTLN_AT(1,"inverse(Df(m))="<<M);
    Vector<ValidatedNumericType> dx=M*fm-(I-M*J)*(x-m);
    ARIADNE_LOG_PRINTLN_AT(1,"dx="<<dx);
    Vector<ValidatedNumericType> nx= m - dx;
    ARIADNE_LOG_PRINTLN_AT(1,"nx="<<nx);
    Vector<ValidatedNumericType> nr(nx);
    ARIADNE_LOG_PRINTLN_AT(1,"nr="<<nr);
    return nr;
}


auto FactoredKrawczykSolver::step(const ValidatedVectorMultivariateFunction& f,
                                  const Vector<ValidatedNumericType>& x) const
    -> Vector<ValidatedNumericType>
{
    ARIADNE_LOG_SCOPE_CREATE;
    Matrix<ValidatedNumericType> I=Matrix<ValidatedNumericType>::identity(x.size(),x.zero_element());
    ARIADNE_LOG_PRINTLN("Testing for root in "<<x);
    ARIADNE_LOG_PRINTLN_AT(1,"e="<<sup_error(x)<<", x="<<x);
    Vector<FloatDPValue> m(cast_exact(x));
    ARIADNE_LOG_PRINTLN_AT(1,"m="<<m);
    Vector<ValidatedNumericType> im(m);
    Vector<ValidatedNumericType> fm=f.evaluate(im);
    ARIADNE_LOG_PRINTLN_AT(1,"f(m)="<<fm);
    Matrix<ValidatedNumericType> J=f.jacobian(x);
    ARIADNE_LOG_PRINTLN_AT(1,"Df(r)="<<J);
    Matrix<ValidatedNumericType> mJ(midpoint(J));
    Matrix<ValidatedNumericType> M=inverse(mJ);
    ARIADNE_LOG_PRINTLN_AT(1,"inverse(Df(m))="<<M);
    Vector<ValidatedNumericType> dx=M*(fm+(J-mJ)*(x-m));
    ARIADNE_LOG_PRINTLN_AT(1,"dx="<<dx);
    Vector<ValidatedNumericType> nx= m - dx;
    ARIADNE_LOG_PRINTLN_AT(1,"nx="<<nx);
    Vector<ValidatedNumericType> nr(nx);
    ARIADNE_LOG_PRINTLN_AT(1,"nr="<<nr);
    return nr;
}


ValidatedVectorMultivariateFunctionModelDP
IntervalNewtonSolver::implicit_step(const ValidatedVectorMultivariateFunction& f,
                                    const ValidatedVectorMultivariateFunctionModelDP& id,
                                    const ValidatedVectorMultivariateFunctionModelDP& h) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    const Nat m=id.size();
    const Nat n=h.size();
    DP pr;

    ARIADNE_LOG_PRINTLN("f="<<f);
    ARIADNE_LOG_PRINTLN("h="<<h);

    ValidatedVectorMultivariateFunctionModelDP mh=h; mh.clobber();
    ARIADNE_LOG_PRINTLN("midpoint(h)="<<mh);

    ValidatedScalarMultivariateFunction zero_function(f.domain());
    Matrix<ValidatedScalarMultivariateFunction> D2f(n,n,zero_function);
    for(Nat i=0; i!=n; ++i) {
        for(Nat j=0; j!=n; ++j) {
            D2f[i][j]=f[i].derivative(m+j);
        }
    }
    ARIADNE_LOG_PRINTLN("D2f="<<D2f);

    ValidatedNumericType zero(0,pr);
    ValidatedScalarMultivariateFunctionModelDP z=h[0]*zero;
    ValidatedVectorMultivariateFunctionModelDP idh=join(id,h);

    Matrix<ValidatedScalarMultivariateFunctionModelDP> J(n,n,z);
    for(Nat i=0; i!=n; ++i) {
        for(Nat j=0; j!=n; ++j) {
            J[i][j]=compose(D2f[i][j],idh);
        }
    }
    ARIADNE_LOG_PRINTLN("J="<<J);

    Matrix<UpperIntervalType> rngJ(n,n,pr);
    for(Nat i=0; i!=n; ++i) {
        for(Nat j=0; j!=n; ++j) {
            UpperIntervalType D2fij=UpperIntervalType(unchecked_evaluate(D2f[i][j],cast_singleton(product(id.range(),h.range()))));
            rngJ[i][j]=intersection(J[i][j].range(),D2fij);
        }
    }
    ARIADNE_LOG_PRINTLN("rngJ="<<rngJ);

    ValidatedVectorMultivariateFunctionModelDP fidmh=compose(f,join(id,mh));
    ARIADNE_LOG_PRINTLN_AT(1,"compose(f,join(id,midpoint(h)))="<<fidmh);

    ValidatedVectorMultivariateFunctionModelDP dh(n,z);
    if(n==1) {
        if(possibly(contains(rngJ[0][0],FloatDPValue(0.0_x,dp)))) {
            ARIADNE_THROW(SingularJacobianException,"IntervalNewtonSolver","D2f(P,X)="<<rngJ[0][0]<<" which contains zero.");
        }
        if(possibly(contains(J[0][0].range(),FloatDPValue(0.0_x,dp)))) {
            dh[0]=fidmh[0]/cast_singleton(rngJ[0][0]);
        } else {
            dh[0]=fidmh[0]/J[0][0];
        }
    } else {
        dh=inverse(cast_singleton(rngJ))*fidmh;
    }
    ARIADNE_LOG_PRINTLN("dh="<<dh);

    ValidatedVectorMultivariateFunctionModelDP nh=mh-dh;
    ARIADNE_LOG_PRINTLN("nh="<<nh);
    return nh;
}


ValidatedVectorMultivariateFunctionModelDP
KrawczykSolver::implicit_step(const ValidatedVectorMultivariateFunction& f,
                              const ValidatedVectorMultivariateFunctionModelDP& p,
                              const ValidatedVectorMultivariateFunctionModelDP& x) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    const Nat np=p.size();
    const Nat nx=x.size();
    Matrix<ValidatedNumericType> I=Matrix<ValidatedNumericType>::identity(nx,dp);
    ARIADNE_LOG_PRINTLN("Contracting x="<<x);
    ARIADNE_LOG_PRINTLN("p="<<p);
    ARIADNE_LOG_PRINTLN("f="<<f);
    ValidatedVectorMultivariateFunctionModelDP mx(x);
    for(Nat i=0; i!=mx.size(); ++i) { mx[i].clobber(); }
    ARIADNE_LOG_PRINTLN("mx="<<mx);
    Vector<FloatDPError> ex(nx,dp);
    for(Nat i=0; i!=nx; ++i) { ex[i]=x[i].error(); }
    Vector<ValidatedNumericType> eix=make_bounds(ex);
    ARIADNE_LOG_PRINTLN("ex="<<ex);
    ValidatedVectorMultivariateFunctionModelDP fm=compose(f,join(p,mx));
    ARIADNE_LOG_PRINTLN_AT(1,"f(p,mx)="<<fm);
    Vector<ValidatedNumericType> rp(np,dp);
    for(Nat i=0; i!=np; ++i) { rp[i]=cast_singleton(p[i].range()); }
    Vector<ValidatedNumericType> rx(nx,dp);
    for(Nat i=0; i!=nx; ++i) { rx[i]=cast_singleton(x[i].range()); }
    Matrix<ValidatedNumericType> J=project(f.jacobian(join(rp,rx)),range(0,nx),range(np,np+nx));
    ARIADNE_LOG_PRINTLN("D2f(r)=J="<<J);
    Matrix<ValidatedNumericType> M=inverse(midpoint(J));
    ARIADNE_LOG_PRINTLN("inverse(D2f(m))=M="<<M);
    ARIADNE_LOG_PRINTLN("M*f(p,mx)="<<M*fm);
    ARIADNE_LOG_PRINTLN("(I-M*J)="<<(I-M*J));
    ARIADNE_LOG_PRINTLN("(I-M*J) * (ex*ValidatedNumericType(-1,+1))="<<(I-M*J)<<"*"<<eix<<"="<<(I-M*J) * eix);
    ValidatedVectorMultivariateFunctionModelDP dx= M*fm - (I-M*J) * eix;
    ARIADNE_LOG_PRINTLN("dx="<<dx);
    ValidatedVectorMultivariateFunctionModelDP nwx= mx - dx;
    ARIADNE_LOG_PRINTLN("nwx="<<nwx);
    return nwx;
}


Void IntervalNewtonSolver::_write(OutputStream& os) const
{
    os << "IntervalNewtonSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}

Void KrawczykSolver::_write(OutputStream& os) const
{
    os << "KrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}

Void FactoredKrawczykSolver::_write(OutputStream& os) const
{
    os << "FactoredKrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}


} // namespace Ariadne
