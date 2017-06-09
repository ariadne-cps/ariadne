/***************************************************************************
 *            differential_inclusion.cpp
 *
 *  Copyright  2008-17  Pieter Collins, Sanja Zivanovic
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

#include "differential_inclusion.hpp"
#include "function/taylor_function.hpp"
#include "solvers/integrator.hpp"

namespace Ariadne {

#define ARIADNE_LOG_PRINT(level, expr) { ARIADNE_LOG(level,#expr << "=" << (expr) << "\n"); }

using ValidatedScalarFunctionModel = ValidatedScalarFunctionModel64;
using ValidatedVectorFunctionModel = ValidatedVectorFunctionModel64;
using ValidatedScalarTaylorFunctionModel = ValidatedScalarTaylorFunctionModel64;
using ValidatedVectorTaylorFunctionModel = ValidatedVectorTaylorFunctionModel64;

Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<ExactIntervalType>& bx) {
    return apply(f,Box<UpperIntervalType>(bx));
}

InclusionIntegratorBase::InclusionIntegratorBase(Sweeper64 sweeper, StepSize step_size)
    : _reconditioner(new LohnerReconditioner(sweeper,number_of_variables_to_keep=4))
    , _sweeper(sweeper)
    , _step_size(step_size)
    , _number_of_steps_between_simplifications(8)
    , _number_of_variables_to_keep(4)
{
}

List<ValidatedVectorFunctionModel> InclusionIntegratorBase::flow(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType X0, Real tmax) const {
    //Solve the differential inclusion dot(x) in f(x)+V for x(0) in X0 up to time T.
    ARIADNE_LOG(2,"\nf:"<<f<<"\nV:"<<V<<"\nX0:"<<X0<<"\ntmax:"<<tmax<<"\n");

    // Ensure all arguments have the correct size;
    auto n=X0.size();
    Precision64 pr;
    assert(f.result_size()==n);
    assert(f.argument_size()==n);
    assert(V.size()==n);

    PositiveFloat64Value hsug(this->_step_size);

    ValidatedVectorFunctionModel evolve_function = ValidatedVectorTaylorFunctionModel::identity(X0,this->_sweeper);
    auto t=PositiveFloat64Value(0.0);

    List<ValidatedVectorFunctionModel> result;

    auto step = 0;
    while (possibly(t<Float64Bounds(tmax,pr))) {
        ARIADNE_LOG(3,"\n");
        ARIADNE_LOG(3,"step:"<<step<<", t:"<<t<<", hsug:"<<hsug);
        if(possibly(t+hsug>Float64Bounds(tmax,pr))) {  //FIXME: Check types for timing;
            hsug=cast_positive(cast_exact((tmax-t).upper()));
        }

        auto D = cast_exact(evolve_function.range());
        UpperBoxType B;
        PositiveFloat64Value h;
        std::tie(h,B)=this->flow_bounds(f,V,D,hsug);
        if(verbosity>=3) { std::clog << ", h:"<<h<<", B:"<<B<<"\n"; }
        auto Phi = this->compute_step(f,V,D,h,B);
        ARIADNE_LOG(7,"Phi="<<Phi<<"\n");
        assert(Phi.domain()[n].upper()==h);
        PositiveFloat64Value new_t=cast_positive(cast_exact((t+h).lower()));

        // Simplify terms in Phi
        ARIADNE_LOG(9,"Phi="<<Phi<<"\n");
        ValidatedVectorTaylorFunctionModel& TPhi = const_cast<ValidatedVectorTaylorFunctionModel&>(dynamic_cast<ValidatedVectorTaylorFunctionModel const&>(*Phi.raw_pointer()));
        ARIADNE_LOG(9,"TPhi="<<TPhi<<"\n");
        TPhi.set_properties(ThresholdSweeper<Float64>(pr,4e-3));
        TPhi.simplify();
        ARIADNE_LOG(9,"TPhi="<<TPhi<<"\n");
        ARIADNE_LOG(9,"Phi="<<Phi<<"\n");
        assert(Phi.domain()[n].upper()==h);
        ARIADNE_LOG(5,"evolve_function.domain()="<<evolve_function.domain()<<"\n");
        ARIADNE_LOG(5,"evolve_function.codomain()="<<evolve_function.codomain()<<"\n");
        ARIADNE_LOG(5,"Phi.domain()="<<Phi.domain()<<"\n");

        assert(evolve_function.result_size()==n);
        SizeType p0=evolve_function.argument_size()-n;
        SizeType p1=Phi.argument_size()-(n+1);

        BoxDomainType A0=evolve_function.domain()[range(n,n+p0)];
        BoxDomainType A1=Phi.domain()[range(n+1,n+1+p1)];

        auto Psi = partial_evaluate(Phi,n,Float64Value(h));
        ARIADNE_LOG(7,"Psi="<<Psi<<"\n");

        // Evolve function is xi(x,a) at s; Flow is phi(x,h,b)
        // Want (x,t,a,b):->phi(xi(x,a),t-s,b))
        auto swp=TPhi.properties();
        auto Tau=IntervalDomainType(t,new_t);
        ARIADNE_LOG(7,"Tau="<<Tau<<"\n");
        BoxDomainType DTA = join(X0,Tau,A0,A1);
        ARIADNE_LOG(7,"DTA="<<DTA<<"\n");
        ValidatedVectorTaylorFunctionModel xf=ValidatedVectorTaylorFunctionModel::projection(DTA,range(0,n),swp);
        ValidatedScalarTaylorFunctionModel tf=ValidatedScalarTaylorFunctionModel::coordinate(DTA,n,swp);
        ARIADNE_LOG(7,"tf="<<tf<<"\n");
        ValidatedScalarTaylorFunctionModel hf=tf-t;
        ARIADNE_LOG(7,"hf="<<hf<<"\n");
        ValidatedVectorTaylorFunctionModel a1f=ValidatedVectorTaylorFunctionModel::projection(DTA,range(n+1+p0,n+1+p0+p1),swp);
        ARIADNE_LOG(7,"a1f="<<join(xf,a1f)<<"\n");
        ValidatedVectorTaylorFunctionModel a0f=ValidatedVectorTaylorFunctionModel::projection(DTA,range(n+1,n+1+p0),swp);
        ARIADNE_LOG(7,"a0f="<<join(xf,a0f)<<"\n");
        ARIADNE_LOG(7,"join(xf,a0f)="<<join(xf,a0f)<<"\n");

        ValidatedVectorTaylorFunctionModel ef=compose(evolve_function,join(xf,a0f));
        ARIADNE_LOG(5,"ef.domain()"<<ef.domain()<<"\n");
        ARIADNE_LOG(5,"ef.range()"<<ef.range()<<"\n");

        ValidatedVectorFunctionModel reach_function=compose(Phi,join(ef,hf,a1f));
        ARIADNE_LOG(7,"reach_function="<<reach_function<<"\n");
        //reach_function=this->_reconditioner->expand_errors(reach_function);

        evolve_function=partial_evaluate(reach_function,n,new_t);
        ARIADNE_LOG(7,"evolve_function="<<evolve_function<<"\n");

        step+=1;
        if (step%this->_number_of_steps_between_simplifications==0) {
            this->_reconditioner->simplify(evolve_function);
        }

        ARIADNE_LOG(5,"new_evolve_function.range():"<<evolve_function.range()<<"\n");
        t=new_t;
        //reach_sets.append(lohner_approximation(reach_function));
        result.append(reach_function);

        ARIADNE_LOG(7,evolve_function<<"\n");
        ARIADNE_LOG(5,"new_evolve_function.errors():"<<evolve_function.errors()<<"\n");
    }
    return result;
}


Pair<PositiveFloat64Value,UpperBoxType> InclusionIntegratorBase::flow_bounds(ValidatedVectorFunction f, UpperBoxType V, BoxDomainType D, PositiveFloat64Approximation hsug) const {
    //! Compute a bound B for the differential inclusion dot(x) in f(x)+V for x(0) in D for step size h;
    ARIADNE_LOG(3,"D:"<<D);

    apply(f,D); //f(D); //image(f,D);

    PositiveFloat64Value h=cast_exact(hsug);
    UpperBoxType B = D + 2*IntervalDomainType(0,h)*(apply(f,D)+V);

    while(not refines(D+h*(apply(f,B)+V),B)) {
        h=hlf(h);
    }
    for(auto i : range(0,4)) {
        B=D+IntervalDomainType(0,h)*(apply(f,B)+V);
    }
    ARIADNE_LOG(3,"B:"<<B);
    return std::make_pair(h,B);
}

Tuple<Float64Error,Float64Error,Float64Error,Float64UpperBound>
InclusionIntegrator3rdOrder::
compute_norms(EffectiveVectorFunction const& f, UpperBoxType const& B) const {
    //! Compute the norms K=|f(B)|, L=|Df(B)|, H=|D2f(B)| and LN=l(Df(B));
    //  Estimate error terms;
    auto Df=f.differential(cast_singleton(B),2);
    ARIADNE_LOG(6,"  Df="<<Df<<"\n");
    Precision64 pr;
    Float64Error ze(pr);
    Float64Error K=ze, L=ze, H=ze; Float64UpperBound LN=ze;
    for (auto i : range(f.result_size())) {
        auto Dfi=Df[i].expansion();
        Float64Error Ki=ze, Li=ze, Hi=ze; Float64UpperBound LNi=ze;
        for (auto ac : Dfi) {
            MultiIndex const& a=ac.index();
            Float64Bounds const& c=ac.coefficient();
            if (a.degree()==0) {
                Ki += mag(c);
            } else if (a.degree()==1) {
                Li += mag(c);
                if (a[i]==1) { LNi += c.upper(); }
                else { LNi += mag(c); }
            } else {
                assert(a.degree()==2);
                Hi += mag(c);
            }
        }
        K=max(K,Ki); L=max(L,Li); H=max(H,Hi); LN=max(LN,LNi);
    }
    return std::tie(K,L,H,LN);
}

ValidatedVectorFunctionModel InclusionIntegrator3rdOrder::
compute_step(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType D, PositiveFloat64Value h, UpperBoxType B) const {
    //! Compute a time-step for the differential inclusion dot(x) in f(x)+V for x(0) in D assuming bound B;
    Precision64 pr;
    auto n=D.size();
    Float64Error K, L, H;
    Float64UpperBound LN;
    std::tie(K,L,H,LN)=this->compute_norms(f,B);
    auto KV=mag(norm(V));
    auto eLN = (possibly(LN>0)) ? Float64Error(dexp(LN*h)) : Float64Error(0u,pr);

    PositiveFloat64Value c2(Float64(7.0/8,pr)); PositiveFloat64Bounds c1(c2/6);
    Float64Error e = (c1*KV*H*(K+KV)+c2*KV*(L*L+H*(K+5u*KV/2u))*eLN)/cast_positive(1u-h*L/2u)*pow(h,3u);
    ARIADNE_LOG(6,"e:"<<e<<", h:"<<h<<", e/h^3:"<<e/pow(h,3u)<<"\n");
    ARIADNE_LOG(6,"  K:"<<K<<", KV:"<<KV<<", L:"<<L<<", LN:"<<LN<<", eLN:"<<eLN<<", H:"<<H<<", e:"<<e<<"\n");

    IntervalDomainType Ht=IntervalDomainType(-h,+h);
    BoxDomainType A=cast_exact_box(3*V);
    BoxDomainType E=error_domain(n,e);

    //  Set up estimate for differential equation;
    //  Use affine wi=ai0+(t-h/2)ai1/h;
    //    with |ai0|<=Vi, |ai1|<=3Vi, |wi(t)|<=5Vi/2, and |w'(t)|<=3Vi/h.;
    //  Alternatively use step inputs wi=ai,0 for t<h/2 and wi=ai,1 for t>h/2;
    //    with |ai0|,|ai1|<= 2Vi , and |wi(t)|<=2Vi.;

    //  The flow is a function of n state variables, 2n parameter variables and 1 time variable;
    //  Assume for now noise in all variables;
    //  We have dot[xi](t)=fi(x(t))+wi(t);
    auto swp=this->_sweeper;
    auto DHVAE=join(D,Ht,V,A,E);
    ARIADNE_LOG(6,"DHVAE:"<<DHVAE<<"\n");
    auto zf=ValidatedScalarTaylorFunctionModel(DHVAE,swp);
    auto xf=ValidatedVectorTaylorFunctionModel::projection(DHVAE,range(0,n),swp);
    auto tf=ValidatedScalarTaylorFunctionModel::coordinate(DHVAE,n,swp);
    auto vf=ValidatedVectorTaylorFunctionModel::projection(DHVAE,range(n+1,n+1+n),swp);
    auto af=ValidatedVectorTaylorFunctionModel::projection(DHVAE,range(n+1+n,n+1+2*n),swp);
    auto ef=ValidatedVectorTaylorFunctionModel::projection(DHVAE,range(n+1+2*n,n+1+3*n),swp);

    auto w=ValidatedVectorTaylorFunctionModel(n,DHVAE,swp);
    ARIADNE_LOG(6,"w:"<<w<<"\n");
    //for (auto i : range(n)) { w[i]=xat[n+i]+xat[2*n+i]*(xat[4*n]-h/2)/h; }
    for (auto i : range(n)) { w[i]=vf[i]+af[i]*(tf-h/2)/h; }
    ARIADNE_LOG(6,"w:"<<w<<"\n");

    auto x0=ValidatedVectorTaylorFunctionModel(n,DHVAE,swp);
    for (auto i : range(n)) { x0[i]=xf[i]; }
    ARIADNE_LOG(6,"x0:"<<x0<<"\n");

    auto phi=ValidatedVectorTaylorFunctionModel(n,DHVAE,swp);
    for (auto i : range(n)) { phi[i]=ValidatedScalarTaylorFunctionModel(DHVAE,swp)+cast_singleton(B[i]); }
    ARIADNE_LOG(6,"phi0:"<<phi<<"\n");

    for (auto i : range(6)) {
        phi=antiderivative(compose(f,phi)+w,n)+x0;
    }
    ARIADNE_LOG(6,"phi.errors():"<<phi.errors()<<"\n");
    ARIADNE_LOG(7,(derivative(phi,n)-(compose(f,phi)+w)).range()<<"\n");

    for (auto i : range(n)) {
        phi[i]=phi[i]+ef[i];
    }

    return phi;
}

Tuple<Float64Error,Float64Error,Float64UpperBound> InclusionIntegrator2ndOrder::
compute_norms(EffectiveVectorFunction const& f, UpperBoxType const& B) const {
    //! Compute the norms K=|f(B)|, L=|Df(B)|, and LN=l(Df(B));
    //  Estimate error terms;
    auto Df=f.differential(cast_singleton(B),1);
    ARIADNE_LOG(6,"  Df="<<Df<<"\n");
    Float64Error ze(0,Precision64());
    Float64Error K=ze, L=ze; Float64UpperBound LN=ze;
    for (auto i : range(f.result_size())) {
        auto Dfi=Df[i].expansion();
        Float64Error Ki=ze, Li=ze; Float64UpperBound LNi=ze;
        for (auto ac : Dfi) {
            MultiIndex const& a=ac.index();
            Float64Bounds const& c=ac.coefficient();
            if (a.degree()==0) {
                Ki += mag(c);
            } else if (a.degree()==1) {
                Li += mag(c);
                if (a[i]==1) { LNi += c.upper(); }
                else { LNi += mag(c); }
            }
        }
        K=max(K,Ki); L=max(L,Li); LN=max(LN,LNi);
    }
    return std::tie(K,L,LN);
}

ValidatedVectorFunctionModel InclusionIntegrator2ndOrder::
compute_step(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType D, PositiveFloat64Value h, UpperBoxType B) const {
    //! Compute a time-step for the differential inclusion dot(x) in f(x)+V for x(0) in D assuming bound B;
    auto n=D.size();
    Float64Error K, L; Float64UpperBound LN;
    std::tie(K,L,LN)=this->compute_norms(f,B);
    PositiveFloat64UpperBound KV=mag(norm(V));
    Float64Error eLN = (possibly(LN>0)) ? Float64Error(dexp(LN*h)) : Float64Error(1u,Precision64());
    Float64Error e = pow(h,2u)*(2u*KV*L*eLN);
    ARIADNE_LOG(6,"e:"<<e<<", h:"<<h<<", e/h^2:"<<e/pow(h,2u)<<"\n");

    ARIADNE_LOG(6,"  K:"<<K<<", KV:"<<KV<<", L:"<<L<<", LN:"<<LN<<", eLN:"<<eLN<<", e:"<<e<<"\n");

    auto E=error_domain(n,e);
    auto Ht=IntervalDomainType(-h,+h);

    //  Set up estimate for differential equation;
    //  Use affine wi=ai0+(t-h/2)ai1/h;
    //    with |ai0|<=Vi, |ai1|<=3Vi, |wi(t)|<=5Vi/2, and |w'(t)|<=3Vi/h.;
    //  Alternatively use step inputs wi=ai,0 for t<h/2 and wi=ai,1 for t>h/2;
    //    with |ai0|,|ai1|<= 2Vi , and |wi(t)|<=2Vi.;

    //  The flow is a function of n state variables, 2n parameter variables and 1 time variable;
    //  Assume for now noise in all variables;
    //  We have dot[xi](t)=fi(x(t))+wi(t);
    auto swp=this->_sweeper;
    auto DHVE=join(D,Ht,V,E);
    // ARIADNE_LOG(6,"DHV["+str(DHV.size())+"]:"<<DHV);
    auto zf=ValidatedScalarTaylorFunctionModel(DHVE,swp);

    auto xf=ValidatedVectorTaylorFunctionModel::projection(DHVE,range(0,n),swp);
    auto tf=ValidatedScalarTaylorFunctionModel::coordinate(DHVE,n,swp);
    auto vf=ValidatedVectorTaylorFunctionModel::projection(DHVE,range(n+1,n+1+n),swp);
    auto ef=ValidatedVectorTaylorFunctionModel::projection(DHVE,range(n+1+n,n+1+2*n),swp);

    auto w=ValidatedVectorTaylorFunctionModel(n,DHVE,swp);

    for (auto i : range(n)) { w[i]=vf[i]; }
    // ARIADNE_LOG(6,"w:"<<w);

    auto x0=ValidatedVectorTaylorFunctionModel(n,DHVE,swp);
    for (auto i : range(n)) { x0[i]=xf[i]; }
    // ARIADNE_LOG(6,"x0:"<<x0);

    auto phi=ValidatedVectorTaylorFunctionModel(n,DHVE,swp);
    for (auto i : range(n)) { phi[i]=zf+cast_singleton(B[i]); }
    // ARIADNE_LOG(6,"phi0:"<<phi);

    for (auto i : range(6)) {
        phi=antiderivative(compose(f,phi)+w,n)+x0;
        // print i,phi.errors();
    }
    ARIADNE_LOG(6,(derivative(phi,n)-(compose(f,phi)+w)).range()<<"\n");

    for (auto i : range(n)) {
        phi[i]=phi[i]+ef[i];
    }
    return phi;
}



LohnerReconditioner::LohnerReconditioner(Sweeper64 sweeper, NumberOfVariablesToKeep number_of_variables_to_keep)
    : _sweeper(sweeper), _number_of_variables_to_keep(number_of_variables_to_keep) {
}

ValidatedVectorFunctionModel LohnerReconditioner::expand_errors(ValidatedVectorFunctionModel Phi) const {
    BoxDomainType domain=Phi.domain();
    BoxDomainType errors=cast_exact(cast_exact(Phi.errors())*Float64UpperInterval(-1,+1)); // FIXME: Avoid cast;
    ARIADNE_LOG(3,"Uniform errors:"<<errors);
    for(SizeType i=0; i!=Phi.result_size(); ++i) { Phi[i].set_error(0); }
    ValidatedVectorFunctionModel error_function=ValidatedVectorTaylorFunctionModel::identity(errors,this->_sweeper);
    return embed(Phi,errors)+embed(domain,error_function);
}

Void LohnerReconditioner::simplify(ValidatedVectorFunctionModel& phi) const {
    ARIADNE_LOG(4,"simplifying\n");
    ARIADNE_LOG(6,"phi="<<phi<<"\n");

    ValidatedVectorTaylorFunctionModel& tphi = dynamic_cast<ValidatedVectorTaylorFunctionModel&>(phi.reference());

    auto m=phi.argument_size();
    auto n=phi.result_size();
    // Compute effect of error terms, but not of original variables;
    Matrix<Float64Approximation> C(m,n);
    for (auto i : range(n)) {
        auto p=tphi[i].model().expansion();

        for (auto ac : p) {
            MultiIndex const& a=ac.index();
            Float64Value& c=ac.coefficient();
            for (auto j : range(m)) {
                if (a[j]!=0) {
                    C[j][i] = C[j][i]+abs(c);
                }
            }
        }
    }

    Vector<Float64Error> e(n,[&phi](SizeType i){return phi[i].error();});
    ARIADNE_LOG(3,"C"<<C<<", e"<<e<<"\n");

    List<Float64> Ce(m,Float64(Precision64()));
    for (auto j : range(m)) {
        for (auto i : range(n)) {
            Ce[j] += C[j][i].raw();
        }
    }
    ARIADNE_LOG(3,"Ce:"<<Ce<<"\n");
    auto SCe=Ce;
    std::sort(SCe.begin(),SCe.end());
    ARIADNE_LOG(3,"SortCe:"<<SCe<<"\n");
    List<SizeType> keep_indices;
    List<SizeType> remove_indices;
    SizeType number_of_variables_to_keep=this->_number_of_variables_to_keep;
    if (m<this->_number_of_variables_to_keep) { number_of_variables_to_keep=m; }
    auto threshold = (SCe[-number_of_variables_to_keep]+SCe[1-number_of_variables_to_keep])/2;
    for (auto j : range(m)) {
        if (Ce[j] < threshold) {
            remove_indices.append(j);
        }
        else {
            keep_indices.append(j);
        }
    }
    ARIADNE_LOG(3,"keep_indices:"<<keep_indices<<"\n");
    ARIADNE_LOG(3,"remove_indices:"<<remove_indices<<"\n");

    auto old_domain=phi.domain();
    auto new_domain=BoxDomainType(Vector<IntervalDomainType>(keep_indices.size(),[&old_domain,&keep_indices](SizeType j){return old_domain[keep_indices[j]];}));
    auto projection=ValidatedVectorTaylorFunctionModel(m,new_domain,this->_sweeper);
    for (auto i : range(new_domain.size())) { projection[keep_indices[i]]=ValidatedScalarTaylorFunctionModel::coordinate(new_domain,i,this->_sweeper); }
    for (auto i : range(remove_indices.size())) {
        auto j=remove_indices[i]; auto cj=cast_singleton(old_domain[j]);
        projection[j]=ValidatedScalarTaylorFunctionModel::constant(new_domain,cj,this->_sweeper); }
    phi=compose(phi,projection);

    phi = this->expand_errors(phi);
}

} // namespace Ariadne;


/*

#include "geometry/zonotope.hpp"

namespace Ariadne {

ValidatedVectorTaylorFunctionModel lohner_approximation(ValidatedVectorTaylorFunctionModel f) {
    auto n=f.result_size();
    auto models=f.models();
    Precision64 pr;
    PositiveFloat64Value zero(pr);
    Vector<Float64Value> b=Vector<Float64Value>(n,zero);
    Vector<Float64Error> e=Vector<Float64Error>(n,zero);
    Matrix<Float64Value> A=Matrix<Float64Value>(n,models[0].argument_size(),zero);
    for (auto i : range(n)) {
        b[i]=models[i].value();
        for (auto j : range(models[0].argument_size())) {
            A[i][j]=models[i].gradient_value(j);
        }
        e[i]=models[i].error();
    }
    auto z=Zonotope(b,A,e);
    // print z.error();
    z=orthogonal_approximation(z);

    b=reinterpret_cast<Vector<Float64Value>const&>(z.centre());
    A=reinterpret_cast<Matrix<Float64Value>const&>(z.generators());
    e=reinterpret_cast<Vector<Float64Error>const&>(z.error());
    auto p=z.number_of_generators();
    Vector<ValidatedTaylorModel64> r(n,ValidatedTaylorModel64(p,f.properties()));
    for (auto i : range(n)) {
        r[i].set_value(b[i]);
        for (auto j : range(p)) {
            r[i].set_gradient(j,A[i][j]);
        }
        r[i].set_error(e[i]);
    }

    return ValidatedVectorTaylorFunctionModel(BoxDomainType(n,IntervalDomainType(-1,+1)),r);
}



} // namespace Ariadne;

*/
