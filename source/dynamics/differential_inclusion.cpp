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

Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<ExactIntervalType>& bx) {
    return apply(f,Box<UpperIntervalType>(bx));
}

InclusionIntegratorBase::InclusionIntegratorBase(SweeperDP sweeper, StepSize step_size)
    : _sweeper(sweeper)
    , _step_size(step_size)
    , _number_of_steps_between_simplifications(8)
    , _number_of_variables_to_keep(4)
{
}

List<ValidatedVectorFunctionModelDP> InclusionIntegratorBase::flow(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType X0, Real tmax) const {
    //Solve the differential inclusion dot(x) in f(x)+V for x(0) in X0 up to time T.
    ARIADNE_LOG(2,"\nf:"<<f<<"\nV:"<<V<<"\nX0:"<<X0<<"\ntmax:"<<tmax<<"\n");

    // Ensure all arguments have the correct size;
    auto n=X0.size();
    DoublePrecision pr;
    assert(f.result_size()==n);
    assert(f.argument_size()==n);
    assert(V.size()==n);

    PositiveFloatDPValue hsug(this->_step_size);

    ValidatedVectorFunctionModelDP evolve_function = ValidatedVectorTaylorFunctionModelDP::identity(X0,this->_sweeper);
    auto t=PositiveFloatDPValue(0.0);

    List<ValidatedVectorFunctionModelDP> result;

    auto step = 0;
    while (possibly(t<FloatDPBounds(tmax,pr))) {
        ARIADNE_LOG(3,"step#:"<<step<<", t:"<<t<<", hsug:"<<hsug << "\n");
        if(possibly(t+hsug>FloatDPBounds(tmax,pr))) {  //FIXME: Check types for timing;
            hsug=cast_positive(cast_exact((tmax-t).upper()));
        }

        ARIADNE_LOG(4,"n. of extra parameters="<<evolve_function.argument_size()-n<<"\n");

        auto D = cast_exact_box(evolve_function.range());
        UpperBoxType B;
        PositiveFloatDPValue h;
        std::tie(h,B)=this->flow_bounds(f,V,D,hsug);
        ARIADNE_LOG(5,"h:"<<h<<", B:"<<B<<"\n");
        auto Phi = this->compute_step(f,V,D,h,B);
        ARIADNE_LOG(5,"Phi="<<Phi<<"\n");
        assert(Phi.domain()[n].upper()==h);
        PositiveFloatDPValue new_t=cast_positive(cast_exact((t+h).lower()));

        ValidatedVectorFunctionModelDP reach_function=compute_reach_function(evolve_function, Phi, t, new_t);
        ARIADNE_LOG(5,"reach_function="<<reach_function<<"\n");

        evolve_function=partial_evaluate(reach_function,n,new_t);
        ARIADNE_LOG(5,"evolve_function="<<evolve_function<<"\n");

        step+=1;

        if (step%this->_number_of_steps_between_simplifications==0) {
            this->_reconditioner->simplify(evolve_function);
            ARIADNE_LOG(5,"new_evolve_function="<<evolve_function<<"\n");
        }

        t=new_t;
        result.append(reach_function);
    }
    return result;
}

ValidatedVectorFunctionModelDP InclusionIntegratorBase::compute_reach_function(ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const {

    // Evolve function is xi(x,a) at s; Flow is phi(x,h,b)
    // Want (x,t,a,b):->phi(xi(x,a),t-s,b))

    SizeType n=evolve_function.result_size();

    SizeType p0=evolve_function.argument_size()-n;
    SizeType p1=Phi.argument_size()-(n+1);

    BoxDomainType Xi=evolve_function.domain()[range(0,n)];
    BoxDomainType A0=evolve_function.domain()[range(n,n+p0)];
    BoxDomainType A1=Phi.domain()[range(n+1,n+1+p1)];

    auto swp=this->_sweeper;
    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType DTA = join(Xi,Tau,A0,A1);
    ValidatedVectorTaylorFunctionModelDP xf=ValidatedVectorTaylorFunctionModelDP::projection(DTA,range(0,n),swp);
    ValidatedScalarTaylorFunctionModelDP tf=ValidatedScalarTaylorFunctionModelDP::coordinate(DTA,n,swp);
    ValidatedScalarTaylorFunctionModelDP hf=tf-t;
    ValidatedVectorTaylorFunctionModelDP a0f=ValidatedVectorTaylorFunctionModelDP::projection(DTA,range(n+1,n+1+p0),swp);
    ValidatedVectorTaylorFunctionModelDP a1f=ValidatedVectorTaylorFunctionModelDP::projection(DTA,range(n+1+p0,n+1+p0+p1),swp);

    ValidatedVectorTaylorFunctionModelDP ef=compose(evolve_function,join(xf,a0f));

    return compose(Phi,join(ef,hf,a1f));
}

Pair<PositiveFloatDPValue,UpperBoxType> InclusionIntegratorBase::flow_bounds(ValidatedVectorFunction f, UpperBoxType V, BoxDomainType D, PositiveFloatDPApproximation hsug) const {

    //! Compute a bound B for the differential inclusion dot(x) in f(x)+V for x(0) in D for step size h;
    ARIADNE_LOG(5,"D:"<<D);

    PositiveFloatDPValue h=cast_exact(hsug);
    UpperBoxType wD = D + (D-D.midpoint());
    UpperBoxType B = wD + 2*IntervalDomainType(0,h)*(apply(f,D)+V);

    while(not refines(D+h*(apply(f,B)+V),B)) {
        h=hlf(h);
    }
    for(auto i : range(0,4)) {
        B=D+IntervalDomainType(0,h)*(apply(f,B)+V);
    }
    ARIADNE_LOG(5,"h:" << h <<", B:"<< B << "\n");

    return std::make_pair(h,B);
}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPUpperBound>
InclusionIntegrator3rdOrder::
compute_norms(ValidatedVectorFunction const& f, UpperBoxType const& B) const {
    //! Compute the norms K=|f(B)|, L=|Df(B)|, H=|D2f(B)| and LN=l(Df(B));
    //  Estimate error terms;
    auto Df=f.differential(cast_singleton(B),2);
    ARIADNE_LOG(6,"  Df="<<Df<<"\n");
    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, L=ze, H=ze; FloatDPUpperBound LN=ze;
    for (auto i : range(f.result_size())) {
        auto Dfi=Df[i].expansion();
        FloatDPError Ki=ze, Li=ze, Hi=ze; FloatDPUpperBound LNi=ze;
        for (auto ac : Dfi) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
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

ValidatedVectorFunctionModelDP InclusionIntegrator3rdOrder::
compute_step(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType D, PositiveFloatDPValue h, UpperBoxType B) const {
    //! Compute a time-step for the differential inclusion dot(x) in f(x)+V for x(0) in D assuming bound B;
    DoublePrecision pr;
    auto n=D.size();
    FloatDPError K, L, H;
    FloatDPUpperBound LN;
    std::tie(K,L,H,LN)=this->compute_norms(f,B);
    auto KV=mag(norm(V));
    auto eLN = (possibly(LN>0)) ? FloatDPError(dexp(LN*h)) : FloatDPError(0u,pr);

    PositiveFloatDPValue c2(FloatDP(7.0/8,pr)); PositiveFloatDPBounds c1(c2/6);
    FloatDPError e = (c1*KV*H*(K+KV)+c2*KV*(L*L+H*(K+5u*KV/2u))*eLN)/cast_positive(1u-h*L/2u)*pow(h,3u);
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
    auto zf=ValidatedScalarTaylorFunctionModelDP(DHVAE,swp);
    auto xf=ValidatedVectorTaylorFunctionModelDP::projection(DHVAE,range(0,n),swp);
    auto tf=ValidatedScalarTaylorFunctionModelDP::coordinate(DHVAE,n,swp);
    auto vf=ValidatedVectorTaylorFunctionModelDP::projection(DHVAE,range(n+1,n+1+n),swp);
    auto af=ValidatedVectorTaylorFunctionModelDP::projection(DHVAE,range(n+1+n,n+1+2*n),swp);
    auto ef=ValidatedVectorTaylorFunctionModelDP::projection(DHVAE,range(n+1+2*n,n+1+3*n),swp);

    auto w=ValidatedVectorTaylorFunctionModelDP(n,DHVAE,swp);
    for (auto i : range(n)) { w[i]=vf[i]+af[i]*(tf-h/2)/h; }
    ARIADNE_LOG(6,"w:"<<w<<"\n");

    auto x0=ValidatedVectorTaylorFunctionModelDP(n,DHVAE,swp);
    for (auto i : range(n)) { x0[i]=xf[i]; }
    ARIADNE_LOG(6,"x0:"<<x0<<"\n");

    auto phi=ValidatedVectorTaylorFunctionModelDP(n,DHVAE,swp);
    for (auto i : range(n)) { phi[i]=ValidatedScalarTaylorFunctionModelDP(DHVAE,swp)+cast_singleton(B[i]); }
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

Tuple<FloatDPError,FloatDPError,FloatDPUpperBound> InclusionIntegrator2ndOrder::
compute_norms(ValidatedVectorFunction const& f, UpperBoxType const& B) const {
    //! Compute the norms K=|f(B)|, L=|Df(B)|, and LN=l(Df(B));
    //  Estimate error terms;
    auto Df=f.differential(cast_singleton(B),1);
    ARIADNE_LOG(6,"  Df="<<Df<<"\n");
    FloatDPError ze(0,dp);
    FloatDPError K=ze, L=ze; FloatDPUpperBound LN=ze;
    for (auto i : range(f.result_size())) {
        auto Dfi=Df[i].expansion();
        FloatDPError Ki=ze, Li=ze; FloatDPUpperBound LNi=ze;
        for (auto ac : Dfi) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
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

ValidatedVectorFunctionModelDP InclusionIntegrator2ndOrder::
compute_step(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType D, PositiveFloatDPValue h, UpperBoxType B) const {
    //! Compute a time-step for the differential inclusion dot(x) in f(x)+V for x(0) in D assuming bound B;
    auto n=D.size();
    FloatDPError K, L; FloatDPUpperBound LN;
    std::tie(K,L,LN)=this->compute_norms(f,B);
    PositiveFloatDPUpperBound KV=mag(norm(V));
    FloatDPError eLN = (possibly(LN>0)) ? FloatDPError(dexp(LN*h)) : FloatDPError(1u,dp);
    FloatDPError e = pow(h,2u)*(2u*KV*L*eLN);
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
    auto zf=ValidatedScalarTaylorFunctionModelDP(DHVE,swp);

    auto xf=ValidatedVectorTaylorFunctionModelDP::projection(DHVE,range(0,n),swp);
    auto tf=ValidatedScalarTaylorFunctionModelDP::coordinate(DHVE,n,swp);
    auto vf=ValidatedVectorTaylorFunctionModelDP::projection(DHVE,range(n+1,n+1+n),swp);
    auto ef=ValidatedVectorTaylorFunctionModelDP::projection(DHVE,range(n+1+n,n+1+2*n),swp);

    auto w=ValidatedVectorTaylorFunctionModelDP(n,DHVE,swp);

    for (auto i : range(n)) { w[i]=vf[i]; }
    ARIADNE_LOG(6,"w:"<<w);

    auto x0=ValidatedVectorTaylorFunctionModelDP(n,DHVE,swp);
    for (auto i : range(n)) { x0[i]=xf[i]; }
    ARIADNE_LOG(6,"x0:"<<x0);

    auto phi=ValidatedVectorTaylorFunctionModelDP(n,DHVE,swp);
    for (auto i : range(n)) { phi[i]=zf+cast_singleton(B[i]); }
    ARIADNE_LOG(6,"phi0:"<<phi);

    for (auto i : range(6)) {
        phi=antiderivative(compose(f,phi)+w,n)+x0;

    }
    ARIADNE_LOG(6,"(derivative(phi,n)-(compose(f,phi)+w)).range():" << (derivative(phi,n)-(compose(f,phi)+w)).range()<<"\n");

    for (auto i : range(n)) {
        phi[i]=phi[i]+ef[i];
    }
    return phi;
}



LohnerReconditioner::LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep)
    : _sweeper(sweeper), _number_of_variables_to_keep(number_of_variables_to_keep) {
    this->verbosity = 0;
}

ValidatedVectorFunctionModelDP LohnerReconditioner::expand_errors(ValidatedVectorFunctionModelDP Phi) const {
    BoxDomainType domain=Phi.domain();
    BoxDomainType errors=cast_exact(cast_exact(Phi.errors())*FloatDPUpperInterval(-1,+1)); // FIXME: Avoid cast;
    ARIADNE_LOG(6,"Uniform errors:"<<errors<<"\n");
    for(SizeType i=0; i!=Phi.result_size(); ++i) { Phi[i].set_error(0); }
    ValidatedVectorFunctionModelDP error_function=ValidatedVectorTaylorFunctionModelDP::identity(errors,this->_sweeper);
    return embed(Phi,errors)+embed(domain,error_function);
}

struct IndexedFloatDP
{
    SizeType index;
    FloatDP value;

    IndexedFloatDP() : index(0), value(FloatDP()) {}
};

OutputStream& operator<<(OutputStream& os, IndexedFloatDP const& ifl) {
    return os << "(" << ifl.index << ":" << ifl.value << ")"; }

struct IndexedFloatDPComparator
{
    inline bool operator() (const IndexedFloatDP& ifl1, const IndexedFloatDP& ifl2)
    {
        return (ifl1.value < ifl2.value);
    }
};

Void LohnerReconditioner::simplify(ValidatedVectorFunctionModelDP& phi) const {
    ARIADNE_LOG(6,"simplifying\n");
    ARIADNE_LOG(6,"phi="<<phi<<"\n");

    ValidatedVectorTaylorFunctionModelDP& tphi = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(phi.reference());

    auto m=phi.argument_size();
    auto n=phi.result_size();

    ARIADNE_LOG(6,"num.parameters="<<m<<", to keep="<< this->_number_of_variables_to_keep <<"\n");
    // Compute effect of error terms, but not of original variables;
    Matrix<FloatDP> C(m,n);
    for (auto i : range(n)) {
        auto p=tphi[i].model().expansion();

        for (auto ac : p) {
            MultiIndex const& a=ac.index();
            FloatDPValue& c=ac.coefficient();
            for (auto j : range(m)) {
                if (a[j]!=0) {
                    C[j][i] = C[j][i]+abs(c).raw();
                }
            }
        }
    }

    ARIADNE_LOG(6,"C"<<C<<"\n");

    Array<IndexedFloatDP> Ce(m);
    for (auto j : range(m)) {
        Ce[j].index = j;
        for (auto i : range(n)) {
            Ce[j].value += C[j][i];
        }
    }
    ARIADNE_LOG(6,"Ce:"<<Ce<<"\n");
    auto SCe=Ce;
    std::sort(SCe.begin(),SCe.end(),IndexedFloatDPComparator());
    ARIADNE_LOG(6,"SortedCe:"<<SCe<<"\n");
    List<SizeType> keep_indices;
    List<SizeType> remove_indices;
    int number_of_variables_to_remove = m - this->_number_of_variables_to_keep;
    ARIADNE_LOG(6, "Number of variables to remove:" << number_of_variables_to_remove<<"\n");
    for (int j : range(m)) {
        if (j < number_of_variables_to_remove) {
            remove_indices.append(SCe[j].index);
        } else {
            keep_indices.append(SCe[j].index);
        }
    }
    ARIADNE_LOG(6,"keep_indices:"<<keep_indices<<"\n");
    ARIADNE_LOG(6,"remove_indices:"<<remove_indices<<"\n");

    for (int i : range(n)) {
        ErrorType error = tphi[i].error();
        for(TaylorModel<ValidatedTag,FloatDP>::ConstIterator iter=tphi[i].model().begin(); iter!=tphi[i].model().end(); ++iter) {
            MultiIndex const& xa=iter->key();
            FloatDPValue const& xv=iter->data();
            for(SizeType k=0; k!=remove_indices.size(); ++k) {
                if(xa[remove_indices[k]]!=0) {
                    error += mag(xv);
                    break;
                }
            }
        }
        tphi[i].set_error(error);
    }

    auto old_domain=phi.domain();
    auto new_domain=BoxDomainType(Vector<IntervalDomainType>(keep_indices.size(),[&old_domain,&keep_indices](SizeType j){return old_domain[keep_indices[j]];}));
    auto projection=ValidatedVectorTaylorFunctionModelDP(m,new_domain,this->_sweeper);
    for (auto i : range(new_domain.size())) { projection[keep_indices[i]]=ValidatedScalarTaylorFunctionModelDP::coordinate(new_domain,i,this->_sweeper); }
    for (auto i : range(remove_indices.size())) {
        auto j=remove_indices[i]; auto cj=old_domain[j].midpoint();
        projection[j]=ValidatedScalarTaylorFunctionModelDP::constant(new_domain,cj,this->_sweeper); }
    phi=compose(phi,projection);

    phi = this->expand_errors(phi);
}

} // namespace Ariadne;


/*

#include "geometry/zonotope.hpp"

namespace Ariadne {

ValidatedVectorTaylorFunctionModelDP lohner_approximation(ValidatedVectorTaylorFunctionModelDP f) {
    auto n=f.result_size();
    auto models=f.models();
    DoublePrecision pr;
    PositiveFloatDPValue zero(pr);
    Vector<FloatDPValue> b=Vector<FloatDPValue>(n,zero);
    Vector<FloatDPError> e=Vector<FloatDPError>(n,zero);
    Matrix<FloatDPValue> A=Matrix<FloatDPValue>(n,models[0].argument_size(),zero);
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

    b=reinterpret_cast<Vector<FloatDPValue>const&>(z.centre());
    A=reinterpret_cast<Matrix<FloatDPValue>const&>(z.generators());
    e=reinterpret_cast<Vector<FloatDPError>const&>(z.error());
    auto p=z.number_of_generators();
    Vector<ValidatedTaylorModelDP> r(n,ValidatedTaylorModelDP(p,f.properties()));
    for (auto i : range(n)) {
        r[i].set_value(b[i]);
        for (auto j : range(p)) {
            r[i].set_gradient(j,A[i][j]);
        }
        r[i].set_error(e[i]);
    }

    return ValidatedVectorTaylorFunctionModelDP(BoxDomainType(n,IntervalDomainType(-1,+1)),r);
}



} // namespace Ariadne;

*/
