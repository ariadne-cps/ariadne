#include "../test/test.h"

#include "vector.h"
#include "polynomial.h"
#include "taylor_model.h"
#include "taylor_function.h"

namespace Ariadne {

Vector<TaylorModel>
my_flow(const Vector<TaylorModel>& vf, const Vector<TaylorModel>& yz, uint order)
{
    uint n=vf.size();
    assert(vf[0].argument_size()==n);
    assert(yz[0].argument_size()==n+1);

    Vector<TaylorModel> y(n,TaylorModel(n+1));
    Vector<TaylorModel> new_y(n,TaylorModel(n+1));

    // Set initial bound for iteration
    for(uint i=0; i!=y.size(); ++i) {
        y[i]=TaylorModel::constant(n+1,yz[i].range()+vf[i].range()*Interval(-1,+1));
    }


    for(uint j=0; j!=order; ++j) {
        new_y=compose(vf,y);
        for(uint i=0; i!=n; ++i) {
            new_y[i].antidifferentiate(n);
            new_y[i]+=yz[i];
            try { new_y[i]=intersection(y[i],new_y[i]); }
            catch(const IntersectionException& e) { std::cerr<<"Warning: "<<e.what()<<"\n"; }
            new_y[i].swap(y[i]);
        }
    }

    for(uint i=0; i!=n; ++i) { y[i].sweep(0.0); }
    //if(h.l==0) { y=split(y,n,true); }
    return y;
}

TaylorFunction
my_flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& h, const uint o)
{
    uint n=vf.size();
    const Vector<Interval>& b=vf.domain();

    std::cerr<<"vf="<<vf<<"\nd="<<d<<"\nh="<<h<<"\n"<<std::endl;
    assert(h.l==-h.u);

    // Scale multiply models by inverse reciprocal of radius
    Vector<TaylorModel> unscaled_vf(vf.size());
    for(uint i=0; i!=vf.size(); ++i) {
        unscaled_vf[i]=vf.models()[i]/rad_ivl(vf.domain()[i]); }
    std::cerr<<"unscaled_vf="<<unscaled_vf<<std::endl;

    Vector<TaylorModel> y0=TaylorModel::scalings(d);
    Vector<TaylorModel> unscaled_y0=unscale(y0,b);
    unscaled_y0=embed(unscaled_y0,1u);
    std::cerr<<"unscaled_y0"<<unscaled_y0<<std::endl;

    Vector<TaylorModel> hscaled_vf=unscaled_vf*h.u;
    std::cerr<<"hscaled_vf"<<hscaled_vf<<std::endl;

    Vector<TaylorModel> xscaled_flow=my_flow(hscaled_vf,unscaled_y0,o);
    std::cerr<<"xscaled_flow="<<xscaled_flow<<std::endl;

    Vector<TaylorModel> model_flow=scale(xscaled_flow,b);
    std::cerr<<"model_flow="<<model_flow<<std::endl;

    TaylorFunction flow(join(d,h),model_flow);
    std::cerr<<"\nflow="<<flow<<"\n"<<std::endl;

    return flow;
}

} // namespace Ariadne

using namespace Ariadne;

Vector<Float> e(uint n, uint i) { return Vector<Float>::unit(n,i); }
Expansion<Float> v(uint n, uint j) { return Expansion<Float>::variable(n,j); }
Polynomial<Float> p(uint n, uint j) { return Polynomial<Float>::variable(n,j); }
TaylorModel t(uint n, uint j) { return TaylorModel::variable(n,j); }
//TaylorVariable t(Vector<Interval> d, uint j) { return TaylorVariable::variable(d,j); }

int main() {
            // Vector field dtx=1.5
        Vector< Polynomial< Float > > vector_field_poly = (1.5+p(1,0)*0.0)*e(1,0);
        Vector<Interval> bounding_box(1, Interval(0.25,0.75));
        Vector<Interval> initial_box(1, Interval(0.375,0.50));
        Float step_size=1./32;
        Interval time_interval(-step_size,+step_size);
        uint order=2;
        TaylorFunction vector_field(bounding_box,vector_field_poly);
        TaylorFunction flow_model=my_flow(vector_field,initial_box,time_interval,order);
        //Vector< Polynomial<Float> > flow_poly = (p(2,0)+1.5*p(2,1))*e(1,0);
        //ARIADNE_TEST_PRINT(vector_field_poly)
        //ARIADNE_TEST_PRINT(flow_model.domain());
        //ARIADNE_TEST_PRINT(flow_model.models());
        //ARIADNE_TEST_EQUAL(flow_model.domain(),join(initial_box,time_interval));
        //ARIADNE_TEST_EQUAL(flow_model.polynomial(),flow_poly);
}