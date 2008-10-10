#include <iostream>
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "sparse_differential.h"
#include "function.h"
#include "approximate_taylor_model.h"

using std::cout; using std::endl; using std::flush;
using namespace Ariadne;

struct Henon {
  static const int result_size=2;
  static const int argument_size=2;
  static const int smoothness=2;
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=p[0]-x[0]*x[0]-p[1]*x[1];
    r[1]=x[0];
  }
};

struct HenonSquared {
  static const int result_size=2;
  static const int argument_size=2;
  static const int smoothness=2;
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[1]=p[0]-x[0]*x[0]-p[1]*x[1];
    r[0]=p[0]-r[1]*r[1]-p[1]*x[0];
  }
};

struct VanDerPol {
  static const int result_size=2;
  static const int argument_size=2;
  static const int smoothness=2;
  template<class R, class A, class P>
  void compute(R& r, const A& x, const P& p) const {
    r[0]=x[1];
    r[1]=p[0]*(1-x[0]*x[0])*x[1]-x[0];
  }
};


int main() {
  {
    // Slicing
    double ax[9]={2,3,5,7,11,13,17,19,23};
    Vector<Float> v(9,ax);
    Slice s=slice(1,2,3);
    Range r(2,5);
    cout << v << "\n" << project(v,s) << "\n";
    project(v,s)=project(v,r);
    cout << v << "\n";
    cout << endl;

    Vector<Float> w(4); cout << w << endl;

    Matrix<Float> A(3,3,ax);
    Slice sr=slice(1,1,2);
    Slice sc=slice(0,2,2);
    cout << A << "\n" << project(A,sr,sc) << endl;
  }

  // Initialise global variables
  double da[4]={-2,2,-2,2};
  double ca[2]={0,0};
  double ea[20]={1.5, 0,-0.375, -1,0,0, 0,0,0,0,  0, 1,0, 0,0,0,  0,0,0,0};
  ApproximateTaylorModel f(box(2,da),point(2,ca),SparseDifferentialVector<Float>(2,2,3,ea));

  cout << f << std::endl;
  cout << f.evaluate(Vector<Float>(2,0.0)) << std::endl;
  cout << f.jacobian(Vector<Float>(2,0.0)) << std::endl;


  {
    // Evaluation
    double xa[2]={0.1,0.2};
    double ixa[4]={0.05,0.15,0.15,0.25};
    Vector<Float> x=point(2,xa);
    Vector<Interval> ix=box(2,ixa);
    cout << "evaluate:\nf("<<x<<")="<<f.evaluate(x)<<"\nf("<<ix<<")="<<f.evaluate(ix)<<"\n\n";
    cout << "jacobian:\nDf("<<x<<")="<<f.jacobian(x)<<"\nDf("<<ix<<")="<<f.jacobian(ix)<<"\n\n";
  }

  {
    // Solution of a system of equations
    double ya[2]={1.6,0.5};
    Vector<Float> y(2,ya);
    Vector<Interval> x=solve(f,y);
    cout << "solve:\nf="<< f << "\ny=" << y << "\nx=" << x << "\nf(x)=" << f.evaluate(x) << "\n\n";
    cout << "f(x)="<< f.evaluate(midpoint(x))<<"\n\n";
    
    Vector<Float> p(2); p[0]=0.5; p[1]=-0.1/0.375;
    cout << "f.evaluate("<<p<<")="<<f.evaluate(p)<<"\n\n";

  }


  {
    // Implicit function
    double fa[20]={0, 1,2,1, 0.3,0.4,0.1,0.0,0.0,0.0,  0, 0,1,1, 0.1,0.2,0.3,0.4,0.5,0.6};
    Vector<Interval> d(3,Interval(-1,1));
    Vector<Float> c(3,0.0);
    ApproximateTaylorModel f(d,c,SparseDifferentialVector<Float>(2,3,2,fa));
    
    c=Vector<Float>(2,0.0);
    Vector<Float> v=Vector<Float>(1,1.0);
    cout << "implicit:\nf="<< f << "\nc=" << c << "\nv=" << v << "\n" << flush;
    ApproximateTaylorModel h=implicit(f);
    cout << "implicit(f)=" << h << "\n\n";
    ApproximateTaylorModel id=ApproximateTaylorModel::identity(project(d,range(0,1)));
    cout << "id=" << id << "\n\n";
    ApproximateTaylorModel g=join(id,h);
    cout << "(id,h)=" << g << endl;
    cout << "compose(f,(id,h))=f(x,h(x))=" << compose(f,g) << endl;
  }
    
  {
    // Flow
    cout << "Scalar flow:"<<endl;
    Vector<Interval> d(1,Interval(-2,2));
    Vector<Float> c(1,0.0);
    double ea[]={0,1,0,0,0,0,0,0};
    SparseDifferentialVector<Float> e(1,1,6,ea); 
    ApproximateTaylorModel am(d,c,e);
    cout << "am="<<am<<endl;
    ApproximateTaylorModel afm=flow(am);
    cout << "afm="<<afm<<endl;
  }

  {
    // Flow
    cout << "Flow:"<<endl;
    Vector<Interval> d(2,Interval(-2,2));
    Vector<Float> c(2,0.0);
    Function<VanDerPol> vdp(Vector<Float>(1,0.25));
    ApproximateTaylorModel vdpm(d,c,vdp,4,1);
    cout << "vdpm="<<vdpm<<endl;
    ApproximateTaylorModel vdpfm=flow(vdpm);
    cout << "vdpfm="<<vdpfm<<endl;
  }

  {
    // Hitting
    cout << "Hitting:"<<endl;
    uint n=1;
    uint d=5;
    Vector<Interval> fd(n,Interval(-2,2));
    Vector<Float> fc(n,0.0);
    double fea[]={1,1,0,0,0,0,0,0};
    SparseDifferentialVector<Float> fe(n,n,d,fea); 
    ApproximateTaylorModel fm(fd,fc,fe);
    cout << "fm=" << fm << endl;
    double gea[]={1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Vector<Interval> gd(n,Interval(-2,2));
    Vector<Float> gc(n,0.0);
    SparseDifferentialVector<Float> ge(1,n,d,gea); 
    ApproximateTaylorModel gm(gd,gc,ge);
    cout << "gm=" << gm << endl;
    const ApproximateTaylorModel hm=hitting(fm,gm);
    cout << "hm=" << hm << endl;
    Slice slc=slice(0,1,1);
    ApproximateTaylorModel hmp=project_model(hm,slc);
    cout << "hmp=" << hmp << endl;
    ApproximateTaylorModel zm=compose(gm,hmp);
    cout << "zm=" << zm << endl;
  }

  return 0;

  {
    // Test constructor from function
    uint order = 3;
    uint smoothness = 2;
    Float apm[2]={1.5,0.375};
    Float apt[2]={1.125,1.25};
    Float abx[4]={0.875,1.375,1.125,1.375};
    Vector<Float> pm=point(2,apm);
    Vector<Float> pt0=point(2,apt);
    Vector<Interval> bx=box(2,abx);

    Function<Henon> h(pm);
    Function<HenonSquared> hsq(pm);

    cout << h.expansion(pt0,3) << endl << endl;

    ApproximateTaylorModel tm1(bx,pt0,h,order,smoothness);
    cout << tm1 << endl;
    
    Vector<Float> pt1=tm1.evaluate(pt0);
    ApproximateTaylorModel tm2(bx,pt1,h,order,smoothness);
    cout << tm2 << endl;
   
    cout << tm1.evaluate(pt0) << endl << endl;
    cout << compose(tm2,tm1) << endl << endl;
    cout << compose(tm1,tm1) << endl << endl;

    ApproximateTaylorModel tmsq(bx,pt0,hsq,order,smoothness);
    cout << tmsq << endl << endl;

    Vector<Float> tv(2); tv[0]=0.03125; tv[1]=0.0625;
    ApproximateTaylorModel tm1tr(bx,Vector<Float>(pt0-tv),h,order,smoothness);
    cout << tm1 << tm1tr << endl;

    cout << tm1.evaluate(pt0) << tm1tr.evaluate(pt0) << endl;
  }   

  
}
