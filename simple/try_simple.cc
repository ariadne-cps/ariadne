#include "simple.h"

#include <iomanip>

std::string demangle(std::string s) {
    char r[100];
    SizeType i=0; SizeType j=9;
    while(s[j]!=0) {
        r[i]=s[j];
        ++i; ++j;
    }
    r[i]=0;
    return std::string(r);
}

template<class T> struct Virtual { ~Virtual() = default; };
template<class T> std::string mangled_class_name() { return typeid(Virtual<T>()).name(); }
template<class T> std::string class_name() { return demangle(mangled_class_name<T>()); }

template<class L> void test_logical_operations(L l) {
    std::cout<<"test_logical_operations(" << class_name<L>() << ")\n";
    L r(l); r=l;
    //r=l&&l; r=l||l; r=!(!l);
    r=conjunction(l,l); r=disjunction(l,l); r=negation(negation(l));
}

template<class X> void test_ring_operations(X x) {
    std::cout<<"test_ring_operations(" << class_name<X>() << ")\n";
    X r=x; r=add(x,x); r=mul(x,x); r=neg(x);
}

template<class X> void test_field_operations(X x) {
    std::cout<<"test_field_operations(" << class_name<X>()  << ")\n";
    X r(x); r=x; r=add(x,x); r=mul(x,x); r=neg(x); x=rec(x);
}

template<class X> void test_elementary_operations(X x) {
    std::cout<<"test_elementary_operations(" << class_name<X>()  << ")\n";
    X r(x); r=x; r=add(x,x); r=mul(x,x); r=neg(x); x=rec(x); x=exp(x); x=log(x);
}

template<class A, class X=NumericType<A>> void test_algebra_operations(A a, X x) {
    std::cout<<"test_algebra_operations(" << class_name<A>()  << "," << class_name<X>() << ")\n";
    A r=a;
    r=a; r=add(a,a); r=mul(a,a); r=neg(a); a=rec(a); a=exp(a); a=log(a);
    r=x; r=add(a,x); r=add(x,a); r=mul(a,x); r=mul(x,a);
}

template<class X1, class X2> void test_mixed_operations(X1 x1, X2 x2) {
    std::cout<<"test_mixed_operations(" << class_name<X1>()  << "," << class_name<X2>() << ")\n";
    auto ra=add(x1,x2); ra=add(x2,x1); auto rm=mul(x1,x2); rm=mul(x2,x1);
};

#define PRINT(x) { std::cout << #x << "=" << x << std::endl; }

void wrap();

void test() {
    Boolean bl=true;
    test_logical_operations(bl);
    Kleenean kl=true;
    test_logical_operations(kl);
    Sierpinskian sl=true;
    test_logical_operations(sl);

    Integer n(3);
    test_ring_operations(n);
    Rational q(1,n);
    test_field_operations(q);
    Real r(q);
    test_elementary_operations(r);
    ValidatedReal vr(r);
    test_elementary_operations(vr);

    Differential<Real> d=Differential<Real>::coordinate(r);
    test_algebra_operations(d,r);
    Differential<ValidatedReal> vd=Differential<ValidatedReal>::coordinate(r);
    test_algebra_operations(vd,r);
    test_algebra_operations(vd,vr);
    Algebra<Real> a(d);
    Algebra<ValidatedReal> va(vd);
    test_algebra_operations(a,r);
    test_algebra_operations(va,r);
    test_algebra_operations(va,vr);

    Variable<Real> v("x");
    Expression<Real> e=v;
    test_elementary_operations(e);
    Expression<ValidatedReal> ve=e;
    test_elementary_operations(ve);
    test_algebra_operations(e,r);
    test_algebra_operations(ve,r);
    test_algebra_operations(ve,vr);

    EuclideanSpace dom(2);
    Function<Real> f=Function<Real>::constant(dom,r);
    Function<ValidatedReal> vf=f;
    test_algebra_operations(f,r);
    test_algebra_operations(vf,r);
    test_algebra_operations(vf,vr);
    test_mixed_operations(vf,f);

}

int main() {
    wrap();
    test();

    Integer n1(3);
    Integer n2(5);
    PRINT(n1);
    Rational q1(1,n1);
    Rational q2(1,n2);
    PRINT(q1);
    Real r1(q1);
    Real r2(q2);
    PRINT(r1);
    ValidatedReal vr1=q1;
    ValidatedReal vr2=r2;
    PRINT(vr1);

    PRINT(add(r1,r2));

    Precision64 pr;
    PRINT(vr1.get(pr));

    Float64Bounds x0(pr);
    PRINT(x0);
    Float64Bounds x1(r1,pr);
    PRINT(x1);
    Float64Bounds x2(r2,pr);
    PRINT(x2);

    x0=add(x1,x2);
    x0=add(x1,r2);
    x0=add(r1,x2);

    PRINT(add(x1,x2));
    PRINT(add(x1,r2));
    PRINT(add(r1,x2));
    PRINT(mul(x1,x2));
    PRINT(mul(x1,r2));
    PRINT(mul(r1,x2));
    PRINT(neg(x2));
    PRINT(rec(x2));
    PRINT(exp(x1));
    PRINT(log(x2));

    vr1=x1;
    PRINT(vr1);
    vr2=x2;
    PRINT(vr2);

    Variable<Real> x("x");
    Variable<Real> y("y");
    PRINT(x);
    Expression<Real> e=log(mul(add(x,1),exp(y)));
    PRINT(e);
    Algebra<Real> a=static_cast<Algebra<Real>>(e);
    PRINT(a);

    e=r1;
    PRINT(e);
    a=r1;
    PRINT(a);

    e=mul(add(x,r1),exp(y));
    PRINT(e);
    Expression<ValidatedReal> ve(e);
    PRINT(ve);


    // FIXME: Should be no need to use explicit cast here!
    ve=mul(add(x,Expression<ValidatedReal>::constant(vr1)),exp(y));
    PRINT(ve);

    e=mul(x,exp(y));
    r1=q1;
    r2=q2;
    PRINT(e);
    std::map<RealVariable,Real> m{{x,r1},{y,r2}};
    Valuation<Real> v=m;
    PRINT(v);
    PRINT(e(v));
    PRINT(evaluate(e,v));
    std::cout << std::endl;
    PRINT(evaluate(x,v));
    PRINT(evaluate(add(x,y),v));
    EuclideanSpace dom(2);
    auto id0=Function<Real>::coordinate(dom,0);
    auto id1=Function<Real>::coordinate(dom,1);
    PRINT(mul(id0,exp(id1)));

    std::map<RealVariable,RealFunction> mf{{x,id0},{y,id1}};
    Valuation<RealFunction> fv(mf);
    PRINT(evaluate_algebra(add(x,y),fv));
    Function<Real> f({x,y},e);
    PRINT(f);

    Vector<RealExpression> ev={x,y};
    PRINT(ev);
    PRINT(fv);
    PRINT(f);
    e=evaluate_algebra(f,ev);
    PRINT(e);
    f=evaluate_algebra(e,fv);
    PRINT(f);

    Float64Bounds third=rec(Float64Bounds(3,3));
    std::cout << std::setprecision(17);
    PRINT(third);
    PRINT(third.error_bound());

    vr1=x1; vr2=r2;
    Vector<Real> rv={r1,r2};
    Vector<ValidatedReal> vrv={vr1,vr2};
    Vector<Float64Bounds> xv={x1,x2};

    PRINT(f(rv));
    PRINT(evaluate(f,rv));

//    PRINT(f(vrv));
    PRINT(evaluate(f,vrv));
    PRINT(f(xv));
    PRINT(evaluate(f,xv));

    ValidatedRealFunction vf(f);
    PRINT(vf(vrv));
    PRINT(evaluate(vf,vrv));
    PRINT(vf(xv));
//    PRINT(evaluate(vf,xv));

    Differential<Real> dr1=Differential<Real>::coordinate(r1);
    Differential<Real> dr2=Differential<Real>::coordinate(r2);
    Vector<Differential<Real>> vdr={dr1,dr2};
    PRINT(f(vdr));
    PRINT(evaluate_algebra(f,vdr));

    Differential<Float64Bounds> dx1=Differential<Float64Bounds>::coordinate(x1);
    Differential<Float64Bounds> dx2=Differential<Float64Bounds>::coordinate(x2);
    Vector<Differential<Float64Bounds>> vdx={dx1,dx2};

//    PRINT(f(vdx));
//    PRINT(evaluate_algebra(f,vdx));
    PRINT(vf(vdx));
    PRINT(evaluate_algebra(vf,vdx));

}
