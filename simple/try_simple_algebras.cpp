#include <type_traits>
#include "simple_algebras.h"

template<template<typename>class ALG> void test_algebra() {
    Real r; Float x;
    ALG<Real> ar; ALG<Float> ax; ax=ar;
    neg(ar); add(ar,ar); add(ar,r); add(r,ar);
    neg(ax); add(ax,ax); add(ax,x); add(x,ax);
    add(ax,ar); add(ar,ax); add(ax,r); add(r,ax);
};

int main() {
    test_algebra<FriendAlgebra>();
    test_algebra<DispatchAlgebra>();
    test_algebra<SfinaeAlgebra>();
    test_algebra<TemplateAlgebra>();
    test_algebra<ExplicitAlgebra>();
    test_algebra<ConcreteFriendAlgebra>();
    test_algebra<ConcreteDispatchAlgebra>();
    test_algebra<ConcreteSfinaeAlgebra>();
    test_algebra<ConcreteTemplateAlgebra>();
    test_algebra<ConcreteExplicitAlgebra>();

//    SfinaeAlgebra<Real> ra; mul(ra,ra);
}
