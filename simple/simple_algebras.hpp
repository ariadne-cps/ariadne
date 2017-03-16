#include <type_traits>
template<class P, class T=bool> using EnableIf = typename std::enable_if<P::value,T>::type;
template<class F, class T> using IsConvertible = std::is_convertible<F,T>;

template<class X> struct ObjectTraits { typedef X SelfType; };
template<class X> using SelfType = typename ObjectTraits<X>::SelfType;
template<class X> using NumericType = typename X::NumericType;
template<class X> using GenericType = typename ObjectTraits<X>::GenericType;



template<class A> struct Operations {
    typedef NumericType<A> X;
    static A _neg(A); static A _add(A,A); static A _add(A,X); static A _add(X,A);
//    static A _add(A, GenericType<A>);
};




template<class A, class X> struct DeclareAlgbraOperators {
    friend A add(A a1, A a2);
    friend A neg(A a);
    friend A add(A a1, X x2);
    friend A add(X x1, A a2);
//    friend A add(A a1, GenericType<A> a2);
};


template<class X> class FriendAlgebra
    : public DeclareAlgbraOperators<FriendAlgebra<X>,X>
{
  public:
    typedef X NumericType;
    FriendAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> FriendAlgebra(FriendAlgebra<R>) { }
};


template<class X> class ConcreteFriendAlgebra
    : public DeclareAlgbraOperators<ConcreteFriendAlgebra<X>,X>
{
  public:
    typedef X NumericType;
    ConcreteFriendAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> ConcreteFriendAlgebra(ConcreteFriendAlgebra<R>) { }
};



template<class A, class X> struct StaticDispatchAlgbraOperators {
    friend A add(A a1, A a2) { return Operations<A>::_add(a1,a2); }
    friend A neg(A a) { return Operations<A>::_neg(a); }
    friend A add(A a1, X x2) { return Operations<A>::_add(a1,x2); }
    friend A add(X x1, A a2) { return Operations<A>::_add(x1,a2); }
//    friend A add(A a1, GenericType<A> a2) { return Operations<A>::_add(a1,a1.create(a2)); }
};


template<class X> class DispatchAlgebra
    : public StaticDispatchAlgbraOperators<DispatchAlgebra<X>,X>
{
  public:
    typedef X NumericType;
    DispatchAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> DispatchAlgebra(DispatchAlgebra<R>) { }
};


template<class X> class ConcreteDispatchAlgebra
    : public StaticDispatchAlgbraOperators<ConcreteDispatchAlgebra<X>,X>
{
  public:
    typedef X NumericType;
    ConcreteDispatchAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> ConcreteDispatchAlgebra(ConcreteDispatchAlgebra<R>) { }
};




struct AlgebraTag { };
template<class A> using IsAlgebra = IsConvertible<A,AlgebraTag>;
template<class A> using IfAlgebra = EnableIf<IsAlgebra<A>,A>;

//template<class A, EnableIf<IsAlgebra<A>> =true> inline A neg(A a) { return Operations<A>::_neg(a); }
template<class A> inline IfAlgebra<A> neg(A a) { return Operations<A>::_neg(a); }
template<class A, EnableIf<IsAlgebra<A>> =true> inline A add(A a1, A a2) { return Operations<A>::_add(a1,a2); }
template<class A, EnableIf<IsAlgebra<A>> =true> inline A add(A a1, NumericType<A> c2) { return Operations<A>::_add(a1,c2); }
template<class A, EnableIf<IsAlgebra<A>> =true> inline A add(NumericType<A> c1, A a2) { return Operations<A>::_add(c1,a2); }
template<class A, EnableIf<IsAlgebra<A>> =true> A mul(A a1, A a2);
template<class A, EnableIf<IsAlgebra<A>> =true> inline A add(A a1, SelfType<A> a2) { return add(a1,a2); }
template<class A, EnableIf<IsAlgebra<A>> =true> inline A add(SelfType<A> a1, A a2) { return add(a1,a2); }
template<class A, EnableIf<IsAlgebra<A>> =true> inline A add(A a1, GenericType<A> a2) { return add(a1,a1.create(a2)); }

template<class X> class SfinaeAlgebra : public AlgebraTag {
  public:
    typedef X NumericType;
    SfinaeAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> SfinaeAlgebra(SfinaeAlgebra<R>) { }
};

template<class X> class ConcreteSfinaeAlgebra : public AlgebraTag {
  public:
    typedef X NumericType;
    ConcreteSfinaeAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> ConcreteSfinaeAlgebra(ConcreteSfinaeAlgebra<R>) { }
};



template<template<typename>class F> struct DeclareTemplateAlgebraOperators {
    template<class X> friend F<X> add(F<X>, F<X>);
    template<class X> friend F<X> neg(F<X>);
    template<class X> friend F<X> add(F<X>, NumericType<F<X>>);
    template<class X> friend F<X> add(NumericType<F<X>>, F<X>);
};

template<template<typename>class F> struct DefineTemplateMixedSelfOperators {
    template<class X> friend F<X> add(F<X> x1, SelfType<F<X>> x2) { return add(x1,x2); }
    template<class X> friend F<X> add(SelfType<F<X>> x1, F<X> x2) { return add(x1,x2); }
    template<class X> friend F<X> add(F<X> x1, GenericType<F<X>> x2) { return add(x1,x1.create(x2)); }
};

template<class X> class TemplateAlgebra
    : public DeclareTemplateAlgebraOperators<TemplateAlgebra>
    , public DefineTemplateMixedSelfOperators<TemplateAlgebra>
{
  public:
    typedef X NumericType;
    TemplateAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> TemplateAlgebra(TemplateAlgebra<R>) { }
};

template<class X> class ConcreteTemplateAlgebra
    : public DeclareTemplateAlgebraOperators<ConcreteTemplateAlgebra>
    , public DefineTemplateMixedSelfOperators<ConcreteTemplateAlgebra>
{
  public:
    typedef X NumericType;
    ConcreteTemplateAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> ConcreteTemplateAlgebra(ConcreteTemplateAlgebra<R>) { }
};



template<class X> class ExplicitAlgebra {
  public:
    typedef X NumericType;
    ExplicitAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> ExplicitAlgebra(ExplicitAlgebra<R>) { }
};
template<class X> ExplicitAlgebra<X> neg(ExplicitAlgebra<X> a);
template<class X> ExplicitAlgebra<X> add(ExplicitAlgebra<X> a1, ExplicitAlgebra<X> a2);
template<class X> ExplicitAlgebra<X> add(ExplicitAlgebra<X> a1, NumericType<ExplicitAlgebra<X>> c2);
template<class X> ExplicitAlgebra<X> add(NumericType<ExplicitAlgebra<X>> c1, ExplicitAlgebra<X> a2);

template<class X> inline ExplicitAlgebra<X> add(ExplicitAlgebra<X> a1, SelfType<ExplicitAlgebra<X>> a2) { return add(a1,a2); }
template<class X> inline ExplicitAlgebra<X> add(SelfType<ExplicitAlgebra<X>> a1, ExplicitAlgebra<X> a2) { return add(a1,a2); }
template<class X> inline ExplicitAlgebra<X> add(ExplicitAlgebra<X> a1, GenericType<ExplicitAlgebra<X>> a2) { return add(a1,a1.create(a2)); }

template<class X> class ConcreteExplicitAlgebra {
  public:
    typedef X NumericType;
    ConcreteExplicitAlgebra() { }
    template<class R, EnableIf<IsConvertible<R,X>> =true> ConcreteExplicitAlgebra(ConcreteExplicitAlgebra<R>) { }
};
template<class X> ConcreteExplicitAlgebra<X> neg(ConcreteExplicitAlgebra<X> a);
template<class X> ConcreteExplicitAlgebra<X> add(ConcreteExplicitAlgebra<X> a1, ConcreteExplicitAlgebra<X> a2);
template<class X> ConcreteExplicitAlgebra<X> add(ConcreteExplicitAlgebra<X> a1, NumericType<ConcreteExplicitAlgebra<X>> c2);
template<class X> ConcreteExplicitAlgebra<X> add(NumericType<ConcreteExplicitAlgebra<X>> c1, ConcreteExplicitAlgebra<X> a2);

template<class X> inline ConcreteExplicitAlgebra<X> add(ConcreteExplicitAlgebra<X> a1, SelfType<ConcreteExplicitAlgebra<X>> a2) { return add(a1,a2); }
template<class X> inline ConcreteExplicitAlgebra<X> add(SelfType<ConcreteExplicitAlgebra<X>> a1, ConcreteExplicitAlgebra<X> a2) { return add(a1,a2); }
template<class X> inline ConcreteExplicitAlgebra<X> add(SelfType<ConcreteExplicitAlgebra<X>> a1, GenericType<ConcreteExplicitAlgebra<X>> a2) { return add(a1,a1.create(a2)); }



class Real { public: Real(){} };
class Float { public: Float(){} Float(Real){} };
