#include "simple_algebras.h"
using SizeType = std::size_t;

template<class X> SizeType instantiate_template_algebra_operations();




FriendAlgebra<Real> add(FriendAlgebra<Real> a1, FriendAlgebra<Real> a2) { return FriendAlgebra<Real>(); }
FriendAlgebra<Real> add(FriendAlgebra<Real> a1, NumericType<FriendAlgebra<Real>> c2) { return FriendAlgebra<Real>(); }
FriendAlgebra<Real> add(NumericType<FriendAlgebra<Real>> c1, FriendAlgebra<Real> a2) { return FriendAlgebra<Real>(); }
FriendAlgebra<Real> neg(FriendAlgebra<Real> a) { return FriendAlgebra<Real>(); }

FriendAlgebra<Float> add(FriendAlgebra<Float> a1, FriendAlgebra<Float> a2) { return FriendAlgebra<Float>(); }
FriendAlgebra<Float> add(FriendAlgebra<Float> a1, NumericType<FriendAlgebra<Float>> c2) { return FriendAlgebra<Float>(); }
FriendAlgebra<Float> add(NumericType<FriendAlgebra<Float>> c1, FriendAlgebra<Float> a2) { return FriendAlgebra<Float>(); }
FriendAlgebra<Float> neg(FriendAlgebra<Float> a) { return FriendAlgebra<Float>(); }

ConcreteFriendAlgebra<Real> add(ConcreteFriendAlgebra<Real> a1, ConcreteFriendAlgebra<Real> a2) { return ConcreteFriendAlgebra<Real>(); }
ConcreteFriendAlgebra<Real> add(ConcreteFriendAlgebra<Real> a1, NumericType<ConcreteFriendAlgebra<Real>> c2) { return ConcreteFriendAlgebra<Real>(); }
ConcreteFriendAlgebra<Real> add(NumericType<ConcreteFriendAlgebra<Real>> c1, ConcreteFriendAlgebra<Real> a2) { return ConcreteFriendAlgebra<Real>(); }
ConcreteFriendAlgebra<Real> neg(ConcreteFriendAlgebra<Real> a) { return ConcreteFriendAlgebra<Real>(); }

ConcreteFriendAlgebra<Float> add(ConcreteFriendAlgebra<Float> a1, ConcreteFriendAlgebra<Float> a2) { return ConcreteFriendAlgebra<Float>(); }
ConcreteFriendAlgebra<Float> add(ConcreteFriendAlgebra<Float> a1, NumericType<ConcreteFriendAlgebra<Float>> c2) { return ConcreteFriendAlgebra<Float>(); }
ConcreteFriendAlgebra<Float> add(NumericType<ConcreteFriendAlgebra<Float>> c1, ConcreteFriendAlgebra<Float> a2) { return ConcreteFriendAlgebra<Float>(); }
ConcreteFriendAlgebra<Float> neg(ConcreteFriendAlgebra<Float> a) { return ConcreteFriendAlgebra<Float>(); }




template<class X> struct Operations<DispatchAlgebra<X>> {
    typedef DispatchAlgebra<X> A;
    static A _neg(A a) { return A(); }
    static A _add(A a1, A a2) { return A(); }
    static A _add(A a1, X x2) { return A(); }
    static A _add(X x1, A a2) { return A(); }
};

template class Operations<DispatchAlgebra<Real>>;
template class Operations<DispatchAlgebra<Float>>;

template<class X> struct Operations<ConcreteDispatchAlgebra<X>> {
    typedef ConcreteDispatchAlgebra<X> A;
    static A _neg(A a) { return A(); }
    static A _add(A a1, A a2) { return A(); }
    static A _add(A a1, X x2) { return A(); }
    static A _add(X x1, A a2) { return A(); }
};

template class Operations<ConcreteDispatchAlgebra<Real>>;
template class Operations<ConcreteDispatchAlgebra<Float>>;



template<class X> struct Operations<SfinaeAlgebra<X>> {
    typedef SfinaeAlgebra<X> A;
    static A _neg(A a) { return A(); }
    static A _add(A a1, A a2) { return A(); }
    static A _add(A a1, X x2) { return A(); }
    static A _add(X x1, A a2) { return A(); }
};

template class Operations<SfinaeAlgebra<Real>>;
template class Operations<SfinaeAlgebra<Float>>;


template<class X> struct Operations<ConcreteSfinaeAlgebra<X>> {
    typedef ConcreteSfinaeAlgebra<X> A;
    static A _neg(A a) { return A(); }
    static A _add(A a1, A a2) { return A(); }
    static A _add(A a1, X x2) { return A(); }
    static A _add(X x1, A a2) { return A(); }
};

template class Operations<ConcreteSfinaeAlgebra<Real>>;
template class Operations<ConcreteSfinaeAlgebra<Float>>;

template<class X> SfinaeAlgebra<X> mul(SfinaeAlgebra<X>, SfinaeAlgebra<X>);

template<class X> TemplateAlgebra<X> add(TemplateAlgebra<X> a1, TemplateAlgebra<X> a2) { return TemplateAlgebra<X>(); }
template<class X> TemplateAlgebra<X> add(TemplateAlgebra<X> a1, NumericType<TemplateAlgebra<X>> c2) { return TemplateAlgebra<X>(); }
template<class X> TemplateAlgebra<X> add(NumericType<TemplateAlgebra<X>> c1, TemplateAlgebra<X> a2) { return TemplateAlgebra<X>(); }
template<class X> TemplateAlgebra<X> neg(TemplateAlgebra<X> a) { return TemplateAlgebra<X>(); }
template SizeType instantiate_template_algebra_operations<TemplateAlgebra<Real>>();
template SizeType instantiate_template_algebra_operations<TemplateAlgebra<Float>>();

template<class X> ConcreteTemplateAlgebra<X> add(ConcreteTemplateAlgebra<X> a1, ConcreteTemplateAlgebra<X> a2) { return ConcreteTemplateAlgebra<X>(); }
template<class X> ConcreteTemplateAlgebra<X> add(ConcreteTemplateAlgebra<X> a1, NumericType<ConcreteTemplateAlgebra<X>> c2) { return ConcreteTemplateAlgebra<X>(); }
template<class X> ConcreteTemplateAlgebra<X> add(NumericType<ConcreteTemplateAlgebra<X>> c1, ConcreteTemplateAlgebra<X> a2) { return ConcreteTemplateAlgebra<X>(); }
template<class X> ConcreteTemplateAlgebra<X> neg(ConcreteTemplateAlgebra<X> a) { return ConcreteTemplateAlgebra<X>(); }
template SizeType instantiate_template_algebra_operations<ConcreteTemplateAlgebra<Real>>();
template SizeType instantiate_template_algebra_operations<ConcreteTemplateAlgebra<Float>>();




template<class X> ExplicitAlgebra<X> add(ExplicitAlgebra<X> a1, ExplicitAlgebra<X> a2) { return ExplicitAlgebra<X>(); }
template<class X> ExplicitAlgebra<X> add(ExplicitAlgebra<X> a1, NumericType<ExplicitAlgebra<X>> c2) { return ExplicitAlgebra<X>(); }
template<class X> ExplicitAlgebra<X> add(NumericType<ExplicitAlgebra<X>> c1, ExplicitAlgebra<X> a2) { return ExplicitAlgebra<X>(); }
template<class X> ExplicitAlgebra<X> neg(ExplicitAlgebra<X> a) { return ExplicitAlgebra<X>(); }
template SizeType instantiate_template_algebra_operations<ExplicitAlgebra<Real>>();
template SizeType instantiate_template_algebra_operations<ExplicitAlgebra<Float>>();


template<class X> ConcreteExplicitAlgebra<X> add(ConcreteExplicitAlgebra<X> a1, ConcreteExplicitAlgebra<X> a2) { return ConcreteExplicitAlgebra<X>(); }
template<class X> ConcreteExplicitAlgebra<X> add(ConcreteExplicitAlgebra<X> a1, NumericType<ConcreteExplicitAlgebra<X>> c2) { return ConcreteExplicitAlgebra<X>(); }
template<class X> ConcreteExplicitAlgebra<X> add(NumericType<ConcreteExplicitAlgebra<X>> c1, ConcreteExplicitAlgebra<X> a2) { return ConcreteExplicitAlgebra<X>(); }
template<class X> ConcreteExplicitAlgebra<X> neg(ConcreteExplicitAlgebra<X> a) { return ConcreteExplicitAlgebra<X>(); }
template SizeType instantiate_template_algebra_operations<ConcreteExplicitAlgebra<Real>>();
template SizeType instantiate_template_algebra_operations<ConcreteExplicitAlgebra<Float>>();




template<class A> SizeType instantiate_template_algebra_operations() {
    typedef NumericType<A> X;
    auto add_cnst_ptr=(A(*)(A,X)) &add;
    auto cnst_add_ptr=(A(*)(X,A)) &add;
    auto add_ptr=(A(*)(A,A)) &add;
    auto neg_ptr=(A(*)(A)) &neg;

    return (SizeType)add_ptr + (SizeType)neg_ptr + (SizeType)add_cnst_ptr + (SizeType)cnst_add_ptr;
};
