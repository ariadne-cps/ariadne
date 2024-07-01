#include <cstddef>
#include <cassert>
#include <cmath>
//#include <type_traits>
#include <concepts>
#include <iostream>
#include <sstream>
#include <functional>
#include <memory>
#include <string>
#include <variant>
#include <vector>
#include <map>

template<class F, class T> concept ConvertibleTo = std::convertible_to<F,T>;
template<class F, class... AS> concept Invocable = std::invocable<F,AS...>;
template<class F, class R, class... AS> concept InvocableReturning = std::invocable<F,AS...> and std::convertible_to<typename std::invoke_result_t<F,AS...>,R>;

template<class T> using RemoveCVType = std::remove_cv_t<T>;
template<class T> using RemoveReferenceType = std::remove_reference_t<T>;
template<class T> using RemoveType = std::remove_cvref_t<T>;

using std::declval;
using Int = int;
using SizeType = std::size_t;
using String = std::string;
using OutputStream = std::ostream;
template<class T> using SharedPointer = std::shared_ptr<T>;
template<class SIG> using Function = std::function<SIG>;
template<class T> using List = std::vector<T>;

struct Neg {
    template<class A> decltype(auto) operator() (A a) const { return neg(a); }
    static constexpr const char name[]="neg";
};
struct Rec {
    template<class A> decltype(auto) operator() (A a) const { return rec(a); }
    static constexpr const char name[]="rec";
};
struct Exp {
    template<class A> decltype(auto) operator() (A a) const { return exp(a); }
    static constexpr const char name[]="exp";
};
struct Add {
    template<class A1, class A2> decltype(auto) operator() (A1 a1, A2 a2) const { return add(a1,a2); }
    static constexpr const char* name="add";
};
struct Mul {
    template<class A1, class A2> decltype(auto) operator() (A1 a1, A2 a2) const { return mul(a1,a2); }
    static constexpr const char* name="mul";
};
struct Pow {
    template<class A, class N> decltype(auto) operator() (A a, N n) const { return pow(a,n); }
    static constexpr const char name[]="pow";
};
struct Get {
    template<class V, class I> decltype(auto) operator() (V const& v, I i) const { return v[i]; }
    static constexpr const char name[]="get";
};


template<class OP, class SIG> concept AnOperator = ConvertibleTo<OP,Function<SIG>> and requires(OP op) {
    op.name;
};

template<class OP> concept HasName = requires(OP op) { { op.name }  -> ConvertibleTo<String>; };

template<HasName OP> OutputStream& operator<<(OutputStream& os, OP op) { return os << op.name; }

template<class... OPS> struct OperatorVariant : public std::variant<OPS...> {
    template<class... AS> decltype(auto) operator() (AS... as) const { return std::visit([&](auto op){return op(as...);}); }
    friend OutputStream& operator<<(OutputStream& os, OperatorVariant<OPS...> ops) {
        return os << std::visit([&](auto op){return op.name;}); }
};

template<class R, class... AS> struct AllowedOperatorsTrait;
template<class R, class... AS> using AllowedOperatorsType = typename AllowedOperatorsTrait<R(AS...)>::Type;

class RealInterface;

class Real {
    SharedPointer<RealInterface> _ptr;
  public:
        Real();
    Real(int n);
    explicit Real(double x) = delete;
    explicit Real(SharedPointer<RealInterface> p) : _ptr(p) { }
    friend Real neg(Real);
    friend Real rec(Real);
    friend Real exp(Real);
    friend Real add(Real,Real);
    friend Real mul(Real,Real);
    friend Real pow(Real, Int);
    friend OutputStream& operator<<(OutputStream& os, Real const& r);
  public:
    static Real make(double x);
    double get_d() const;
};

template<> struct AllowedOperatorsTrait<Real(Real)> { typedef OperatorVariant<Neg,Rec,Exp> Type; };
template<> struct AllowedOperatorsTrait<Real(Real,Real)> { typedef OperatorVariant<Add,Mul> Type; };

class RealInterface {
  public:
    virtual ~RealInterface() = default;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};
struct RealConstant : public RealInterface {
    double _x;
    RealConstant(double x) : _x(x) { }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_x; }

};
Real::Real() : Real(0) { }
Real::Real(int n) : _ptr(std::make_shared<RealConstant>(n)) { }
OutputStream& operator<<(OutputStream& os, Real const& r) { return r._ptr->_write(os); }

double Real::get_d() const { auto p=dynamic_cast<RealConstant const*>(this->_ptr.operator->()); assert(p); return p->_x; }
Real Real::make(double x) { return Real(std::make_shared<RealConstant>(x)); }
Real neg(Real r) { return Real::make(-r.get_d()); }
Real add(Real r1, Real r2) { return Real::make(r1.get_d()+r2.get_d()); }
Real mul(Real r1, Real r2) { return Real::make(r1.get_d()*r2.get_d()); }
Real rec(Real r) { return Real::make(1/r.get_d()); }
Real exp(Real r) { return Real::make(std::exp(r.get_d())); }
Real pow(Real r, Int n) { return Real::make(std::pow(r.get_d(),n)); }

static_assert(AnOperator<Add,Real(Real,Real)>);
static_assert(not AnOperator<Add,Real(Real)>);

template<class T> using Scalar = T;

template<class T> class Vector {
    List<T> _v;
  public:
    template<class G> requires InvocableReturning<G,T,SizeType> Vector(SizeType n, G const& g) {
        this->_v.reserve(n); for (SizeType i=0; i!=n; ++i) { this->_v.push_back(g(i)); } }
    SizeType size() const { return this->_v.size(); }
    T const& operator[] (SizeType i) const { return this->_v[i]; }
    T const& get(SizeType i) const { return this->_v[i]; }
    friend OutputStream& operator<<(OutputStream& os, Vector<T> const& v) {
        os << "{"; if (v.size()!=0) { os << v[0u]; for (SizeType i=1u; i!=v.size(); ++i) { os << "," << v[i]; } } os << "}"; return os; }
};
using RealVector = Vector<Real>;
Real get(RealVector const& v, SizeType i) { return v[i]; }
RealVector neg(RealVector);
RealVector add(RealVector v1, RealVector v2) {
    assert(v1.size()==v2.size()); return RealVector(v1.size(),[&v1,&v2](SizeType i){return add(v1[i],v2[i]);}); }
RealVector mul(RealVector v1, Real s2) {
    return RealVector(v1.size(),[&v1,&s2](SizeType i){return mul(v1[i],s2);}); }
RealVector mul(Real s1, RealVector v2) {
    return RealVector(v2.size(),[&s1,&v2](SizeType i){return mul(s1,v2[i]);}); }
RealVector lmul(RealVector v1, Real s2) {
    return RealVector(v1.size(),[&v1,&s2](SizeType i){return mul(v1[i],s2);}); }
RealVector rmul(Real s1, RealVector v2) {
    return RealVector(v2.size(),[&s1,&v2](SizeType i){return mul(s1,v2[i]);}); }

template<class T> class Constant {
    String _name; T _value;
  public:
    Constant(String name, T const& value) : _name(name), _value(value) { }
    //Constant(T const& value) : Constant(to_str(value),value) { }
    Constant(T const& value) : Constant(""  ,value) { }
    operator T const& () const { return this->_value; }
    friend OutputStream& operator<<(OutputStream& os, Constant<T> const& c) {
        if (!c._name.empty()) { os << c._name << "="; } os << c._value; return os; }
  private:
    String to_str(T const& value) { std::stringstream ss; ss << value; return ss.str(); }
};


/*
OutputStream& operator<<(OutputStream& os, Function<Real(Real)> op) {
    //if (op == Function<Real(Real)>(static_cast<Real(Real)>(&operator-))) { return os << 'sub'; }
    return os << "unknown";
}

OutputStream& operator<<(OutputStream& os, Function<Real(Real,Real)> op) {
    if (op == Function<Real(Real,Real)>(&add)) { return os << "add"; }
    return os << "unknown";
}
*/

template<class T> class Variable {
    String _name;
  public:
    explicit Variable(String n) : _name(n) { }
    String const& name() const { return this->_name; }
    friend bool operator<(Variable<T> const& v1, Variable<T> const& v2) { return v1.name()<v2.name(); }
    friend OutputStream& operator<<(OutputStream& os, Variable<T> const& v) { return os << v.name(); }
};

template<class T, class X=T> class Valuation {
    X _nul;
    std::map<Variable<T>,X> _values;
    std::map<Variable<Vector<T>>,Vector<X>> _vector_values;
    template<class K, class V> static V& at(std::map<K,V>& m, K const& k) {
        return m[k]; }
    template<class K, class V> static V const& get(std::map<K,V> const& m, K const& k) {
        auto p=m.find(k); assert(p!=m.end()); return p->second; }
  public:
    Valuation() : _values(), _nul() { }
    X& operator[](Variable<T> var) { return at(this->_values,var); }
    Vector<X> const& operator[](Variable<Vector<T>> var) const { return get(this->_values,var); }
    Vector<X>& operator[](Variable<Vector<T>> var) { return at(this->_values,var); }
    X const& operator[](Variable<T> var) const { return get(this->_values,var); }
    X create_constant(T const& c) const { X r(this->_nul); r=c; return r; }
    friend OutputStream& operator<<(OutputStream& os, Valuation<T,X> const& v) {
        os << "{"; for (auto t : v._values) { os << t.first << ":" << t.second << ", "; } os << "}"; return os; }
};

template<class OP, class...> struct Symbolic;

template<class OP, class A> struct Symbolic<OP,A> {
    OP _op; A _a;
    Symbolic(OP op, A a) : _op(op), _a(a) { }
    friend OutputStream& operator<<(OutputStream& os, Symbolic<OP,A> const& sym) { return os << sym._op << "(" << sym._a << ")"; }
};
template<class OP, class A1, class A2> struct Symbolic<OP,A1,A2> {
    OP _op; A1 _a1; A2 _a2;
    Symbolic(OP op, A1 a1, A2 a2) : _op(op), _a1(a1), _a2(a2) { }
    friend OutputStream& operator<<(OutputStream& os, Symbolic<OP,A1,A2> const& sym) { return os << sym._op << "(" << sym._a1 << "," << sym._a2 << ")"; }
};

class UntypedExpressionNode {
  public:
    virtual ~UntypedExpressionNode() = default;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};
template<class T> class ExpressionNode;

template<class SIG> struct ArgumentTraits;
template<class T, class A1> struct ArgumentTraits<T(A1)> {
    typedef A1 FirstArgumentType; };
template<class T, class A1, class A2> struct ArgumentTraits<T(A1,A2)> {
    typedef A1 FirstArgumentType; typedef A2 SecondArgumentType; };
template<class SIG> using FirstArgumentType = typename ArgumentTraits<SIG>::FirstArgumentType;
template<class SIG> using SecondArgumentType = typename ArgumentTraits<SIG>::SecondArgumentType;

template<class T> struct SelfTrait { typedef T Type; };
template<class T> using SelfType = typename SelfTrait<T>::Type;

template<class T> class Expression;

template<class T> struct IsExpression : std::false_type { };
template<class T> struct IsExpression<Expression<T>> : std::true_type { };
template<class T> struct IsExpression<Variable<T>> : std::true_type { };
template<class T> struct IsExpression<Constant<T>> : std::true_type { };
template<class T> concept AnExpression = IsExpression<T>::value;

template<class T> struct TermTrait;
template<class T> struct TermTrait<Constant<T>> { typedef T Type; };
template<class T> struct TermTrait<Variable<T>> { typedef T Type; };
template<class T> struct TermTrait<Expression<T>> { typedef T Type; };
template<class T> using TermType = typename TermTrait<T>::Type;

template<class OP, class E> concept AUnaryExpressionOperator
    = AnExpression<E> and Invocable<OP,TermType<E>>;

template<class OP, class E1, class E2> concept ABinaryExpressionOperator
    = AnExpression<E1> and AnExpression<E2> and Invocable<OP,TermType<E1>,TermType<E2>>;

template<class OP, class E1, class A2> concept AGradedExpressionOperator
    = AnExpression<E1> and (not AnExpression<A2>) and Invocable<OP,TermType<E1>,A2>;


class ExpressionOperators {
    template<class T1, class T2> friend decltype(auto) add(Expression<T1> const& e1, Expression<T2> const& e2) {
        return Expression<decltype(add(declval<T1>(),declval<T2>()))>(Add(),e1,e2); }
    template<class T1, class T2> friend decltype(auto) mul(Expression<T1> const& e1, Expression<T2> const& e2) {
        return Expression<decltype(mul(declval<T1>(),declval<T2>()))>(Mul(),e1,e2); }
    template<class T1, class T2> friend decltype(auto) mul(Expression<T1> const& e1, T2 const& c2) {
        return mul(e1,Expression<T2>(Constant<T2>(c2))); }
};

template<class T> class Expression : public ExpressionOperators {
    SharedPointer<ExpressionNode<T>> _node;
  public:
    explicit Expression(SharedPointer<ExpressionNode<T>> node) : _node(node) { }
    const ExpressionNode<T>* node_ptr() const { return this->_node.operator->(); }

    Expression(Constant<T> c);
    Expression(Variable<T> v);

//    template<class OP, class A> requires AnOperator<OP,T(A)> Expression(OP op, Expression<A> const& e);
//    template<class OP, class A1, class A2> requires AnOperator<OP,T(A1,A2)> Expression(OP op, Expression<A1> const& e1, Expression<A2> const& e2);

    template<class A> Expression(T op(A), Expression<RemoveType<A>> const& e);
    template<class A1, class A2> Expression(T op(A1,A2), Expression<RemoveType<A1>> const& e1, Expression<RemoveType<A2>> const& e2);
    template<class A1, class A2> Expression(T op(A1,A2), Expression<RemoveType<A1>> const& e1, SelfType<A2> const& c2);

    template<class OP, class E> requires AUnaryExpressionOperator<OP,E> Expression(OP op, E const& e);
    template<class OP, class E1, class E2> requires ABinaryExpressionOperator<OP,E1,E2> Expression(OP op, E1 const& e1, E2 const& e2);
    template<class OP, class E1, class A2> requires AGradedExpressionOperator<OP,E1,A2> Expression(OP op, E1 const& e1, A2 const& a2);

//    template<class A, class OP, class E> requires AnOperator<OP,T(A)> and AnExpression<E,A> Expression(OP op, E const& e);
//    template<class A1, class A2, class OP, class E1, class E2> requires AnOperator<OP,T(A1,A2)> and AnExpression<E1,A1> and AnExpression<E2,A2> Expression(OP op, E1 const& e1, E2 const& e2);

//    template<class A> Expression(Function<T(A)> op, Expression<SelfType<A>> const& e);
//    template<class A1, class A2> Expression(Function<T(A1,A2)> op, Expression<SelfType<A1>> const& e1, Expression<SelfType<A2>> const& e2);

    friend OutputStream& operator<<(OutputStream& os, Expression<T> const& e) {
        return e._node->_write(os); }


};

template<class SIG> OutputStream& operator<<(OutputStream& os, Function<SIG> const& f) { return os << "op"; }


template<class T> class ExpressionNode : public UntypedExpressionNode { };

template<class T> class ConstantExpressionNode : public ExpressionNode<T>, public Constant<T> {
  public:
    ConstantExpressionNode(T const& value) : Constant<T>(value) { }
    ConstantExpressionNode(String name, T const& value) : Constant<T>(name,value) { }
    ConstantExpressionNode(Constant<T> const& cnst) : Constant<T>(cnst) { }
    virtual OutputStream& _write(OutputStream& os) const override { return os << static_cast<Constant<T>const&>(*this); }
};

template<class T> class VariableExpressionNode : public ExpressionNode<T>, public Variable<T> {
  public:
    VariableExpressionNode(String name) : Variable<T>(name) { }
    VariableExpressionNode(Variable<T> var) : Variable<T>(var) { }
    virtual OutputStream& _write(OutputStream& os) const override { return os << static_cast<Variable<T>const&>(*this); }
};

template<class T, class OP, class... AS> class SymbolicExpressionNode : public ExpressionNode<T>, public Symbolic<OP,AS...> {
  public:
    SymbolicExpressionNode(Symbolic<OP,AS...> s) : Symbolic<OP,AS...>(s) { }
    SymbolicExpressionNode(OP op, AS... as) : Symbolic<OP,AS...>(op,as...) { }
    virtual OutputStream& _write(OutputStream& os) const override { return os << static_cast<Symbolic<OP,AS...>const&>(*this); }
};

template<class T> Expression<T>::Expression(Constant<T> c) : _node(std::make_shared<ConstantExpressionNode<T>>(c)) { }
template<class T> Expression<T>::Expression(Variable<T> v) : _node(std::make_shared<VariableExpressionNode<T>>(v)) { }

template<class T, class A> auto
make_expression_node(Function<T(A)> op, Expression<RemoveType<A>> e) -> SharedPointer<ExpressionNode<T>> {
    using O = Function<T(A)>; using E=Expression<RemoveType<A>>;
    return SharedPointer<ExpressionNode<T>>(new SymbolicExpressionNode<T,O,E>(op,e)); }
template<class T, class A1, class A2> auto
make_expression_node(Function<T(A1,A2)> op, Expression<RemoveType<A1>> e1, Expression<RemoveType<A2>> e2) -> SharedPointer<ExpressionNode<T>> {
    using O = Function<T(A1,A2)>; using E1=Expression<RemoveType<A1>>; using E2=Expression<RemoveType<A2>>;
    return SharedPointer<ExpressionNode<T>>(new SymbolicExpressionNode<T,O,E1,E2>(op,e1,e2)); }
template<class T, class A1, class A2> auto
make_expression_node(Function<T(A1,A2)> op, Expression<RemoveType<A1>> e1, A2 a2) -> SharedPointer<ExpressionNode<T>> {
    using O = Function<T(A1,A2)>; using E1=Expression<RemoveType<A1>>;
    return SharedPointer<ExpressionNode<T>>(new SymbolicExpressionNode<T,O,E1,A2>(op,e1,a2)); }

template<class T, class OP, class E> requires AUnaryExpressionOperator<OP,E> auto
make_expression_node(OP op, E e) -> SharedPointer<ExpressionNode<T>> {
    using A = TermType<E>;
    return SharedPointer<ExpressionNode<T>>(new SymbolicExpressionNode<T,OP,Expression<A>>(op,e)); }
template<class T, class OP, class E1, class E2> requires ABinaryExpressionOperator<OP,E1,E2> auto
make_expression_node(OP op, E1 e1, E2 e2) -> SharedPointer<ExpressionNode<T>> {
    using A1=TermType<E1>; using A2=TermType<E2>;
    return SharedPointer<ExpressionNode<T>>(new SymbolicExpressionNode<T,OP,Expression<A1>,Expression<A2>>(op,e1,e2)); }
template<class T, class OP, class E1, class A2> requires AGradedExpressionOperator<OP,E1,A2> auto
make_expression_node(OP op, E1 e1, A2 a2) -> SharedPointer<ExpressionNode<T>> {
    using A1=TermType<E1>;
    return SharedPointer<ExpressionNode<T>>(new SymbolicExpressionNode<T,OP,Expression<A1>,A2>(op,e1,a2)); }


template<class T> template<class OP, class E> requires AUnaryExpressionOperator<OP,E>
Expression<T>::Expression(OP op, E const& e)
    :  _node(make_expression_node<T>(op,e)) { }
template<class T> template<class OP, class E1, class E2> requires ABinaryExpressionOperator<OP,E1,E2>
Expression<T>::Expression(OP op, E1 const& e1, E2 const& e2)
    :  _node(make_expression_node<T>(op,e1,e2)) { }
template<class T> template<class OP, class E1, class A2> requires AGradedExpressionOperator<OP,E1,A2>
Expression<T>::Expression(OP op, E1 const& e1, A2 const& a2)
    :  _node(make_expression_node<T>(op,e1,a2)) { }

//template<class T> template<class A1, class A2>
//Expression<T>::Expression(Function<T(A1,A2)> op, Expression<SelfType<A1>> const& e1, Expression<SelfType<A2>> const& e2)
//    :  _node(make_expression_node<T>(op,e1,e2)) { }

template<class T> template<class A>
Expression<T>::Expression(T op(A), Expression<RemoveType<A>> const& e)
    :  _node(make_expression_node<T>(Function<T(A)>(op),e)) { }
template<class T> template<class A1, class A2>
Expression<T>::Expression(T op(A1,A2), Expression<RemoveType<A1>> const& e1, Expression<RemoveType<A2>> const& e2)
    :  _node(make_expression_node<T>(Function<T(A1,A2)>(op),e1,e2)) { }
template<class T> template<class A1, class A2>
Expression<T>::Expression(T op(A1,A2), Expression<RemoveType<A1>> const& e1, SelfType<A2> const& a2)
    :  _node(make_expression_node<T>(Function<T(A1,A2)>(op),e1,a2)) { }

template<class T> Expression<T> get(Expression<Vector<T>> const& ve, SizeType i) { return Expression<T>(Get(),ve,i); }

template<class T, class X> X evaluate(Expression<T> const& e, Valuation<T,X> const& v);

template<class T, class X> struct ExpressionEvaluator {
    static X evaluate(T const& c, Valuation<T,X> const& v) { return v.create_constant(c); }
    static X evaluate(Constant<T> const& c, Valuation<T,X> const& v) { return v.create_constant(static_cast<T const&>(c)); }
    static X evaluate(Variable<T> const& x, Valuation<T,X> const& v) { return v[x]; }
    template<class OP, class A> static decltype(auto) evaluate(Symbolic<OP,A> const& s, Valuation<T,X> const& v) { return s._op(evaluate(s._a,v)); }
    template<class OP, class A1, class A2> static decltype(auto) evaluate(Symbolic<OP,A1,A2> const& s, Valuation<T,X> const& v) { return s._op(evaluate(s._a1,v),evaluate(s._a2,v)); }

    static X evaluate(Symbolic<Get,Expression<Vector<T>>,SizeType> const& s, Valuation<T,X> const& v) { return evaluate(s._a1,v)[s._a2]; }

    static Vector<X> evaluate(Constant<Vector<T>> const& c, Valuation<T,X> const& v) { return c; }

    static X evaluate(Expression<T> const& e, Valuation<T,X> const& v) {
        ExpressionNode<T> const* p=e.node_ptr();
        if (auto cp = dynamic_cast<ConstantExpressionNode<T>const*>(p)) {
            return evaluate(*cp,v); }
        else if (auto vp = dynamic_cast<VariableExpressionNode<T>const*>(p)) {
            return evaluate(*vp,v); }
        else if (auto sp = dynamic_cast<SymbolicExpressionNode<T,Rec,Expression<T>>const*>(p)) {
            return evaluate(*sp,v); }
        else if (auto sp = dynamic_cast<SymbolicExpressionNode<T,Add,Expression<T>,Expression<T>>const*>(p)) {
            return evaluate(*sp,v); }
        else if (auto sp = dynamic_cast<SymbolicExpressionNode<T,Mul,Expression<T>,Expression<T>>const*>(p)) {
            return evaluate(*sp,v); }
        else if (auto gp = dynamic_cast<SymbolicExpressionNode<T,Get,Expression<Vector<T>>,SizeType>const*>(p)) {
            return evaluate(*gp,v); }
        else {
            std::cerr<<"Unevaluable expression "<<e<<"\n"; assert(false);
        }
    }

    static Vector<X> evaluate(Expression<Vector<T>> const& e, Valuation<T,X> const& v) {
        //std::cerr<<"Vector<X> Evaluator::evaluate(Expression<Vector<T>> const& e, Valuation<T,X> const& v)" << std::endl;
        ExpressionNode<Vector<T>> const* p=e.node_ptr();
        //std::cerr<<"p="<<(void*)p<<"; "; p->_write(std::cerr); std::cerr << std::endl;
        if (auto cp = dynamic_cast<ConstantExpressionNode<Vector<T>>const*>(p)) {
            return evaluate(*cp,v); }
        else if (auto vp = dynamic_cast<VariableExpressionNode<Vector<T>>const*>(p)) {
            return evaluate(*vp,v); }
        else if (auto sp = dynamic_cast<SymbolicExpressionNode<Vector<T>,Add,Expression<Vector<T>>,Expression<Vector<T>>>const*>(p)) {
            return evaluate(*sp,v); }
        else if (auto sp = dynamic_cast<SymbolicExpressionNode<Vector<T>,Mul,Expression<T>,Expression<Vector<T>>>const*>(p)) {
            return evaluate(*sp,v); }
        else if (auto sp = dynamic_cast<SymbolicExpressionNode<Vector<T>,Mul,Expression<Vector<T>>,Expression<T>>const*>(p)) {
            return evaluate(*sp,v); }
        else {
            std::cerr<<"Unevaluable expression "<<e<<"\n"; assert(false);
        }
    }

};

template<class T, class X> X evaluate(Expression<T> const& e, Valuation<T,X> const& v) {
    return ExpressionEvaluator<T,X>::evaluate(e,v); }
template<class T, class X> Vector<X> evaluate(Expression<Vector<T>> const& e, Valuation<T,X> const& v) {
    return ExpressionEvaluator<T,X>::evaluate(e,v); }

#define PRINT(e) { std::cout << #e << ": " << (e) << std::endl; }

int main() {
    Constant<Real> five(Real(5));
    //Constant<Real> five("five",Real(5));
    Variable<Real> x("x");
    auto e=Expression<Real>(five);
/*
    e=Expression<Real>(&rec,five);
    e=Expression<Real>(&rec,x);
    e=Expression<Real>(&rec,e);
    e=Expression<Real>(&add,five,x);
    e=Expression<Real>(&mul,five,x);
*/
    e=Expression<Real>(Rec(),five);
    e=Expression<Real>(Rec(),x);
    e=Expression<Real>(Rec(),e);
    e=Expression<Real>(Add(),five,x);
    e=Expression<Real>(Mul(),five,x);

    Valuation<Real> val; val[x]=Real(2);
    evaluate(e,val);
//    PRINT(e);
//    PRINT(val);
//    PRINT(evaluate(e,val));

//    e=Expression<Real>(Mul(),five,x);
//    std::cout << e << "\n";

    const int primes[] = {2,3,5,7,11};
    Constant<Vector<Real>> cv(RealVector(3,[&](SizeType i){return primes[i];}));
    auto ev=Expression<Vector<Real>>(cv);
    PRINT(e);
    PRINT(ev);
    PRINT(evaluate(ev,val));
    PRINT(evaluate(add(ev,ev),val));
    PRINT(evaluate(mul(ev,e),val));
//    PRINT(evaluate(mul(ev,3),val));
/*
    ev=Expression<RealVector>(&add,ev,ev);
    ev=Expression<RealVector>(&lmul,ev,e);
    ev=Expression<RealVector>(&rmul,e,ev);
    e=Expression<Real>(&get,ev,(SizeType)2u);
*/
    SizeType i=1u;
    PRINT(evaluate(get(mul(ev,e),i),val));
    return 0;
    ev=Expression<RealVector>(Add(),ev,ev);
    ev=Expression<RealVector>(Mul(),e,ev);
    ev=Expression<RealVector>(Mul(),ev,e);
    e=Expression<Real>(Get(),ev,i);
    PRINT(ev);
    PRINT(evaluate(ev,val));
    PRINT(e);
    PRINT(val);
    PRINT(evaluate(e,val));


}
