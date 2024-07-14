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
#include <tuple>
#include <vector>
#include <set>
#include <map>

template<class T1, class T2> concept Same = std::same_as<T1,T2>;
template<class F, class T> concept ConvertibleTo = std::convertible_to<F,T>;
template<class T, class... FS> concept ConstructibleFrom = std::constructible_from<T,FS...>;
template<class F, class... AS> concept Invocable = std::invocable<F,AS...>;
template<class F, class R, class... AS> concept InvocableReturning = std::invocable<F,AS...> and std::convertible_to<typename std::invoke_result_t<F,AS...>,R>;
template<class D, class B> concept DerivedFrom = std::derived_from<D,B>;

template<class T> using RemoveCVType = std::remove_cv_t<T>;
template<class T> using RemoveReferenceType = std::remove_reference_t<T>;
template<class T> using RemoveType = std::remove_cvref_t<T>;

using std::declval;
using Void = void;
using Bool = bool;
using Int = int;
using SizeType = std::size_t;
using String = std::string;
using OutputStream = std::ostream;
template<class T> using SharedPointer = std::shared_ptr<T>;
template<class SIG> using Function = std::function<SIG>;
template<class T> using List = std::vector<T>;
template<class T> using Set = std::set<T>;
template<class K, class V> using Map = std::map<K,V>;
template<class T> using InitializerList = std::initializer_list<T>;
template<class... TS> using Variant = std::variant<TS...>;
template<class... TS> using Tuple = std::tuple<TS...>;

template<class T> String to_str(T const& t) { std::stringstream ss; ss<<t; return ss.str(); }

template<class OS, class TUP, std::size_t N> inline void write_tuple(OS& os, TUP const& tup, std::integral_constant<std::size_t,N>) {
    if constexpr (N!=1) { write_tuple(os,tup,std::integral_constant<std::size_t,N-1>()); os << ','; } os << std::get<N-1>(tup); }
template<class OS, class... TS> OS& operator<<(OS& os, Tuple<TS...> const& tup) {
    typename std::tuple_size<std::tuple<TS...>>::type sz; os << "("; write_tuple(os,tup,sz); os << ")"; return os;}

template<class T> OutputStream& operator<<(OutputStream& os, List<T> const& l) {
    bool fst=true; os << "["; for (auto t : l) { if (fst) { fst=false; } else { os << ","; } os << t; } os << "]"; return os; }

template<class K, class V> OutputStream& operator<<(OutputStream& os, Map<K,V> const& m) {
    bool fst=true; os << "{"; for (auto t : m) { if (fst) { fst=false; } else { os << ","; } os << t.first << ":" << t.second; } os << "}"; return os; }

struct Pos {
    template<class A> decltype(auto) operator() (A a) const { return pos(a); }
    static constexpr const char name[]="pos"; static constexpr const char symbol = '+';
};
struct Neg {
    template<class A> decltype(auto) operator() (A a) const { return neg(a); }
    static constexpr const char name[]="neg"; static constexpr const char symbol = '-';
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
    static constexpr const char* name="add"; static constexpr const char symbol = '+';
};
struct Sub {
    template<class A1, class A2> decltype(auto) operator() (A1 a1, A2 a2) const { return sub(a1,a2); }
    static constexpr const char* name="sub"; static constexpr const char symbol = '-';
};
struct Mul {
    template<class A1, class A2> decltype(auto) operator() (A1 a1, A2 a2) const { return mul(a1,a2); }
    static constexpr const char* name="mul"; static constexpr const char symbol = '*';
};
struct Div {
    template<class A1, class A2> decltype(auto) operator() (A1 a1, A2 a2) const { return div(a1,a2); }
    static constexpr const char* name="div"; static constexpr const char symbol = '/';
};
struct Pow {
    template<class A, class N> decltype(auto) operator() (A a, N n) const { return pow(a,n); }
    static constexpr const char name[]="pow";
};
struct Get {
    template<class V, class I> decltype(auto) operator() (V const& v, I i) const { return v[i]; }
    static constexpr const char name[]="get";
};

template<class T1, class T2> using AddType = decltype(add(declval<T1>(),declval<T2>()));
template<class T1, class T2> using SubType = decltype(sub(declval<T1>(),declval<T2>()));

template<class OP, class SIG> concept AnOperator = ConvertibleTo<OP,Function<SIG>> and requires(OP op) { op.name; };

template<class... OPS> struct OperatorVariant : public Variant<OPS...> {
    template<class OP> requires ConvertibleTo<OP,Variant<OPS...>> OperatorVariant(OP op) : Variant<OPS...>(op) { }
//    template<class... ROPS> requires ConvertibleTo<Variant<ROPS...>,Variant<OPS...>> OperatorVariant(OperatorVariant<ROPS...> op) : Variant<OPS...>(op) { }
    template<class... AS> decltype(auto) operator() (AS... as) const { return std::visit([&](auto op){return op(as...);},*this); }
    friend OutputStream& operator<<(OutputStream& os, OperatorVariant<OPS...> ops) {
        std::visit([&](auto op){os<<op.name;},ops); return os; }
};


class RealInterface;

class Real {
    SharedPointer<RealInterface> _ptr;
  public:
        Real();
    Real(int n);
    explicit Real(double x) = delete;
    explicit Real(SharedPointer<RealInterface> p) : _ptr(p) { }
    friend Real pos(Real);
    friend Real neg(Real);
    friend Real rec(Real);
    friend Real exp(Real);
    friend Real add(Real,Real);
    friend Real sub(Real,Real);
    friend Real mul(Real,Real);
    friend Real div(Real,Real);
    friend Real pow(Real, Int);
    friend OutputStream& operator<<(OutputStream& os, Real const& r);
  public:
    static Real make(double x);
    double get_d() const;
};

Real operator ""_dec(long double x) { return Real::make(x); }

template<class R, class... AS> struct AllowedOperatorsTrait;
template<class R, class... AS> using AllowedOperatorsType = typename AllowedOperatorsTrait<R(AS...)>::Type;
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
double Real::get_d() const { auto p=dynamic_cast<RealConstant const*>(this->_ptr.operator->()); assert(p); return p->_x; }
Real Real::make(double x) { return Real(std::make_shared<RealConstant>(x)); }
OutputStream& operator<<(OutputStream& os, Real const& r) { return r._ptr->_write(os); }
Real pos(Real r) { return Real::make(+r.get_d()); }
Real neg(Real r) { return Real::make(-r.get_d()); }
Real add(Real r1, Real r2) { return Real::make(r1.get_d()+r2.get_d()); }
Real sub(Real r1, Real r2) { return Real::make(r1.get_d()-r2.get_d()); }
Real mul(Real r1, Real r2) { return Real::make(r1.get_d()*r2.get_d()); }
Real div(Real r1, Real r2) { return Real::make(r1.get_d()/r2.get_d()); }
Real rec(Real r) { return Real::make(1/r.get_d()); }
Real exp(Real r) { return Real::make(std::exp(r.get_d())); }
Real pow(Real r, Int n) { return Real::make(std::pow(r.get_d(),n)); }

static_assert(AnOperator<Add,Real(Real,Real)>);
static_assert(not AnOperator<Add,Real(Real)>);

/*--------------- Vector ----------------------------------------*/

template<class T> using Scalar = T;

template<class T> class Vector {
    List<T> _v;
  public:
    Vector(Vector<T>const& v)
        : Vector(v.size(), [&v](SizeType i){return static_cast<T>(v[i]);}) { }
    Vector(InitializerList<T> v) : _v(v) { }
    template<class TT> requires (not Same<TT,T>) and ConvertibleTo<TT,T> Vector(Vector<TT>const& v)
        : Vector(v.size(), [&v](SizeType i){return static_cast<T>(v[i]);}) { }
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
template<class T> Vector<T> pos(Vector<T> v) {
    return Vector<T>(v.size(),[&v](SizeType i){return pos(v[i]);}); }
template<class T> Vector<T> neg(Vector<T> v) {
    return Vector<T>(v.size(),[&v](SizeType i){return neg(v[i]);}); }
template<class T1, class T2> Vector<AddType<T1,T2>> add(Vector<T1> v1, Vector<T2> v2) {
    assert(v1.size()==v2.size()); return Vector<AddType<T1,T2>>(v1.size(),[&v1,&v2](SizeType i){return add(v1[i],v2[i]);}); }
template<class T1, class T2> Vector<SubType<T1,T2>> sub(Vector<T1> v1, Vector<T2> v2) {
    assert(v1.size()==v2.size()); return Vector<SubType<T1,T2>>(v1.size(),[&v1,&v2](SizeType i){return sub(v1[i],v2[i]);}); }
template<class T> Vector<T> mul(Vector<T> v1, T s2) {
    return Vector<T>(v1.size(),[&v1,&s2](SizeType i){return mul(v1[i],s2);}); }
template<class T> Vector<T> mul(T s1, Vector<T> v2) {
    return Vector<T>(v2.size(),[&s1,&v2](SizeType i){return mul(s1,v2[i]);}); }
template<class T> Vector<T> div(Vector<T> v1, T s2) {
    return Vector<T>(v1.size(),[&v1,&s2](SizeType i){return div(v1[i],s2);}); }
template<class T> Vector<T> lmul(Vector<T> v1, T s2) {
    return Vector<T>(v1.size(),[&v1,&s2](SizeType i){return mul(v1[i],s2);}); }
template<class T> Vector<T> rmul(T s1, Vector<T> v2) {
    return Vector<T>(v2.size(),[&s1,&v2](SizeType i){return mul(s1,v2[i]);}); }

template<class X> class Expression;

template<class X> concept AScalar = Same<X,Real> or Same<X,Expression<Real>>;

template<class X> requires AScalar<X> inline Vector<X> vectorise(X s) { return Vector<X>{s}; }
template<class X> requires AScalar<X> inline Vector<X> const& vectorise(Vector<X> const& v) { return v; }

template<class X> requires AScalar<X> Vector<X> join(X s) { return vectorise(s); }
template<class X> requires AScalar<X> Vector<X> join(Vector<X> v) { return v; }
template<class X> Vector<X> join(Vector<X> v1, Vector<X> v2) {
    auto n1=v1.size(); auto n2=v2.size(); return Vector<X>( n1+n2, [n1,&v1,&v2](SizeType i) { return i<n1 ? v1[i] : v2[i-n1]; } ); }

template<class T0, class T1, class... TS> decltype(auto) join(T0 t0, T1 t1, TS... ts) { return join(join(vectorise(t0),vectorise(t1)),ts...); }

using RealScalar = Scalar<Real>;
using RealVector = Vector<Real>;


/**************** SYMBOLIC *****************************************/



template<class T> class Constant;
template<class T> class Variable;
using RealVariable = Variable<Real>;
using RealScalarVariable = Variable<Scalar<Real>>;
using RealVectorVariable = Variable<Vector<Real>>;
template<class T> class Expression;
using RealExpression = Expression<Real>;
using RealVectorExpression = Expression<Vector<Real>>;

/*
class TemplatedExpressionOperators {
    template<class T1, class T2> friend decltype(auto) add(Expression<T1> const& e1, Expression<T2> const& e2) {
        return Expression<decltype(add(declval<T1>(),declval<T2>()))>(Add(),e1,e2); }
    template<class T1, class T2> friend decltype(auto) mul(Expression<T1> const& e1, Expression<T2> const& e2) {
        return Expression<decltype(mul(declval<T1>(),declval<T2>()))>(Mul(),e1,e2); }
    template<class T1, class T2> friend decltype(auto) mul(Expression<T1> const& e1, T2 const& c2) {
        return mul(e1,Expression<T2>(Constant<T2>(c2))); }
    template<class T1, class T2> friend decltype(auto) mul(T1 const& c1, Expression<T2> const& e2) {
        return mul(Expression<T1>(Constant<T1>(c1)),e2); }
};
*/

class TemplatedOperators {
/*
    template<class T> friend decltype(auto) operator+(T t) { return t; }
    template<class T> friend decltype(auto) operator-(T t) { return neg(t); }
    template<class T1, class T2> friend decltype(auto) operator+(T1 t1, T2 t2) { return add(t1,t2); }
    template<class T1, class T2> friend decltype(auto) operator-(T1 t1, T2 t2) { return sub(t1,t2); }
    template<class T1, class T2> friend decltype(auto) operator*(T1 t1, T2 t2) { return mul(t1,t2); }
    template<class T1, class T2> friend decltype(auto) operator/(T1 t1, T2 t2) { return div(t1,t2); }
*/
};

template<class T> class ExpressionOperations;

template<> class ExpressionOperations<Real> : public TemplatedOperators {
    friend RealExpression neg(RealExpression);
    friend RealExpression rec(RealExpression);
    friend RealExpression exp(RealExpression);
    friend RealExpression add(RealExpression, RealExpression);
    friend RealExpression sub(RealExpression, RealExpression);
    friend RealExpression mul(RealExpression, RealExpression);
    friend RealExpression div(RealExpression, RealExpression);
    friend RealExpression operator+(RealExpression, RealExpression);
    friend RealExpression operator-(RealExpression, RealExpression);
    friend RealExpression operator*(RealExpression, RealExpression);
    friend RealExpression operator/(RealExpression, RealExpression);
};

template<> class ExpressionOperations<Vector<Real>> : public TemplatedOperators {
/*
    friend RealExpression neg(RealExpression);
    friend RealExpression rec(RealExpression);
    friend RealExpression exp(RealExpression);
    friend RealExpression add(RealExpression, RealExpression);
    friend RealExpression sub(RealExpression, RealExpression);
    friend RealExpression mul(RealExpression, RealExpression);
    friend RealExpression div(RealExpression, RealExpression);
*/
    friend RealVectorExpression operator+(RealVectorExpression, RealVectorExpression);
    friend RealVectorExpression operator-(RealVectorExpression, RealVectorExpression);
    friend RealVectorExpression operator*(RealExpression, RealVectorExpression);
    friend RealVectorExpression operator*(RealVectorExpression, RealExpression);
    friend RealVectorExpression operator/(RealVectorExpression, RealExpression);
};


template<class T> class Constant {
    String _name; T _value;
  public:
    Constant(String name, T const& value) : _name(name), _value(value) { }
    //Constant(T const& value) : Constant(to_str(value),value) { }
    Constant(T const& value) : Constant(""  ,value) { }
    operator T const& () const { return this->_value; }
    String const& name() const { return this->_name; }
    T const& value() const { return this->_value; }
    friend OutputStream& operator<<(OutputStream& os, Constant<T> const& c) {
        if (!c._name.empty()) { os << c._name << "="; } os << c._value; return os; }
  private:
    String to_str(T const& value) { std::stringstream ss; ss << value; return ss.str(); }
};

template<class T> class Constant<Vector<T>> {
    String _name; Vector<T> _value;
  public:
    Constant(String name, Vector<T> const& value) : _name(name), _value(value) { }
    //Constant(T const& value) : Constant(to_str(value),value) { }
    Constant(Vector<T> const& value) : Constant(""  ,value) { }
    operator Vector<T> const& () const { return this->_value; }
    String const& name() const { return this->_name; }
    Vector<T> const& value() const { return this->_value; }
    Constant<T> operator[](SizeType i) const {
        if (this->_name=="") { return Constant<T>("",this->_value[i]); }
        else { return Constant<T>(this->_name+'['+to_str(i)+']',this->_value[i]); } }
    friend OutputStream& operator<<(OutputStream& os, Constant<Vector<T>> const& c) {
        if (!c._name.empty()) { os << c._name << "="; } os << c._value; return os; }
  private:
    static String to_str(T const& value) { std::stringstream ss; ss << value; return ss.str(); }
    static String to_str(SizeType i) { std::stringstream ss; ss << i; return ss.str(); }
};

template<class T> class Expression;

template<class T> class Variable
    : public ExpressionOperations<T>
{
    String _name;
  public:
    explicit Variable(String n) : _name(n) { }
    String const& name() const { return this->_name; }
    friend bool operator<(Variable<T> const& v1, Variable<T> const& v2) { return v1.name()<v2.name(); }
    friend OutputStream& operator<<(OutputStream& os, Variable<T> const& v) { return os << v.name(); }
};
template<class T> class Variable<Vector<T>>
    : public ExpressionOperations<Vector<T>>
{
    String _name; SizeType _size;
  public:
    explicit Variable(String n, SizeType s) : _name(n), _size(s) { }
    String const& name() const { return this->_name; }
    SizeType size() const { return this->_size; }
    Expression<T> operator[](SizeType i) const;
    friend bool operator<(Variable<Vector<T>> const& v1, Variable<Vector<T>> const& v2) { return v1.name()<v2.name(); }
    //friend OutputStream& operator<<(OutputStream& os, Variable<Vector<T>> const& v) { return os << v.name() << "⟨" << v.size() << "⟩"; }
    friend OutputStream& operator<<(OutputStream& os, Variable<Vector<T>> const& v) { return os << v.name(); }
};

template<class K, class V> struct MapAt {
    std::map<K,V>& m; K k;
    Void operator=(V const& v) { auto p=m.find(k); if (p==m.end()) { m.insert(std::pair<K,V>(k,v)); } else { p->second=v; } }
    operator V const& () const { auto p=m.find(k); assert(p!=m.end()); return p->second; }
};

template<class K, class V> struct MapGet {
    std::map<K,V> const& m; K k;
    operator V const& () const { auto p=m.find(k); assert(p!=m.end()); return p->second; }
};

template<class T, class X=T> class Valuation {
    std::map<Variable<T>,X> _values;
    std::map<Variable<Vector<T>>,Vector<X>> _vector_values;
    template<class K, class V> static MapAt<K,V> at(std::map<K,V>& m, K const& k) {
        return MapAt<K,V>{m,k}; }
    template<class K, class V> static V const& get(std::map<K,V> const& m, K const& k) {
        auto p=m.find(k); if (p==m.end()) { std::cerr << "\nm=" << m << "\n  k=" << k << "\n"; } assert(p!=m.end()); return p->second; }
  public:
    Valuation() : _values() { }
    MapAt<Variable<T>,X> operator[](Variable<T> var) { return at(this->_values,var); }
    Vector<X> const& operator[](Variable<Vector<T>> var) const { return get(this->_vector_values,var); }
    MapAt<Variable<Vector<T>>,Vector<X>> operator[](Variable<Vector<T>> var) { return at(this->_vector_values,var); }
    X const& operator[](Variable<T> var) const { return get(this->_values,var); }
    X create_constant(T const& c) const {
        if constexpr (ConstructibleFrom<X,T>) { return X(c); } else { X r=this->_values.begin()->second; r=c; return r; } }
    friend OutputStream& operator<<(OutputStream& os, Valuation<T,X> const& v) {
        os << "{";
        for (auto t : v._values) { os << t.first << ":" << t.second << ", "; }
        for (auto t : v._vector_values) { os << t.first << ":" << t.second << ", "; }
        os << "}"; return os; }
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

template<class OP, class... AS> Symbolic(OP, AS...) -> Symbolic<OP,AS...>;

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



template<class T> class Expression;
using RealExpression = Expression<Real>;
using RealVectorExpression = Expression<Vector<Real>>;
using RealExpressionVector = Vector<Expression<Real>>;

template<class T> struct ExpressionVariantTrait;
template<class T> using ExpressionVariantType = typename ExpressionVariantTrait<T>::Type;
template<> struct ExpressionVariantTrait<Real> {
    typedef Variant<   Constant<Real>
                    ,  Variable<Real>
                    ,  Symbolic<OperatorVariant<Pos,Neg,Rec,Exp>,Expression<Real>>
                    ,  Symbolic<OperatorVariant<Add,Sub,Mul,Div>,Expression<Real>,Expression<Real>>
                    ,  Symbolic<OperatorVariant<Get>,Expression<Vector<Real>>,SizeType>
                     > Type;
};
template<> struct ExpressionVariantTrait<Vector<Real>> {
    typedef Variant< Constant<Vector<Real>>,
                     Variable<Vector<Real>>,
                     Vector<Expression<Real>>,
                     Symbolic<OperatorVariant<Pos,Neg>,Expression<Vector<Real>>>,
                     Symbolic<OperatorVariant<Add,Sub>,Expression<Vector<Real>>,Expression<Vector<Real>>>,
                     Symbolic<OperatorVariant<Mul>,Expression<Real>,Expression<Vector<Real>>>,
                     Symbolic<OperatorVariant<Mul,Div>,Expression<Vector<Real>>,Expression<Real>>
                    > Type;
};

template<class T> class ExpressionNode : public virtual UntypedExpressionNode, public ExpressionVariantType<T> {
  public:
    ExpressionNode(ExpressionVariantType<T> opv) : ExpressionVariantType<T>(opv) { }
    virtual OutputStream& _write(OutputStream& os) const;
};

template<class T> class Expression;

/*template<class T> class Expression : public ExpressionOperators {
    SharedPointer<ExpressionNode<T>> _node;
  public:
    explicit Expression(SharedPointer<ExpressionNode<T>> node) : _node(node) { }
    const ExpressionNode<T>* node_ptr() const { return this->_node.operator->(); }

    Expression(Constant<T> c);
    Expression(Variable<T> v);


    friend OutputStream& operator<<(OutputStream& os, Expression<T> const& e) {
        return e._node->_write(os); }
};
*/

template<class T> using ScalarOrVectorVariable = Variant<Variable<Scalar<T>>,Variable<Vector<T>>>;

template<> class Expression<Real> : public ExpressionOperations<Real> {
    typedef Real T;
    SharedPointer<ExpressionNode<T>> _node;
  public:
    explicit Expression(SharedPointer<ExpressionNode<T>> node) : _node(node) { }
    const ExpressionNode<T>* node_ptr() const { return this->_node.operator->(); }
    const ExpressionNode<T>& node_ref() const { return this->_node.operator*(); }

    Expression();
    Expression(T const& c);
    Expression<T>& operator=(T const& c);

    Expression(Constant<T> c);
    Expression(Variable<T> v);
    Expression(OperatorVariant<Pos,Neg,Rec,Exp>, Expression<T> const&);
    Expression(OperatorVariant<Add,Sub,Mul,Div>, Expression<T> const&, Expression<T> const&);

    Expression(OperatorVariant<Get>, Expression<Vector<T>>, SizeType);
    friend OutputStream& operator<<(OutputStream& os, Expression<T> const& e);

    friend Expression<T> pos(Expression<T> e) { return Expression<T>(Pos(),e); }
    friend Expression<T> neg(Expression<T> e) { return Expression<T>(Neg(),e); }
    friend Expression<T> rec(Expression<T> e) { return Expression<T>(Rec(),e); }
    friend Expression<T> exp(Expression<T> e) { return Expression<T>(Exp(),e); }
    friend Expression<T> add(Expression<T> e1, Expression<T> e2) { return Expression<T>(Add(),e1,e2); }
    friend Expression<T> sub(Expression<T> e1, Expression<T> e2) { return Expression<T>(Sub(),e1,e2); }
    friend Expression<T> mul(Expression<T> e1, Expression<T> e2) { return Expression<T>(Mul(),e1,e2); }
    friend Expression<T> div(Expression<T> e1, Expression<T> e2) { return Expression<T>(Div(),e1,e2); }
    friend Expression<T> operator+(Expression<T> e1, Expression<T> e2) { return add(e1,e2); }
    friend Expression<T> operator*(Expression<T> e1, Expression<T> e2) { return mul(e1,e2); }

     Set<ScalarOrVectorVariable<Real>> variables() const;
};

/*
Expression<Real> neg(Expression<Real> e) { return Expression<Real>(Neg(),e); }
Expression<Real> rec(Expression<Real> e) { return Expression<Real>(Rec(),e); }
Expression<Real> exp(Expression<Real> e) { return Expression<Real>(Exp(),e); }
Expression<Real> add(Expression<Real> e1, Expression<Real> e2) { return Expression<Real>(Add(),e1,e2); }
Expression<Real> sub(Expression<Real> e1, Expression<Real> e2) { return Expression<Real>(Sub(),e1,e2); }
Expression<Real> mul(Expression<Real> e1, Expression<Real> e2) { return Expression<Real>(Mul(),e1,e2); }
Expression<Real> div(Expression<Real> e1, Expression<Real> e2) { return Expression<Real>(Div(),e1,e2); }
*/


/*
template<class T, class OP, class E> requires AUnaryExpressionOperator<OP,E> auto
make_expression_node(OP op, E e) -> SharedPointer<ExpressionNode<T>> {
    using A = TermType<E>;
    return SharedPointer<ExpressionNode<T>>(new ExpressionNode<T>(Symbolic<OperatorVariant<Neg,Rec,Exp>,Expression<A>>(op,e))); }
template<class T, class OP, class E1, class E2> requires ABinaryExpressionOperator<OP,E1,E2> auto
make_expression_node(OP op, E1 e1, E2 e2) -> SharedPointer<ExpressionNode<T>> {
    using A1=TermType<E1>; using A2=TermType<E2>;
    return SharedPointer<ExpressionNode<T>>(new ExpressionNode<T>(Symbolic<OperatorVariant<Add,Sub,Mul,Div>,Expression<A1>,Expression<A2>>(op,e1,e2))); }
*/

Expression<Real>::Expression() : Expression(Constant<T>(0)) { }
Expression<Real>::Expression(T const& c) : Expression(Constant<T>(c)) { }
auto Expression<Real>::operator=(T const& c) -> Expression<T>& { return *this=Expression(Constant<T>(c)); }

Expression<Real>::Expression(Constant<T> c) : _node(std::make_shared<ExpressionNode<T>>(c)) { }

Expression<Real>::Expression(Variable<T> v) : _node(std::make_shared<ExpressionNode<T>>(v)) { }

Expression<Real>::Expression(OperatorVariant<Pos,Neg,Rec,Exp> op, Expression<T> const& e)
    :  _node(std::make_shared<ExpressionNode<T>>(Symbolic<OperatorVariant<Pos,Neg,Rec,Exp>,Expression<T>>(op,e))) { }

Expression<Real>::Expression(OperatorVariant<Add,Sub,Mul,Div> op, Expression<T> const& e1, Expression<T> const& e2)
    :  _node(std::make_shared<ExpressionNode<T>>(Symbolic<OperatorVariant<Add,Sub,Mul,Div>,Expression<T>,Expression<T>>(op,e1,e2))) { }



template<> class Expression<Vector<Real>> : public ExpressionOperations<Vector<Real>> {
    typedef Vector<Real> T; typedef Real S;
    SharedPointer<ExpressionNode<T>> _node;
  public:
    explicit Expression(SharedPointer<ExpressionNode<T>> node) : _node(node) { }
    const ExpressionNode<T>* node_ptr() const { return this->_node.operator->(); }
    const ExpressionNode<T>& node_ref() const { return this->_node.operator*(); }


    Expression(T const& c);
    Expression<T>& operator=(T const& c);

    Expression(Constant<T> c);
    Expression(Variable<T> v);
    Expression(Vector<Expression<S>> ev);
    Expression(OperatorVariant<Pos,Neg>, Expression<T>);
    Expression(OperatorVariant<Add,Sub>, Expression<T> ,Expression<T>);
    Expression(OperatorVariant<Mul>, Expression<S> ,Expression<T>);
    Expression(OperatorVariant<Mul,Div>, Expression<T> ,Expression<S>);

    SizeType size() const;
    Expression<Real> operator[](SizeType i) const;

    friend Expression<T> neg(Expression<T> e) { return Expression<T>(Neg(),e); }
    friend Expression<T> add(Expression<T> e1, Expression<T> e2) { return Expression<T>(Add(),e1,e2); }
    friend Expression<T> sub(Expression<T> e1, Expression<T> e2) { return Expression<T>(Sub(),e1,e2); }
    friend Expression<T> mul(Expression<S> e1, Expression<T> e2) { return Expression<T>(Mul(),e1,e2); }
    friend Expression<T> mul(Expression<T> e1, Expression<S> e2) { return Expression<T>(Mul(),e1,e2); }
    friend Expression<T> div(Expression<T> e1, Expression<S> e2) { return Expression<T>(Div(),e1,e2); }
    friend Expression<T> operator+(Expression<T> e1, Expression<T> e2) { return add(e1,e2); }
    friend Expression<T> operator-(Expression<T> e1, Expression<T> e2) { return sub(e1,e2); }
    friend Expression<T> operator*(Expression<S> e1, Expression<T> e2) { return mul(e1,e2); }
    friend Expression<T> operator*(Expression<T> e1, Expression<S> e2) { return mul(e1,e2); }
    friend Expression<T> operator/(Expression<T> e1, Expression<S> e2) { return div(e1,e2); }

    Set<ScalarOrVectorVariable<Real>> variables() const;
    friend OutputStream& operator<<(OutputStream& os, Expression<T> const& e);
};

Expression<Vector<Real>>::Expression(T const& c) : Expression(Constant<T>(c)) { }
auto Expression<Vector<Real>>::operator=(T const& c) -> Expression<T>& { return *this=Expression<T>(c); }

Expression<Vector<Real>>::Expression(Constant<T> c)
    : _node(std::make_shared<ExpressionNode<T>>(c)) { }
Expression<Vector<Real>>::Expression(Variable<T> v)
    : _node(std::make_shared<ExpressionNode<T>>(v)) { }
Expression<Vector<Real>>::Expression(Vector<Expression<S>> ev)
    : _node(std::make_shared<ExpressionNode<T>>(ev)) { }
Expression<Vector<Real>>::Expression(OperatorVariant<Pos,Neg> op, Expression<T> e)
    :  _node(std::make_shared<ExpressionNode<T>>(Symbolic<OperatorVariant<Pos,Neg>,Expression<T>>(op,e))) { }
Expression<Vector<Real>>::Expression(OperatorVariant<Add,Sub> op, Expression<T> e1, Expression<T> e2)
    :  _node(std::make_shared<ExpressionNode<T>>(Symbolic<OperatorVariant<Add,Sub>,Expression<T>,Expression<T>>(op,e1,e2))) { }
Expression<Vector<Real>>::Expression(OperatorVariant<Mul> op, Expression<S> e1, Expression<T> e2)
    :  _node(std::make_shared<ExpressionNode<T>>(Symbolic<OperatorVariant<Mul>,Expression<S>,Expression<T>>(op,e1,e2))) { }
Expression<Vector<Real>>::Expression(OperatorVariant<Mul,Div> op, Expression<T> e1, Expression<S> e2)
    :  _node(std::make_shared<ExpressionNode<T>>(Symbolic<OperatorVariant<Mul,Div>,Expression<T>,Expression<S>>(op,e1,e2))) { }
Expression<Real>::Expression(OperatorVariant<Get> op, Expression<Vector<T>> ve, SizeType i)
    :  _node(std::make_shared<ExpressionNode<T>>(Symbolic<OperatorVariant<Get>,Expression<Vector<T>>,SizeType>(op,ve,i))) { }


template<class T> Expression<T> get(Expression<Vector<T>> const& ve, SizeType i) {
    return Expression<T>(Get(),ve,i); }

template<> auto Variable<Vector<Real>>::operator[] (SizeType i) const -> Expression<Real> {
//    return Variable<Real>(this->name()+'['+to_str(i)+']'); }
    return Expression<Vector<Real>>(*this)[i]; }


Vector<Expression<Real>> vectorise(Expression<Real> e) { return Vector<Expression<Real>>{e}; }
Vector<Expression<Real>> vectorise(Expression<Vector<Real>> ev) { return Vector<Expression<Real>>(ev.size(),[&ev](SizeType i){return get(ev,i);}); }



template<class T, class X> X evaluate(Expression<T> const& e, Valuation<T,X> const& v);

template<class S> struct ExpressionSize {
    using V=Vector<S>;
    static SizeType _size(Constant<V> const& c) { return c.value().size(); }
    static SizeType _size(Variable<V> const& v) { return v.size(); }
    static SizeType _size(Vector<Expression<S>> const& v) { return v.size(); }
    template<class OP> static SizeType _size(Symbolic<OP,Expression<V>> const& s) { return size(s._a); }
    template<class OP> static SizeType _size(Symbolic<OP,Expression<V>,Expression<V>> const& s) { auto n1=size(s._a1); auto n2=size(s._a2); assert(n1==n2); return n1; }
    template<class OP> static SizeType _size(Symbolic<OP,Expression<V>,Expression<S>> const& s) { return size(s._a1); }
    template<class OP> static SizeType _size(Symbolic<OP,Expression<S>,Expression<V>> const& s) { return size(s._a2); }

    static SizeType size(Expression<V> const& ve) {
        return std::visit([](auto en){return _size(en);},ve.node_ref()); }
};

auto Expression<Vector<Real>>::size() const -> SizeType {
    return ExpressionSize<Real>::size(*this); }


template<class T, class X> struct ExpressionEvaluator {
    static X _evaluate(T const& c, Valuation<T,X> const& v) { return v.create_constant(c); }
    static X _evaluate(Constant<T> const& c, Valuation<T,X> const& v) { return v.create_constant(static_cast<T const&>(c)); }
    static X _evaluate(Variable<T> const& x, Valuation<T,X> const& v) { return v[x]; }
    template<class OP, class A> static decltype(auto) _evaluate(Symbolic<OP,A> const& s, Valuation<T,X> const& v) { return s._op(evaluate(s._a,v)); }
    template<class OP, class A1, class A2> static decltype(auto) _evaluate(Symbolic<OP,A1,A2> const& s, Valuation<T,X> const& v) { return s._op(evaluate(s._a1,v),evaluate(s._a2,v)); }

    template<class S> static X _evaluate(Symbolic<OperatorVariant<Get>,Expression<Vector<S>>,SizeType> const& s, Valuation<T,X> const& v) { return evaluate(s._a1,v)[s._a2]; }

    static Vector<X> _evaluate(Constant<Vector<T>> const& c, Valuation<T,X> const& v) {
        return Vector<X>(c.value().size(),[&c](SizeType i){return static_cast<X>(c.value()[i]);}); }
    static Vector<X> _evaluate(Variable<Vector<T>> const& x, Valuation<T,X> const& v) { return v[x]; }
    static Vector<X> _evaluate(Vector<Expression<T>> const& ev, Valuation<T,X> const& v) {
        return Vector<X>(ev.size(),[&ev,&v](SizeType i){return evaluate(ev[i],v);}); }

    static X evaluate(Expression<T> const& e, Valuation<T,X> const& v) {
        return std::visit([&v](auto sym){return _evaluate(sym,v);},e.node_ref()); }
    static Vector<X> evaluate(Expression<Vector<T>> const& e, Valuation<T,X> const& v) {
        return std::visit([&v](auto sym){return _evaluate(sym,v);},e.node_ref()); }
};

template<class T, class X> X evaluate(Expression<T> const& e, Valuation<T,X> const& v) {
    return ExpressionEvaluator<T,X>::evaluate(e,v); }
template<class T, class X> Vector<X> evaluate(Expression<Vector<T>> const& e, Valuation<T,X> const& v) {
    return ExpressionEvaluator<T,X>::evaluate(e,v); }
template<class T, class X> Vector<X> evaluate(Vector<Expression<T>> const& e, Valuation<T,X> const& v) {
    return ExpressionEvaluator<T,X>::evaluate(e,v); }


struct ExpressionSimplifier {
    using T=Real; using ST=Scalar<T>; using VT=Vector<T>; using VE=Expression<Vector<Real>>; using SE=Expression<Scalar<Real>>;
    static SE _simplify(ST const& c) { return SE(c); }
    static SE _simplify(Constant<T> const& c) { return c; }
    static SE _simplify(Variable<T> const& x) { return x; }
    template<class OP> static SE _simplify(Symbolic<OP,SE> const& s) { return SE(s._op,s._a); }
    template<class OP> static SE _simplify(Symbolic<OP,SE,SE> const& s) { return SE(s._op,s._a1,s._a2); }

    static VE _simplify_vector(VT const& c) { return c; }
    static VE _simplify_vector(Constant<VT> const& c) { return c; }
    static VE _simplify_vector(Variable<VT> const& x) { return x; }
    static VE _simplify_vector(Vector<SE> const& v) { return v; }
    static VE _simplify_vector(Symbolic<OperatorVariant<Pos,Neg>,VE> const& s) { return VE(s._op,s._a); }
    static VE _simplify_vector(Symbolic<OperatorVariant<Add,Sub>,VE,VE> const& s) { return VE(s._op,s._a1,s._a2); }
    static VE _simplify_vector(Symbolic<OperatorVariant<Mul,Div>,VE,SE> const& s) { return VE(s._op,s._a1,s._a2); }
    static VE _simplify_vector(Symbolic<OperatorVariant<Mul>,SE,VE> const& s) { return VE(s._op,s._a1,s._a2); }

    static SE _simplify_get(Vector<T> const& c, SizeType i) {
        return c[i]; }
    static SE _simplify_get(Constant<Vector<T>> const& c, SizeType i) {
        return c[i]; }
    static SE _simplify_get(Variable<Vector<T>> const& v, SizeType i) {
        return get(VE(v),i); }
    static SE _simplify_get(Vector<SE> const& v, SizeType i) {
        return v[i]; }

    static SE _simplify_get(Symbolic<OperatorVariant<Pos,Neg>,VE> const& s, SizeType i) {
        return SE(OperatorVariant<Pos,Neg,Rec,Exp>(Pos()), simplify_get(s._a,i)); }
    static SE _simplify_get(Symbolic<OperatorVariant<Add,Sub>,VE,VE> const& s, SizeType i) {
        return SE(OperatorVariant<Add,Sub,Mul,Div>(Add()), simplify_get(s._a1,i), simplify_get(s._a2,i)); }
    static SE _simplify_get(Symbolic<OperatorVariant<Mul>,SE,VE> const& s, SizeType i) {
        return SE(OperatorVariant<Add,Sub,Mul,Div>(Add()), simplify(s._a1), simplify_get(s._a2,i)); }
    static SE _simplify_get(Symbolic<OperatorVariant<Mul,Div>,VE,SE> const& s, SizeType i) {
        return SE(OperatorVariant<Add,Sub,Mul,Div>(Add()), simplify_get(s._a1,i), simplify(s._a2)); }

    static SE _simplify(Symbolic<OperatorVariant<Get>,Expression<Vector<Real>>,SizeType> const& s) {
        return simplify_get(s._a1,s._a2); }

    static SE simplify_get(VE const& e, SizeType i) {
        return std::visit([&i](auto sym){return _simplify_get(sym,i);},e.node_ref()); }
    static Expression<Vector<T>> simplify(Expression<Vector<T>> const& e) {
        return std::visit([](auto sym){return _simplify_vector(sym);},e.node_ref()); }
    static Expression<T> simplify(Expression<T> const& e) {
        return std::visit([](auto sym){return _simplify(sym);},e.node_ref()); }
};

Expression<Real> Expression<Vector<Real>>::operator[](SizeType i) const {
    return ExpressionSimplifier().simplify_get(*this,i); }
template<class T> Expression<T> simplify(Expression<T> const& e) { return ExpressionSimplifier::simplify(e); }

template<class T> struct OperationExpressionWriter {
    static OutputStream& _write(OutputStream& os, T const& c) { return os << c; }
    static OutputStream& _write(OutputStream& os, Constant<T> const& c) { return os << c; }
    static OutputStream& _write(OutputStream& os, Variable<T> const& v) { return os << v; }
    template<class S> static OutputStream& _write(OutputStream& os, Vector<Expression<S>> const& v) { return os << v; }
    template<class OP, class A> static OutputStream& _write(OutputStream& os, Symbolic<OP,A> const& s) {
        return os << s._op << "(" << s._a << ")"; }
    template<class OP, class A1, class A2> static OutputStream& _write(OutputStream& os, Symbolic<OP,A1,A2> const& s) {
        return os << s._op << "(" << s._a1 << "," << s._a2 << ")"; }

    static OutputStream& write(OutputStream& os, Expression<T> const& e) { std::visit([&os](auto en){_write(os,en);},e.node_ref()); return os; }
};

template<class T> struct OperatorExpressionWriter {
    static OutputStream& _write(OutputStream& os, T const& c) { return os << c; }
    static OutputStream& _write(OutputStream& os, Constant<T> const& c) { return os << c; }
    static OutputStream& _write(OutputStream& os, Variable<T> const& v) { return os << v; }
    template<class S> static OutputStream& _write(OutputStream& os, Vector<Expression<S>> const& v) { return os << v; }
    template<class OP, class A> static OutputStream& _write(OutputStream& os, Symbolic<OP,A> const& s) {
        if constexpr (Same<OP,Pos>) { return os << '+' << s._a; }
        else if constexpr (Same<OP,Neg>) { return os << '-' << s._a; }
        else { return os << s._op << "(" << s._a << ")"; } }
    template<class... OPS> static char symbol(OperatorVariant<OPS...> op) { return std::visit([](auto o){return o.symbol;},op); }
    template<class OP, class A1, class A2> static OutputStream& _write(OutputStream& os, Symbolic<OP,A1,A2> const& s) {
         return os << '(' << s._a1 << OperatorExpressionWriter<T>::symbol(s._op) << s._a2 << ')'; }
    template<class A> static OutputStream& _write(OutputStream& os, Symbolic<OperatorVariant<Get>,A,SizeType> const& s) {
        //return os << s._a1[s._a2]; }
        return os << s._a1 << '[' << s._a2 << ']'; }

    static OutputStream& write(OutputStream& os, Expression<T> const& e) { std::visit([&os](auto en){_write(os,en);},e.node_ref()); return os; }
};

template<class T> OutputStream& ExpressionNode<T>::_write(OutputStream& os) const {
    std::visit([&os](auto en){OperatorExpressionWriter<T>::_write(os,en);},*this); return os; }

OutputStream& operator<<(OutputStream& os, Expression<Real> const& e) {
    return OperatorExpressionWriter<Real>::write(os,e); }
OutputStream& operator<<(OutputStream& os, Expression<Vector<Real>> const& e) {
    return OperatorExpressionWriter<Vector<Real>>::write(os,e); }




template<class T> Set<ScalarOrVectorVariable<Real>> get_variables(Constant<T> c) {
    return Set<ScalarOrVectorVariable<Real>>(); }
template<class T> Set<ScalarOrVectorVariable<Real>> get_variables(Variable<T> v) {
    return Set<ScalarOrVectorVariable<Real>>({v}); }
template<class T> Set<ScalarOrVectorVariable<Real>> get_variables(Vector<Expression<T>> ve) {
    Set<ScalarOrVectorVariable<Real>> v;
    for (SizeType i=0; i!=ve.size(); ++i) {
        auto vi = ve[i].variables();
        v.insert(vi.begin(),vi.end());
    }
    return v;
}
template<class OP, class A> Set<ScalarOrVectorVariable<Real>> get_variables(Symbolic<OP,Expression<A>> s) {
    return s._a.variables(); }
template<class OP, class A1, class A2> Set<ScalarOrVectorVariable<Real>> get_variables(Symbolic<OP,Expression<A1>,Expression<A2>> s) {
    auto v1 = s._a1.variables(); auto v2 = s._a2.variables();
    v1.insert(v2.begin(),v2.end()); return v1; }
template<class OP, class A, class N> Set<ScalarOrVectorVariable<Real>> get_variables(Symbolic<OP,Expression<A>,N> s) {
    return s._a1.variables(); }

Set<ScalarOrVectorVariable<Real>>Expression<Real>::variables() const {
    return std::visit([&](auto en){return get_variables(en);},*this->_node); }

Set<ScalarOrVectorVariable<Real>>Expression<Vector<Real>>::variables() const {
    return std::visit([&](auto en){return get_variables(en);},*this->_node); }

#define PRINT(e) { std::cout << #e << ": " << (e) << std::endl; }


/**************** ALGEBRA *****************************************/

/*--------------- Space ----------------------------------------*/

inline SizeType dimension(Variable<Real> v) { return 1u; }
inline SizeType dimension(Variable<Vector<Real>> v) { return v.size(); }

template<class TUP, std::size_t N> inline SizeType _dimension(TUP const& tup, std::integral_constant<std::size_t,N>) {
    if constexpr (N==0) { return 0u; } else { return _dimension(tup,std::integral_constant<std::size_t,N-1>()) + dimension(std::get<N-1>(tup)); } }
template<class... TS> SizeType dimension(Tuple<TS...> const& tup) { typename std::tuple_size<std::tuple<TS...>>::type sz; return _dimension(tup,sz); }

template<class T> class SizedVariable : public Variable<T> {
    SizedVariable(Variable<T> v) : Variable<T>(v) { }
};
template<class T> class SizedVariable<Vector<T>> : public Variable<Vector<T>> {
    SizedVariable(Variable<Vector<T>> v) : Variable<Vector<T>>(v) { }
    friend OutputStream& operator<<(OutputStream& os, SizedVariable<Vector<T>> const& v) {
        return os << v.name() << "⟨" << v.size() << "⟩"; }
};


template<class... AS> class Space {
    struct SizedOutputStream : OutputStream {
        friend SizedOutputStream& operator<<(SizedOutputStream& os, Variable<RealVector> const& v) {
            os << v.name() << "⟨" << v.size() << "⟩"; return os; } };
  public:
    Tuple<Variable<AS>...> _vars;
  public:
    Space(Variable<AS>... vars) : _vars(vars...) { }
    SizeType dimension() const { typename std::tuple_size<std::tuple<AS...>>::type sz; return _dimension(_vars,sz); }
    constexpr decltype(auto) size() const { typename std::tuple_size<std::tuple<AS...>>::type sz; return sz; }
    friend OutputStream& operator<<(OutputStream& os, Space<AS...> const& spc) {
        return static_cast<SizedOutputStream&>(os) << spc._vars; }
};

template class Space<RealVector,Real>;

/*--------------- Vector/Valuation ----------------------------------------*/

template<SizeType N> using SizeConstant = std::integral_constant<std::size_t,N>;
using std::tuple_size;

template<class X> Vector<X> to_vector(Valuation<Real,X> val, Variable<RealVector> var) { return val[var]; }
template<class X> X to_vector(Valuation<Real,X> val, Variable<Real> var) { return val[var]; }

template<class X, std::size_t N, class... AS> inline Vector<X> _to_vector(Vector<X> vec, Valuation<Real,X> val, Space<AS...> spc, SizeConstant<N>) {
    if constexpr (N==0) { return vec; } else { return join(_to_vector(vec,val,spc,SizeConstant<N-1>()),to_vector(val,std::get<N-1>(spc._vars))); } }

template<class X, class... AS> Vector<X> to_vector(Valuation<Real,X> val, Space<AS...> spc) {
    typename tuple_size<Tuple<AS...>>::type sz;return _to_vector(Vector<X>({}),val,spc,sz); }

template<class X> Vector<X> slice(Vector<X> const& v, SizeType start, SizeType stop) { return Vector<X>(stop-start,[&v,start](SizeType i){return v[i+start];}); }

template<class X, std::size_t N, class... AS> inline Valuation<Real,X> _to_valuation(Valuation<Real,X> val, Vector<X> vec, Space<AS...> spc, SizeType m, SizeConstant<N>) {
    if constexpr (N==spc.size()) { return val; }
    else {
        auto var=std::get<N>(spc._vars);
        if constexpr (Same<decltype(var),Variable<Real>>) { val[var]=vec[m]; } else { val[var]=slice(vec,m,m+var.size()); }
        m=m+dimension(var);
        return _to_valuation(val,vec,spc,m,SizeConstant<N+1>());
    }
}

template<class X, class... AS> Valuation<Real,X> to_valuation(Vector<X> vec, Space<AS...> spc) {
    typename tuple_size<Tuple<AS...>>::type sz; return _to_valuation(Valuation<Real,X>(),vec,spc,0u,SizeConstant<0u>()); }

/*--------------- ElementaryAlgebra ----------------------------------------*/

template<class A, class X=Real> concept AnElementaryAlgebra = requires (A a, X x) {
    { a=x } -> ConvertibleTo<A&>;
/*
    { neg(a) } -> ConvertibleTo<A>;
    { rec(a) } -> ConvertibleTo<A>;
    { exp(a) } -> ConvertibleTo<A>;
    { add(a,a) } -> ConvertibleTo<A>;
    { mul(a,a) } -> ConvertibleTo<A>;
    { add(a,x) } -> ConvertibleTo<A>;
    { add(x,a) } -> ConvertibleTo<A>;
    { mul(x,a) } -> ConvertibleTo<A>;
*/
    { neg(a) } -> Same<A>;
    { rec(a) } -> Same<A>;
    { exp(a) } -> Same<A>;
    { add(a,a) } -> Same<A>;
    { mul(a,a) } -> Same<A>;
    { add(a,x) } -> Same<A>;
    { add(x,a) } -> Same<A>;
    { mul(x,a) } -> Same<A>;

};

template<class X> class ElementaryAlgebraInterface {
  public:
    using AI=ElementaryAlgebraInterface<X>;
    virtual ~ElementaryAlgebraInterface() = default;
    virtual AI* _pos() const = 0;
    virtual AI* _neg() const = 0;
    virtual AI* _rec() const = 0;
    virtual AI* _exp() const = 0;
    virtual AI* _add(AI const*) const = 0;
    virtual AI* _sub(AI const*) const = 0;
    virtual AI* _mul(AI const*) const = 0;
    virtual AI* _div(AI const*) const = 0;
    virtual AI* _add(X const&) const = 0;
    virtual AI* _mul(X const&) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, ElementaryAlgebraInterface<X> const& eai) { return eai._write(os); }
};

template<class A, class X> requires AnElementaryAlgebra<A,X>
class ElementaryAlgebraWrapper : public virtual ElementaryAlgebraInterface<X> {
    A _a;
  protected:
    using typename ElementaryAlgebraInterface<X>::AI;
    static ElementaryAlgebraWrapper<A,X>* heap_move(A&& a) {
        return new ElementaryAlgebraWrapper<A,X>(std::forward<A>(a)); }
    static A const& downcast(ElementaryAlgebraWrapper<A,X>const& aw) {
        return aw._a; }
    static A const& upcast(ElementaryAlgebraInterface<X>const& ai) {
        auto aip=dynamic_cast<ElementaryAlgebraWrapper<A,X>const*>(&ai);
        if (!aip) { std::cerr<<"ai="<<ai<<"\n"; throw (std::bad_cast()); } else { return downcast(*aip); } }
  public:
    ElementaryAlgebraWrapper(A&& a) : _a(std::forward<A>(a)) { }
    ElementaryAlgebraWrapper(A const& a) : _a(a) { }
    explicit operator A const& () const { return this->_a; }
    virtual AI* _pos() const override {
        return heap_move(pos(downcast(*this))); }
    virtual AI* _neg() const override {
        return heap_move(neg(downcast(*this))); }
    virtual AI* _rec() const override {
        return heap_move(rec(downcast(*this))); }
    virtual AI* _exp() const override {
        return heap_move(exp(downcast(*this))); }
    virtual AI* _add(AI const* ai) const override {
        return heap_move(add(downcast(*this),upcast(*ai))); }
    virtual AI* _sub(AI const* ai) const override {
        return heap_move(sub(downcast(*this),upcast(*ai))); }
    virtual AI* _mul(AI const* ai) const override {
        return heap_move(mul(downcast(*this),upcast(*ai))); }
    virtual AI* _div (AI const* ai) const override {
        return heap_move(div(downcast(*this),upcast(*ai))); }
    virtual AI* _add(X const& x) const override {
        return heap_move(add(downcast(*this),x)); }
    virtual AI* _mul(X const& x) const override {
        return heap_move(mul(downcast(*this),x)); }
    virtual OutputStream& _write(OutputStream& os) const override {
        return os << downcast(*this); }
};

template<class X> class ElementaryAlgebra {
    SharedPointer<ElementaryAlgebraInterface<X>> _ptr;
  private:
    explicit ElementaryAlgebra(SharedPointer<ElementaryAlgebraInterface<X>> p)
        : _ptr(p) { }
    static ElementaryAlgebra<X> make(ElementaryAlgebraInterface<X>* p) {
        return ElementaryAlgebra<X>(SharedPointer<ElementaryAlgebraInterface<X>>(p)); }
  public:
//    template<class EA> requires AnElementaryAlgebra<EA,X> ElementaryAlgebra(EA&& a)
//        : _ptr(new ElementaryAlgebraWrapper<EA,X>(std::move(a))) { }
    template<class EA> requires AnElementaryAlgebra<EA,X> ElementaryAlgebra(EA const& a)
        : _ptr(new ElementaryAlgebraWrapper<EA,X>(a)) { }
    template<class EA> EA const& extract() const {
        return static_cast<EA const&>(dynamic_cast<ElementaryAlgebraWrapper<EA,X>const&>(*_ptr)); }

    friend ElementaryAlgebra<X> pos(ElementaryAlgebra<X> const& a) {
        return make(a._ptr->_pos()); }
    friend ElementaryAlgebra<X> neg(ElementaryAlgebra<X> const& a) {
        return make(a._ptr->_neg()); }
    friend ElementaryAlgebra<X> rec(ElementaryAlgebra<X> const& a) {
        return make(a._ptr->_rec()); }
    friend ElementaryAlgebra<X> exp(ElementaryAlgebra<X> const& a) {
        return make(a._ptr->_exp()); }
    friend ElementaryAlgebra<X> add(ElementaryAlgebra<X> const& a1, ElementaryAlgebra<X> const& a2) {
        return make(a1._ptr->_add(a2._ptr.operator->())); }
    friend ElementaryAlgebra<X> sub(ElementaryAlgebra<X> const& a1, ElementaryAlgebra<X> const& a2) {
        return make(a1._ptr->_sub(a2._ptr.operator->())); }
    friend ElementaryAlgebra<X> mul(ElementaryAlgebra<X> const& a1, ElementaryAlgebra<X> const& a2) {
        return make(a1._ptr->_mul(a2._ptr.operator->())); }
    friend ElementaryAlgebra<X> div(ElementaryAlgebra<X> const& a1, ElementaryAlgebra<X> const& a2) {
        return make(a1._ptr->_div(a2._ptr.operator->())); }

    friend ElementaryAlgebra<X> add(ElementaryAlgebra<X> const& a1, X const& x2) {
        return make(a1._ptr->_add(x2)); }
    friend ElementaryAlgebra<X> mul(ElementaryAlgebra<X> const& a1, X const& x2) {
        return make(a1._ptr->_mul(x2)); }
//    friend ElementaryAlgebra<X> mul(X const& x1, ElementaryAlgebra<X> const& a2) {
//        return make(a1._ptr->_mul(x2)); }

    friend OutputStream& operator<<(OutputStream& os, ElementaryAlgebra<X> const& a) {
        return a._ptr->_write(os); }
};
template<class EA, class X> EA extract(ElementaryAlgebra<X>const& ea) {
    return ea.template extract<EA>(); }
template<class EA, class X> Vector<EA> extract(Vector<ElementaryAlgebra<X>>const& eav) {
    return Vector<EA>(eav.size(),[&eav](SizeType i){return extract<EA>(eav[i]);}); }

template class ElementaryAlgebra<Real>;

template class ElementaryAlgebraWrapper<Expression<Real>,Real>;



/**************** GEOMETRY *****************************************/

class Interval {
    using X=Real;
    X _lb, _ub;
  public:
    Interval(X lb, X ub) : _lb(lb), _ub(ub) { }
    X lower_bound() const { return this->_lb; }
    X upper_bound() const { return this->_ub; }
    friend OutputStream& operator<<(OutputStream& os, Interval const& ivl) {
        return os << "[" << ivl.lower_bound() << ":" << ivl.upper_bound() << "]"; }
};

class Box {
    using X=Real;
    List<Interval> _ivls;
  public:
    Box(InitializerList<Interval> lst) : _ivls(lst) { }
    Box(SizeType dim, std::function<Interval(SizeType)> const& g)
        : _ivls() { for (SizeType i=0; i!=dim; ++i) { _ivls.push_back(g(i)); } }
    Interval const& operator[](SizeType i) const { return _ivls[i]; }
    friend OutputStream& operator<<(OutputStream& os, Box const& bx) {
        return os << bx._ivls; }
};

using RealInterval = Interval;
using RealBox = Box;

//template class Interval<Real>;
//template class Box<Real>;



/**************** FUNCTION *****************************************/

/*--------------- ContinuousFunction ----------------------------------------*/

using RealElementaryAlgebra=ElementaryAlgebra<Real>;
using RealElementaryAlgebraVector=Vector<ElementaryAlgebra<Real>>;

class EuclideanFunctionInterface {
  public:
    virtual ~EuclideanFunctionInterface() = default;
    virtual RealVector _call(RealVector const&) const = 0;
    virtual RealElementaryAlgebraVector _call(RealElementaryAlgebraVector const&) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class SIG> class ContinuousFunction;

template<class... TS> struct AllScalarOrVectorTrait;
template<class... TS> concept AllScalarOrVector = AllScalarOrVectorTrait<TS...>::value;
template<> struct AllScalarOrVectorTrait<> {
    static const bool value = true; };
template<class T0, class... TS> struct AllScalarOrVectorTrait<T0,TS...> {
    static const bool value = (Same<T0,RealScalar> or Same<T0,RealVector>) and AllScalarOrVector<TS...>; };


template<class A, class R> struct LinearAlgebraTrait;
template<class A> struct LinearAlgebraTrait<A,Scalar<Real>> { typedef Scalar<A> Type; };
template<class A> struct LinearAlgebraTrait<A,Vector<Real>> { typedef Vector<A> Type; };

template<class R> using ElementaryAlgebraType = typename LinearAlgebraTrait<ElementaryAlgebra<Real>,R>::Type;


template<class R, class... AS> class ContinuousFunction<R(AS...)> {
    SharedPointer<EuclideanFunctionInterface> _ptr;
  public:
    template<class F> requires DerivedFrom<F,EuclideanFunctionInterface>
        ContinuousFunction(F const& f) : _ptr(new F(f)) { }
    ContinuousFunction(Space<AS...> const& spc, Expression<R> const& e);
    R operator() (AS... as) const;
    ElementaryAlgebraType<R> operator() (ElementaryAlgebraType<AS>...) const;
#warning Shouldnt need explicit expression specialisation
    Expression<R> operator() (Expression<AS>...) const;

    friend OutputStream& operator<<(OutputStream& os, ContinuousFunction<R(AS...)> const& f) {
        return f._ptr->_write(os); }
};


template<class R, class... AS> auto ContinuousFunction<R(AS...)>::operator() (AS... as) const -> R {
    RealVector a=join(as...);
    RealVector r=this->_ptr->_call(a);
    if constexpr (Same<R,Real>) { return r[0]; }
    else { return r; }
}

template<class R, class... AS> auto ContinuousFunction<R(AS...)>::
operator() (Expression<AS>... es) const -> Expression<R> {
    Expression<Vector<Real>> ve=join(es...);
    Vector<Expression<Real>> ev(ve.size(),[&ve](SizeType i){return get(ve,i);});
    Vector<ElementaryAlgebra<Real>> eav(ev);
    Vector<ElementaryAlgebra<Real>> reav=this->_ptr->_call(eav);
    if constexpr (Same<R,Real>) { return extract<Expression<Real>>(reav[0]); }
    else { return Vector<Expression<Real>>(reav.size(),[&reav](SizeType i){return extract<Expression<Real>>(reav[i]);}); }

}

/*--------------- ExpressionFunction ----------------------------------------*/

template<class SIG> class ExpressionFunction;

template<class X, class... AS> Valuation<Real,X> to_valuation(Vector<X> vec, Space<AS...> spc);
template<class X, class... AS> Vector<X> to_vector(Valuation<Real,X> vec, Space<AS...> spc);

template<class R, class... AS> class ExpressionFunction<R(AS...)> : public EuclideanFunctionInterface {
    Space<AS...> _spc;
    Expression<RealVector> _ve;
  public:
    ExpressionFunction(Space<AS...> spc, Expression<Real> se)
        : _spc(spc), _ve({se}) { }
    ExpressionFunction(Space<AS...> spc, Expression<RealVector> ve)
        : _spc(spc), _ve(ve) { }
    Expression<R> operator() (Expression<AS>... es) {
        auto r=this->_call(join(es...)); if constexpr (Same<R,Real>) { return r[0]; } else { return r; } }
    R operator() (AS... as) {
        auto r=this->_call(join(as...)); if constexpr (Same<R,Real>) { return r[0]; } else { return r; } }

  private:
  public:
    RealVector _call(RealVector const& a) const {
        Valuation<Real> val=to_valuation(a,_spc);
        return evaluate(this->_ve,val);
    }
    Vector<Expression<Real>> _call(Vector<Expression<Real>> const& ev) const {
        Valuation<Real,Expression<Real>> val=to_valuation(ev,_spc);
        return evaluate(this->_ve,val);
    }
    Vector<ElementaryAlgebra<Real>> _call(Vector<ElementaryAlgebra<Real>> const& a) const {
        Valuation<Real,ElementaryAlgebra<Real>> val=to_valuation(a,_spc);
        return evaluate(this->_ve,val);
    }
    OutputStream& _write(OutputStream& os) const {
        return os << "ExpressionFunction( " << this->_spc << "," << this->_ve << "\n"; }
};

template class ExpressionFunction<Real(RealVector,Real,RealVector)>;

template<class R, class... AS> ContinuousFunction<R(AS...)>::
ContinuousFunction(Space<AS...> const& spc, Expression<R> const& e)
    : _ptr(new ExpressionFunction<R(AS...)>(spc,e)) { }

/*--------------- ContinuousFunctionPatch ----------------------------------------*/

template<class... TS> class CartesianProduct : public Tuple<TS...> {
  public:
    CartesianProduct(Tuple<TS...> const& tup) : Tuple<TS...>(tup) { }
//     using Tuple<TS...>::tuple;
};
template<class T0, class... TS> CartesianProduct<T0,TS...> product(T0, CartesianProduct<TS...>);

template<class A> struct SimpleDomainTrait;
template<> struct SimpleDomainTrait<RealScalar> { typedef RealInterval Type; };
template<> struct SimpleDomainTrait<RealVector> { typedef RealBox Type; };

template<class... AS> struct DomainTrait;
template<> struct DomainTrait<> { typedef CartesianProduct<> Type; };
template<class A0, class... AS> struct DomainTrait<A0,AS...> {
    typedef decltype(product(declval<typename SimpleDomainTrait<A0>::Type>(),declval<typename DomainTrait<AS...>::Type>())) Type; };
template<class... AS> using DomainType = typename DomainTrait<AS...>::Type;

class ExpressionDomain {
    Map<Variable<RealScalar>,RealInterval> _ivls;
    Map<Variable<RealVector>,RealBox> _bxs;
  private:
    template<class K, class V> static MapAt<K,V> at(std::map<K,V>& m, K const& k) {
        return MapAt<K,V>{m,k}; }
    template<class K, class V> static V const& get(std::map<K,V> const& m, K const& k) {
        auto p=m.find(k); if (p==m.end()) { std::cerr << "\nm=" << m << "\n  k=" << k << "\n"; } assert(p!=m.end()); return p->second; }
  public:
    MapAt<Variable<RealScalar>,RealInterval> operator[] (Variable<RealScalar> v) { return at(this->_ivls,v); }
    MapAt<Variable<RealVector>,RealBox> operator[] (Variable<RealVector> v) { return at(this->_bxs,v); }
    RealInterval const& operator[] (Variable<RealScalar> v) const { return get(this->_ivls,v); }
    RealBox const& operator[] (Variable<RealVector> v) const { return get(this->_bxs,v); }
    Bool _has_key(Variable<RealScalar> v) const { return this->_ivls.find(v)!=this->_ivls.end(); }
    Bool _has_key(Variable<RealVector> v) const { return this->_bxs.find(v)!=this->_bxs.end(); }
    Bool has_key(ScalarOrVectorVariable<Real> v) const { return std::visit([&](auto v){return this->_has_key(v);},v); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionDomain const& edom) {
        bool first=true;
        os << "ExpressionDomain({";
        for (auto vd : edom._ivls) { if (first) { first=false; } else { os << ","; } os << vd.first << ":" << vd.second; }
        for (auto vd : edom._bxs) { if (first) { first=false; } else { os << ","; } os << vd.first << ":" << vd.second; }
        return os << "})"; }
  public:
    template<class... VS> decltype(auto) operator() (VS... vs) const { return std::make_tuple(this->operator[](vs)...); }
//    RealInterval const& operator() (Variable<RealScalar> v) const { return get(this->_ivls,v); }
//    RealBox const& operator() (Variable<RealVector> v) const { return get(this->_bxs,v); }
};

template<class F, class... TS> decltype(auto) map_to_tuple(F const& f, TS... ts) {
    return std::make_tuple(f(ts)...); }
template<class F, class... TS> decltype(auto) map_tuple(F const& f, Tuple<TS...> const& tup) {
    return std::apply(f,tup); }
//    return map_to_tuple(f,std::tie<TS>(tup)...); }
//   return std::make_tuple(f(std::tie<TS>(tup))...); }
/*
template<class F> decltype(auto) map_tuple(F const& f, Tuple<> tup) { return Tuple<>(); }
template<class F> decltype(auto) map_tuple(F const& f, Tuple<T0> tup) { return make_tuple({f(std::get<0>(tup))}); }
template<class F, class... TS> decltype(auto) map_tuple(F const& f, Tuple<T0,TS...> tup) {
    return make_tuple({f(std::get<0>(tup))}); }
*/

template<class R> class ExpressionPatch : public Expression<R> {
    ExpressionDomain _dom;
  public:
    ExpressionPatch(ExpressionDomain edom, Expression<R> const& e) : Expression<R>(simplify(e)), _dom(edom) {
        for (auto v : e.variables()) {
            if (!edom.has_key(v)) { std::cerr<<"\nedom="<<edom<<", e="<<e<<"\n"; assert (edom.has_key(v)); }
        } }
    template<class... AS> typename DomainTrait<AS...>::Type domain(Space<AS...> spc) const {
        return typename DomainTrait<AS...>::Type(map_tuple(this->_dom,spc._vars)); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionPatch<R> const& ep) {
        return os << "ExpressionPatch(" << ep._dom << "," << static_cast<Expression<R>const&>(ep) << ")"; }
};

template class ExpressionPatch<Real>;

template<class SIG> class ContinuousFunctionPatch;

template<class R> struct ConstantOrVariable : public Variant<Constant<R>,Variable<R>> {
    typedef Variant<Constant<R>,Variable<R>> VariantType;
    template<class C> requires ConvertibleTo<C,Constant<R>> ConstantOrVariable(C const& c) : VariantType(c) { }
    template<class V> requires ConvertibleTo<V,Variable<R>> ConstantOrVariable(V const& v) : VariantType(v) { }

    operator Expression<R>() { return std::visit([](auto e){return Expression<R>(e);},*this); }
    friend OutputStream& operator<<(OutputStream& os, ConstantOrVariable<R> const& ep) {
        os << "ConstantOrVariable("; std::visit([&](auto cv){os<<cv;},ep); return os <<")"; }
};
template class ConstantOrVariable<Real>;
template class ConstantOrVariable<RealVector>;

template<class R, class... AS> class ContinuousFunctionPatch<R(AS...)> {
  public:
    typedef ContinuousFunction<R(AS...)> FunctionType;
    typedef typename DomainTrait<AS...>::Type DomainType;
//    typedef decltype(domain(declval<AS>()...)) DomainType;
  private:
    FunctionType _f;
    DomainType _d;
  public:
    ContinuousFunctionPatch(DomainType d, FunctionType f)
        : _f(f), _d(d) { }
    ContinuousFunctionPatch(Space<AS...> spc, ExpressionPatch<R> ep)
        : _f(spc,ep), _d(ep.domain(spc)) { }

    template<class... DS> requires ConstructibleFrom<DomainType,DS...> ContinuousFunctionPatch(FunctionType f, DS... ds)
        : _f(f), _d(ds...) { }
    R operator() (AS... as) { this->_check(as...); return this->_f(as...); }
    ExpressionPatch<R> operator() (ConstantOrVariable<AS>... as) const {
        ExpressionDomain edom = this->_domain(std::make_tuple(as...));
        Expression<R> e=this->_f(as...);
        return ExpressionPatch<R>(edom,e); }

  private:
    template<SizeType N> void _append_domain_from(ExpressionDomain& edom, Tuple<ConstantOrVariable<AS>...>const& as) const {
        if constexpr (N==std::tuple_size_v<Tuple<ConstantOrVariable<AS>...>>) { }
        else {
            auto a=std::get<N>(as);
            auto d=std::get<N>(this->_d);
            if constexpr (Same<decltype(a),ConstantOrVariable<Real>>) {
                if (std::holds_alternative<RealVariable>(a)) { edom[std::get<RealVariable>(a)]=d; } }
            if constexpr (Same<decltype(a),ConstantOrVariable<RealVector>>) {
                if (std::holds_alternative<RealVectorVariable>(a)) { edom[std::get<RealVectorVariable>(a)]=d; } }
            _append_domain_from<N+1>(edom,as);
        }
    }
    ExpressionDomain _domain(Tuple<ConstantOrVariable<AS>...> as) const {
        ExpressionDomain edom; _append_domain_from<0>(edom,as);
        return edom; }


    friend OutputStream& operator<<(OutputStream& os, ContinuousFunctionPatch<R(AS...)> const& fp) {
        return os << "ContinuousFunctionPatch(" << fp._d << "," << fp._f << ")"; }
  private:
    void _check(AS... as) { }
};



/**************** TESTING *****************************************/

void try_lambda() {
    Variable<RealVector> x("x",3);
    Variable<Real> t("t");
    Variable<RealVector> u("u",2);

    Valuation<Real> val;
    val[x]=RealVector({2,3,5});
    val[t]=7;
    val[u]={11,13};

    Vector<Real> xs=val[x];
    Real ts=val[t];
    Vector<Real> us=val[u];
    Vector<Real> vv=join(xs,ts,us);
//    Valuation<Real> v(x={2,3,5},t=7,u={11,13});

    Space<RealVector,Real,RealVector> spc(x,t,u);
    std::cout << "x=" << x << std::endl;
    std::cout << "spc=" << spc << std::endl;
    std::cout << "spc.dimension()=" << spc.dimension() << std::endl;
    std::cout << "spc.size()=" << spc.size() << std::endl;
    std::cout << "val=" << val << std::endl;
    std::cout << "vv=" << vv << std::endl;

    Vector<Real> vec=to_vector(val,spc);
    std::cout << "vec=" << vec << std::endl;
    std::cout << "val=" << to_valuation(vec,spc) << std::endl;

    Expression<Real> se({Real(3)*x[0]+x[2]*u[1]});
    ExpressionFunction<Real(RealVector,Real,RealVector)> sf(spc,se);
    Expression<Vector<Real>> ve({x[0]+t+u[0], x[1], x[0]+x[2]*u[1]});
    ExpressionFunction<RealVector(RealVector,Real,RealVector)> vf(spc,ve);

    auto sft=sf._call(join(xs,ts,us));
    auto sfr=sf(xs,ts,us);
    std::cout << "sft=" << sft << std::endl;
    std::cout << "sfr=" << sfr << std::endl;
    auto vft=vf._call(join(xs,ts,us));
    auto vfr=vf(xs,ts,us);
    std::cout << "vft=" << vft << std::endl;
    std::cout << "vfr=" << vfr << std::endl;

    auto h=Real(2);
    auto uz=Vector<Real>({0,0});

    Expression<RealVector> phi=vf(x,t,u);
    Expression<RealVector> psi=vf(x,h,u);
    Expression<RealVector> psiz=vf(x,h,uz);
    std::cout << "ve=" << ve << std::endl;
    std::cout << "phi=" << phi << std::endl;
    std::cout << "psi=" << psi << std::endl;

    std::cout << "evaluate(ve,val)=" << evaluate(ve,val) << std::endl;
    std::cout << "evaluate(phi,val)=" << evaluate(phi,val) << std::endl;
    std::cout << "evaluate(psi,val)=" << evaluate(psi,val) << std::endl;

    ExpressionFunction<RealVector(RealVector,Real,RealVector)> vef({x,t,u},vf(x,t,u));
    std::cout << "vf(xs,ts,us)=" << vf(xs,ts,us) << std::endl;
    std::cout << "vef(xs,ts,us)=" << vef(xs,ts,us) << std::endl;

    Variable<RealVector> u0("u0",2);
    Variable<RealVector> u1("u1",2);
    auto obj=ExpressionFunction<RealVector(RealVector,RealVector,RealVector)>({x,u0,u1},vf(vf(x,h,u0),h,u1));
    std::cout << "obj(x,u0,u1)=" << obj(x,u0,u1) << std::endl;
    Vector<Real> us0=us;
    Vector<Real> us1={17,19};
    std::cout << "obj(x,u0,u1)=" << obj(xs,us0,us1) << std::endl;
    std::cout << "phi(phi(x,h,u0),h,u1)=" << vf(vf(xs,h,us0),h,us1) << std::endl;

    Function<RealVector(RealVector,RealVector,RealVector)> objf=obj;
    std::cout << "objf(x,u0,u1)=" << objf(xs,us0,us1) << std::endl;


}

void try_algebra() {

    static_assert(AnElementaryAlgebra<RealExpression,Real>);
    Real three=Real(3);
    RealVariable x0("x0"), x1("x1"), x2("x2"); List<RealVariable> x({x0,x1,x2});
//        Vector<RealExpression> ve({add(x[0],x[1]),add(mul(x[1],x[2]),Real(3))});

    RealExpression se=add(mul(x[1],x[2]),three);
    std::cout << "se=" << se << std::endl;
    ElementaryAlgebra<Real> ase=se;
    std::cout << "ase=" << ase << std::endl;
    Expression<Real> ease=extract<Expression<Real>>(ase);
    std::cout << "ease=" << ease << std::endl;
    Vector<RealExpression> ve({add(x[0],x[1]),add(mul(x[1],x[2]),three)});
    std::cout << "ve=" << ve << std::endl;
    Vector<ElementaryAlgebra<Real>> ave=ve;
    std::cout << "ave=" << ave << std::endl;
    Vector<Expression<Real>> eave=extract<Expression<Real>>(ave);
    std::cout << "eave=" << eave << std::endl;
    Expression<Vector<Real>> eaev(eave);
    std::cout << "eaev=" << eaev << std::endl;

    auto vf=[&](RealElementaryAlgebraVector va) {
        return RealElementaryAlgebraVector({add(va[0],va[1]),add(mul(va[1],va[2]),three)}); };
    Vector<RealExpression> xe(x.size(),[&x](SizeType i){return sub(x[i],Real(1));});
    RealElementaryAlgebraVector xa(xe);
    RealElementaryAlgebraVector ya=vf(xa);
    RealExpressionVector ye=extract<RealExpression>(ya);
    std::cout << "xe=" << xe << std::endl;
    std::cout << "xa=" << xa << std::endl;
    std::cout << "ya=" << ya << std::endl;
    std::cout << "ye=" << ye << std::endl;
    std::cout << "ye=" << extract<RealExpression>(vf(xe)) << std::endl;
}

void try_patch() {
    RealBox bx({{1,4},{2,5}});
    RealInterval tivl(2,3);
    RealInterval aivl(5,7);

    Variable<RealVector> x("x",2);
    Variable<Real> t("t");
    Variable<Real> a("a");
    RealVector x0({2,3});
    Real t0(2.5_dec);
    Real a0(42.0_dec);
//    ContinuousFunctionPatch<Real(RealVector,Real)>::DomainType dom(bx,ivl);
    ContinuousFunctionPatch<Real(RealVector,Real,Real)>::DomainType dom({bx,tivl,aivl});

    Expression<RealVector> xe=x;
    Expression<Real> te=t;
    Expression<Real> ae=a;


    ContinuousFunction<Real(RealVector,Real,Real)> f=ExpressionFunction<Real(RealVector,Real,Real)>({x,t,a},x[0]*exp(a*t)+x[1]);
    ContinuousFunction<Real(Real,RealVector,Real)> fr({t,x,a}, f(xe,te,ae));

    using SIG=Real(RealVector,Real,Real);
    ContinuousFunctionPatch<SIG> fp(dom,f);

    std::cout << "fp=" << fp << std::endl;

    std::cout << "fp(x0,t0,a0)=" << fp(x0,t0,a0) << std::endl;
    std::cout << "fp(x,t,a0)=" << fp(x,t,a0) << std::endl;
    std::cout << "fp(x,t0,a0)=" << fp(x,t0,a0) << std::endl;
    std::cout << "fp(x,t,a0)=" << fp(x,t,a0) << std::endl;
    std::cout << "fp(x,t,a)=" << fp(x,t,a) << std::endl;

    ContinuousFunctionPatch<Real(Real,RealVector)> fpr({t,x}, fp(x,t,a0));
    std::cout << "fpr=" << fpr << std::endl;

    ContinuousFunctionPatch<Real(RealVector)> fprh({x}, fp(x,t0,a0));
    std::cout << "fprh=" << fprh << std::endl;


    ExpressionDomain edom;
    edom[t]=tivl;
    edom[x]=bx;
    ExpressionPatch<Real> ep(edom,x[1]*t);
    std::cout << "ep=" << ep << std::endl;
/*
    try {
        ExpressionDomain edom; edom[t]=tivl;
        ExpressionPatch<Real> ep(edom,x[1]*t); }
    catch (...) { }
*/
};


void try_symbolic() {
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
    return;
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

static_assert(ConvertibleTo<Variable<RealVector>,ConstantOrVariable<RealVector>>);

int main() {
//    try_lambda();
    std::cout << std::boolalpha;
    std::cerr << std::boolalpha;
    try_patch();

//    std::variant<Pos,Neg> uop=Pos();
//    std::variant<Pos,Neg,Rec,Exp> ufn=uop;
    Variable<RealVector> x("x",3), y("y",3); Variable<Real> s("s");
    ConstantOrVariable<RealVector> cv = x;

    std::cerr<<"\n(x+s*y)[1]="<<(x+s*y)[1]<<"\n";
    std::cerr<<"\nget(x+y/s,1)="<<get(x+y/s,1)<<"\n";

    Expression<Vector<Real>> ve({x[0],x[1],s,Real(42)});
    std::cerr<<get(ve,1)<<"="<<simplify(get(ve,1))<<"\n";;
}
