/***************************************************************************
 *            korteweg_devries.cpp
 *
 *  Copyright  2019  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>

#include "config.hpp"

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"
#include "numeric/numeric.hpp"
#include "numeric/complex.hpp"
#include "algebra/expansion.inl.hpp"
//#include "algebra/expansion.tpl.hpp"
#include "function/fourier_polynomial.hpp"
#include "function/fourier_polynomial.tpl.hpp"
#include "function/taylor_function.hpp"

#include "ariadne.hpp"

#include "solvers/integral.hpp"

namespace Ariadne {

struct Coordinate { };

template<class T> using PropertiesOfType = decltype(properties(declval<T>()));

class UnivariateFunctionClass {
  public:
    typedef ValidatedTag Paradigm;
    UnivariateFunctionClass(Coordinate);
    friend UnivariateFunctionClass operator+(UnivariateFunctionClass,UnivariateFunctionClass);
    friend UnivariateFunctionClass operator*(UnivariateFunctionClass,UnivariateFunctionClass);
    friend UnivariateFunctionClass operator*(ValidatedNumber,UnivariateFunctionClass);
    friend UnivariateFunctionClass operator/(UnivariateFunctionClass,ValidatedNumber);

    friend UnivariateFunctionClass add(UnivariateFunctionClass,UnivariateFunctionClass);
    friend UnivariateFunctionClass sub(UnivariateFunctionClass,UnivariateFunctionClass);
    friend UnivariateFunctionClass mul(UnivariateFunctionClass,UnivariateFunctionClass);
    friend UnivariateFunctionClass div(UnivariateFunctionClass,UnivariateFunctionClass);
    friend UnivariateFunctionClass sqr(UnivariateFunctionClass);
    friend UnivariateFunctionClass sqrt(UnivariateFunctionClass);
    friend UnivariateFunctionClass exp(UnivariateFunctionClass);
    friend UnivariateFunctionClass sin(UnivariateFunctionClass);
    friend UnivariateFunctionClass cos(UnivariateFunctionClass);
    friend UnivariateFunctionClass atan(UnivariateFunctionClass);

};

ValidatedScalarMultivariateTaylorFunctionModelDP cast_positive(ValidatedScalarMultivariateTaylorFunctionModelDP f) { return f; }
decltype(auto) mag(ValidatedScalarMultivariateTaylorFunctionModelDP f) { return mag(norm(f)); }

//template<class MF>
class ValidatedScalarUnivariateTaylorFunctionModel {
    typedef ValidatedScalarMultivariateTaylorFunctionModelDP MF;

    typedef ValidatedScalarMultivariateFunction UnivariateFunctionType;
    //typedef ValidatedScalarUnivariateFunction UnivariateFunctionType;
    typedef ValidatedScalarMultivariateTaylorFunctionModelDP MultivariateFunctionModelType;
    typedef ValidatedScalarUnivariateTaylorFunctionModel UnivariateFunctionModelType;
    MF _mf;
  public:
    typedef IntervalDomainType DomainType;
    typedef ValidatedTag Paradigm;
    typedef DoublePrecision PrecisionType;
    typedef MultivariateFunctionModelType::PropertiesType SweeperType;
    typedef Tuple<DomainType,SweeperType> PropertiesType;

    ValidatedScalarUnivariateTaylorFunctionModel() : ValidatedScalarUnivariateTaylorFunctionModel(IntervalDomainType(-1,1),SweeperDP()) { }
    ValidatedScalarUnivariateTaylorFunctionModel(PropertiesType prp) : ValidatedScalarUnivariateTaylorFunctionModel(std::get<0>(prp),std::get<1>(prp)) { }
    ValidatedScalarUnivariateTaylorFunctionModel(IntervalDomainType dom, SweeperType prp) : _mf({dom},prp) { }
    ValidatedScalarUnivariateTaylorFunctionModel(IntervalDomainType dom, Coordinate, SweeperType prp)
        : ValidatedScalarUnivariateTaylorFunctionModel(MultivariateFunctionModelType::coordinate({dom},0u,prp)) { }

    static UnivariateFunctionModelType constant(IntervalDomainType dom, ValidatedNumber c, SweeperType prp) {
        return MF::constant({dom},{c,prp.precision()},prp); }
    static UnivariateFunctionModelType coordinate(IntervalDomainType dom, SweeperType prp) {
        return MF::coordinate({dom},0u,prp); }

    PrecisionType precision() const { return _mf.precision(); }
    SweeperType sweeper() const { return _mf.properties(); }
    PropertiesType properties() const { return make_tuple(this->domain(),this->sweeper()); }

    decltype(auto) terms() const { return _mf.expansion(); }

//    UnivariateFunctionWrapper<MF>(MultivariateFunctionModelType mf) : MF(mf) { }
    ValidatedScalarUnivariateTaylorFunctionModel(MultivariateFunctionModelType mf) : _mf(mf) { }

    friend UnivariateFunctionModelType operator+(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(f1._mf+f2._mf); }
    friend UnivariateFunctionModelType operator-(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(f1._mf-f2._mf); }
    friend UnivariateFunctionModelType operator*(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(f1._mf*f2._mf); }
    friend UnivariateFunctionModelType operator/(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(f1._mf/f2._mf); }
    friend UnivariateFunctionModelType& operator+=(UnivariateFunctionModelType& f1,UnivariateFunctionModelType f2) {
        f1._mf+=f2._mf; return f1; }
    friend UnivariateFunctionModelType& operator*=(UnivariateFunctionModelType& f1,UnivariateFunctionModelType f2) {
        f1._mf*=f2._mf; return f1; }
    friend UnivariateFunctionModelType operator+(UnivariateFunctionModelType f1, ValidatedNumber c2) {
        return cast_univariate(f1._mf+c2); }
    friend UnivariateFunctionModelType operator-(UnivariateFunctionModelType f1, ValidatedNumber c2) {
        return cast_univariate(f1._mf-c2); }
    friend UnivariateFunctionModelType operator*(UnivariateFunctionModelType f1, ValidatedNumber c2) {
        return cast_univariate(f1._mf*c2); }
    friend UnivariateFunctionModelType operator/(UnivariateFunctionModelType f1, ValidatedNumber c2) {
        return cast_univariate(f1._mf/c2); }
    friend UnivariateFunctionModelType operator+(ValidatedNumber c1, UnivariateFunctionModelType f2) {
        return cast_univariate(c1+f2._mf); }
    friend UnivariateFunctionModelType operator-(ValidatedNumber c1, UnivariateFunctionModelType f2) {
        return cast_univariate(c1-f2._mf); }
    friend UnivariateFunctionModelType operator*(ValidatedNumber c1, UnivariateFunctionModelType f2) {
        return cast_univariate(c1*f2._mf); }

    friend UnivariateFunctionModelType operator+(UnivariateFunctionType const& f1, UnivariateFunctionModelType f2) {
        return cast_univariate(f1+f2._mf); }
    friend UnivariateFunctionModelType operator+(UnivariateFunctionModelType f1, UnivariateFunctionType const& f2) {
        return cast_univariate(f1._mf+f2); }
    friend UnivariateFunctionModelType operator*(UnivariateFunctionType const& f1, UnivariateFunctionModelType f2) {
        return cast_univariate(f1*f2._mf); }
    friend UnivariateFunctionModelType operator*(UnivariateFunctionModelType f1, UnivariateFunctionType const& f2) {
        return cast_univariate(f1._mf*f2); }
    friend UnivariateFunctionModelType& operator+=(UnivariateFunctionModelType& f1, UnivariateFunctionType const& f2) {
        f1._mf+=f2; return f1; }
    friend UnivariateFunctionModelType& operator*=(UnivariateFunctionModelType& f1, UnivariateFunctionType const& f2) {
        f1._mf*=f2; return f1; }

    static UnivariateFunctionModelType cast_univariate(MultivariateFunctionModelType f) { return UnivariateFunctionModelType(std::move(f)); }

    //friend UnivariateFunctionModelType nul(UnivariateFunctionModelType f) { return cast_univariate(nul(f._mf)); }
    //friend UnivariateFunctionModelType nul(UnivariateFunctionModelType f);
    friend UnivariateFunctionModelType nul(UnivariateFunctionModelType f) { return f*0; }
    friend UnivariateFunctionModelType sqr(UnivariateFunctionModelType f) { return cast_univariate(sqr(f._mf)); }
    friend UnivariateFunctionModelType add(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(add(f1._mf,f2._mf)); }
    friend UnivariateFunctionModelType sub(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(sub(f1._mf,f2._mf)); }
    friend UnivariateFunctionModelType mul(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(mul(f1._mf,f2._mf)); }
    friend UnivariateFunctionModelType div(UnivariateFunctionModelType f1,UnivariateFunctionModelType f2) {
        return cast_univariate(div(f1._mf,f2._mf)); }
    friend UnivariateFunctionModelType sqrt(UnivariateFunctionModelType f) { return cast_univariate(sqrt(f._mf)); }
    friend UnivariateFunctionModelType exp(UnivariateFunctionModelType f) { return cast_univariate(exp(f._mf)); }
    friend UnivariateFunctionModelType log(UnivariateFunctionModelType f) { return cast_univariate(log(f._mf)); }
    friend UnivariateFunctionModelType sin(UnivariateFunctionModelType f) { return cast_univariate(sin(f._mf)); }
    friend UnivariateFunctionModelType cos(UnivariateFunctionModelType f) { return cast_univariate(cos(f._mf)); }
    friend UnivariateFunctionModelType tan(UnivariateFunctionModelType f) { return cast_univariate(tan(f._mf)); }
    friend UnivariateFunctionModelType atan(UnivariateFunctionModelType f) { return cast_univariate(atan(f._mf)); }

    friend decltype(auto) mag(UnivariateFunctionModelType const& f) { return mag(f._mf); }
    friend ValidatedKleenean operator>(UnivariateFunctionModelType f, Int c) { return f._mf.range()>c; }
    friend ValidatedKleenean operator<(UnivariateFunctionModelType f, Int c) { return f._mf.range()<c; }
    friend ValidatedNegatedSierpinskian operator==(UnivariateFunctionModelType,UnivariateFunctionModelType);


    friend ValidatedScalarUnivariateTaylorFunctionModel const& cast_positive(ValidatedScalarUnivariateTaylorFunctionModel const& f) { return f; }
    friend PropertiesType properties(ValidatedScalarUnivariateTaylorFunctionModel const& f) { return f.properties(); }

    IntervalDomainType domain() const { return this->_mf.domain()[0]; }
    template<class X> decltype(auto) operator() (X const& x) { return _mf(Vector<X>({x})); }
    template<class X> friend decltype(auto) evaluate(UnivariateFunctionModelType f, X const& x) { return f(x); }

    friend OutputStream& operator<<(OutputStream& os, ValidatedScalarUnivariateTaylorFunctionModel const& f) {
        Array<String> names={"t"}; os << "TF["<<f.domain()<<"]("; f._mf.expansion().write(os,names); return os<<")"; }
        //return os << f._mf; }
};

//typedef UnivariateFunctionWrapper<ValidatedScalarMultivariateTaylorFunctionModelDP> ValidatedScalarUnivariateTaylorFunctionModel;


using UFP = UnivariateFourierPolynomial<Complex<FloatDPBounds>>;

Complex<FloatDPBounds> operator/(Complex<FloatDPBounds> const& z1, Int const& x2) { return div(z1,FloatDPBounds(x2,dp)); }

template class Complex<ValidatedScalarUnivariateTaylorFunctionModel>;

//template class UnivariateFourierPolynomial<Complex<ValidatedScalarUnivariateTaylorFunctionModel>>;

using FXT = UnivariateFourierPolynomial<Complex<ValidatedScalarUnivariateTaylorFunctionModel>>;

template<class X> using Fourier = UnivariateFourierPolynomial<X>;

template<class X> decltype(auto) properties(Fourier<X> const& f) { return properties(f.zero_coefficient()); }


template<class T, class X> Fourier<X> B2(T t, Fourier<X> v1, Fourier<X> v2) {
    auto prp=properties(v1);
    Fourier<X> r(prp);
    auto prec=min(v1.precision(),v2.precision());
    Complex<FloatDPBounds> i(0,1,prec);

    for (auto t1 : v1.terms()) {
        short k1 = t1.index();
        assert(k1!=0);
        for (auto t2 : v2.terms()) {
            short k2 = t2.index();
            assert(k2!=0);
            short k = k1+k2;
            r[k] += exp(i*(3*k*k1*k2)*t)/(k1*k2)*t1.coefficient()*t2.coefficient();
        }
    }
    return r;
}



template<class T, class X> Fourier<X> R3(T t, Fourier<X> v1, Fourier<X> v2, Fourier<X> v3) {
    auto prop=properties(v1);
    auto prec=min(min(v1.precision(),v2.precision()),v3.precision());
    X z = X(prop);
    Fourier<X> r(z);
    Complex<FloatDPBounds> i(0,1,prec);

    std::cerr<<"v1="<<v1<<"\n";
    std::cerr<<"v1.terms().number_of_nonzeros()="<<v1.terms().number_of_nonzeros()<<"\n";
    for (auto t1 : v1.terms()) {
        short k1 = t1.index();
        assert(k1!=0);
        for (auto t2 : v2.terms()) {
            auto k2 = t2.index();
            for (auto t3 : v3.terms()) {
                auto k3 = t3.index();
                auto k = k1+k2+k3;

                std::cerr<<"k1,2,3="<<k1<<","<<k2<<","<<k3<<";"<<k<<"\n";
                auto p = exp(i*(3*(k1+k2)*(k2+k3)*(k3*k1))*t)/(k1)*((t1.coefficient()*t2.coefficient())*t3.coefficient());
                std::cerr<<"  p.nnz="<<p.real_part().terms().number_of_nonzeros()<<" "<<std::flush;
                r[k] += p;
                std::cerr<<"  nnz="<<r.terms().number_of_nonzeros()<<std::endl;
            }
        }
    }
    std::cerr<<"Done!\n";
    return r;
}

template<class T, class X> Fourier<X> W(T t, Fourier<X> u) {
    typedef BiUniIndex I;
    using Constants::i;
    Fourier<X>& v=u;
    for (auto ut : u.terms()) {
        I k = ut.index();
        X& uk = ut.coefficient();
        uk *= exp(i*(k*k*k)*t);
    }
    return std::move(v);
}

} // namespace Ariadne


#define ARIADNE_PRINT(expr) { std::cout << #expr << "=" << (expr) << "\n"; }

int main() {
    using namespace Ariadne;

    Complex<Integer> i(0,1);

    IntervalDomainType dom(-1,+1);
    DoublePrecision pr;
    SweeperDP swp=ThresholdSweeperDP(pr,1e-3);

    ValidatedScalarUnivariateTaylorFunctionModel akz(dom,swp);
    ARIADNE_PRINT(akz)
    ValidatedScalarUnivariateTaylorFunctionModel akx(dom,Coordinate(),swp);
    ARIADNE_PRINT(akx)
    Complex<ValidatedScalarUnivariateTaylorFunctionModel> ukz(akz,akz);
    ARIADNE_PRINT(ukz)

    UnivariateFourierPolynomial<Complex<ValidatedScalarUnivariateTaylorFunctionModel>> u(akz);
    ARIADNE_PRINT(u)
    ARIADNE_PRINT(akx*akx-1+i*akx)

    auto uk0=akz;
    auto uk1=akx;
    auto uk2=akx*akx-1+i*akx;
    auto uk3=akx*akx*akx/2-5*akx/3;
    //u.terms().append(0,uk0);
    u.terms().append(-1,uk1);
    ARIADNE_PRINT(u)
    //u[2]=uk2;
    u.terms().append(1,uk2);
    ARIADNE_PRINT(u)
    ARIADNE_PRINT(uk3)
    ARIADNE_PRINT(u.terms().at(2))
    u.terms().append(2,akx);
/*
    ARIADNE_PRINT(u.terms().at(3))
    BiUniIndex a=3;
    ARIADNE_PRINT(a)
    ARIADNE_PRINT(uk3)
    auto u3=u[3];
    ARIADNE_PRINT(u3._a);
    auto ua=u[a];
    ARIADNE_PRINT(ua._a);
    u3=uk3;
    u[a]=uk3;
    ARIADNE_PRINT(u.terms().at(3))
*/

    ARIADNE_PRINT(u)
    u.terms().sort();
    ARIADNE_PRINT(u)

    //ValidatedScalarUnivariateFunction t=ValidatedScalarUnivariateFunction::coordinate();
    ValidatedScalarMultivariateFunction t=ValidatedScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedScalarUnivariateTaylorFunctionModel tm=ValidatedScalarUnivariateTaylorFunctionModel::coordinate(dom,swp);

    tm*t;
    ARIADNE_PRINT(tm*t)

    ARIADNE_PRINT(W(t,u));
    ARIADNE_PRINT(B2(t,u,u))
    ARIADNE_PRINT(R3(t,u,u,u));

    return 0;

}
