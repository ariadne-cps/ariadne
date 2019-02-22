/***************************************************************************
 *            test_fourier_polynomial.cpp
 *
 *  Copyright 2008--18  Pieter Collins
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

#include <iostream>
#include <iomanip>
#include "config.hpp"

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"
#include "numeric/numeric.hpp"
#include "numeric/complex.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/differential.hpp"
#include "function/function.hpp"
#include "function/polynomial.hpp"

#include "function/fourier_polynomial.hpp"
#include "function/fourier_polynomial.tpl.hpp"

#include "algebra/sweeper.hpp"

#include "../test.hpp"

namespace Ariadne {
template<> String class_name<Complex<FloatDPApproximation>>() { return "ComplexFloatDPApproximation"; }
template<> String class_name<Complex<FloatMPApproximation>>() { return "ComplexFloatMPApproximation"; }
}

using std::cout; using std::cerr; using std::endl;
using namespace Ariadne;

extern template Ariadne::Nat Ariadne::Error<Ariadne::FloatMP>::output_places;
extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatMP>::output_places;

template<class I, class X> decltype(auto) sup_norm(Polynomial<I,X> const& p) {
    auto r=mag(p.expansion().zero_coefficient());
    for (auto term : p) { r += mag(term.coefficient()); }
    return r;
}

template<class X> class TestFourierPolynomial
{
    typedef PrecisionType<X> PR;
    PR pr;
  public:
    TestFourierPolynomial(PrecisionType<X> prec);
    Void test();
  private:
    Void test_concept();
    Void test_univariate();
    Void test_multivariate();
};

template<class X> TestFourierPolynomial<X>::TestFourierPolynomial(PrecisionType<X> prec)
    : pr(prec)
{
}

template<class X> Void TestFourierPolynomial<X>::test()
{
    FloatValue<PR>::set_output_places(18);
    FloatBounds<PR>::set_output_places(18);
    FloatApproximation<PR>::set_output_places(8);

    ARIADNE_TEST_CALL(test_univariate());
    ARIADNE_TEST_CALL(test_multivariate());
}

template<class X> Void TestFourierPolynomial<X>::test_concept()
{
    X c;
    UnivariateFourierPolynomial<X>* xp; UnivariateFourierPolynomial<X>& x=*xp;
    x=+x; x=-x; x=x+x; x=x-x; x=x*x;
    x=x+c; x=c+x; x=x-c; x=c-x; x=x*c; x=c*x; x=x/c;
    // x=sqr(x);
    x=pow(x,2u);


}


template<class X> Void TestFourierPolynomial<X>::test_univariate() {
#warning
    X eps(1e-12);
    typedef Complex<X> CX;
    ARIADNE_TEST_NAMED_CONSTRUCT(UnivariateFourierPolynomial<CX>,one,constant(1,pr));
    ARIADNE_TEST_NAMED_CONSTRUCT(UnivariateFourierPolynomial<CX>,x,basis(1,pr));
    ARIADNE_TEST_NAMED_CONSTRUCT(UnivariateFourierPolynomial<CX>,x2,basis(2,pr));
    ARIADNE_TEST_NAMED_CONSTRUCT(UnivariateFourierPolynomial<CX>,conj_x,basis(-1,pr));

    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_PRINT(x*x);
    ARIADNE_TEST_PRINT(x*x*x);
    ARIADNE_TEST_PRINT(x*x*x*x);

    ARIADNE_TEST_EQUAL(x*x,x2);

    ARIADNE_TEST_CONSTRUCT(X,y,(-0.75_dy,pr));

    ARIADNE_TEST_WITHIN(x(y),CX(cos(y),sin(y)),2*eps);
    ARIADNE_TEST_WITHIN((x*x)(y),CX(cos(2*y),sin(2*y)),2*eps);
    ARIADNE_TEST_WITHIN((x*x*x)(y),CX(cos(3*y),sin(3*y)),2*eps);
    ARIADNE_TEST_WITHIN((x*x*x*x)(y),CX(cos(4*y),sin(4*y)),2*eps);

    ARIADNE_TEST_WITHIN((2*x*x-one)(y),CX(2*cos(2*y)-1,2*sin(2*y)),2*eps);
    ARIADNE_TEST_WITHIN((2*x*x-CX(1,0,pr))(y),CX(2*cos(2*y)-1,2*sin(2*y)),2*eps);
    ARIADNE_TEST_WITHIN((2*x*x-1)(y),CX(2*cos(2*y)-1,2*sin(2*y)),2*eps);
    ARIADNE_TEST_WITHIN((3*x*x*x-5*x)(y),CX(3*cos(3*y)-5*cos(y),3*sin(3*y)-5*sin(y)),2*eps);
    ARIADNE_TEST_WITHIN((8*(sqr(x)-1)*sqr(x)+1)(y),CX(8*cos(4*y)-8*cos(2*y)-1,8*sin(4*y)-8*sin(2*y)),2*eps);

    UnivariatePolynomial<X> p=UnivariatePolynomial<X>::coordinate();
    ARIADNE_TEST_PRINT(p);

    ARIADNE_TEST_WITHIN(two_norm(3*x*x*x-4*x),5,2*eps);
    ARIADNE_TEST_WITHIN(sup_norm(3*x*x*x-4*x),7,2*eps);

}


template<class X> Void TestFourierPolynomial<X>::test_multivariate() {
/*
    ARIADNE_TEST_NAMED_CONSTRUCT(MultivariateFourierPolynomial<X>,one,constant(2u,1,pr));
    ARIADNE_TEST_NAMED_CONSTRUCT(MultivariateFourierPolynomial<X>,x,coordinate(2u,0u,pr));
    ARIADNE_TEST_NAMED_CONSTRUCT(MultivariateFourierPolynomial<X>,y,coordinate(2u,1u,pr));

    Vector<X> v({0.5,-0.75},pr);
    ARIADNE_TEST_EQUALS(x(v),v[0]);
    ARIADNE_TEST_EQUALS((x*x)(v),v[0]*v[0]);
    ARIADNE_TEST_EQUALS((x*x*x)(v),v[0]*v[0]*v[0]);
    ARIADNE_TEST_EQUALS((x*x*x*x)(v),v[0]*v[0]*v[0]*v[0]);
    ARIADNE_TEST_EQUALS(pow(x,3)(v),pow(v[0],3));
    ARIADNE_TEST_EQUALS(pow(x,4)(v),pow(v[0],4));
    ARIADNE_TEST_EQUALS((x*x*x*x-x*x)(v),v[0]*v[0]*v[0]*v[0]-v[0]*v[0]);
    ARIADNE_TEST_EQUALS((pow(x,4)-pow(x,2))(v),pow(v[0],4)-pow(v[0],2));

    ARIADNE_TEST_EQUALS(y(v),v[1]);
    ARIADNE_TEST_EQUALS((y*y)(v),v[1]*v[1]);
    ARIADNE_TEST_EQUALS((y*y*y)(v),v[1]*v[1]*v[1]);
    ARIADNE_TEST_EQUALS((y*y*y*y)(v),v[1]*v[1]*v[1]*v[1]);

    assert(v.size()==2);
    ARIADNE_TEST_EQUALS((x*y)(v),v[0]*v[1]);
    ARIADNE_TEST_EQUALS((x*x*y)(v),v[0]*v[0]*v[1]);
    ARIADNE_TEST_EQUALS((x*y*y)(v),v[0]*v[1]*v[1]);
*/
}



Int main() {
    MultiplePrecision mp(128);
    TestFourierPolynomial<FloatDPApproximation>(dp).test();
    TestFourierPolynomial<FloatMPApproximation>(mp).test();

    ThresholdSweeper<FloatDP> sweeper_dp(dp,1e-8);
    ThresholdSweeper<FloatMP> sweeper_mp(mp,std::pow(2.0,-64));
    ARIADNE_TEST_PRINT(sweeper_mp.precision());
    ARIADNE_TEST_PRINT(sweeper_mp);

    std::deque<int> deq;
    deq.push_back(0);
    deq.push_back(1);
    deq.push_back(2);
    deq.push_front(-1);
    std::cerr<<"deq[0]="<<deq[0]<<"\n";

    return ARIADNE_TEST_FAILURES;
}



