/***************************************************************************
 *            test_eigenvalues.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
#include <fstream>

#include "config.hpp"

#include "../test.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/eigenvalues.hpp"

using namespace std;
using namespace Ariadne;


template<class X> OutputStream& write_matrix(OutputStream& os, Matrix<X> const& A, SizeType wdth) {
    os << "\n";
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            os << (j==0?"    [":",")<<std::setw(wdth)<<A[i][j];
        }
        os << "]\n";
    }
    return os;
}
OutputStream& operator<<(OutputStream& os, Matrix<FloatDPApproximation> const& A) { return write_matrix(os,A,10); }
OutputStream& operator<<(OutputStream& os, Matrix<FloatMPApproximation> const& A) { return write_matrix(os,A,12); }
OutputStream& operator<<(OutputStream& os, Matrix<FloatDPBounds> const& A) { return write_matrix(os,A,10); }
OutputStream& operator<<(OutputStream& os, Matrix<FloatMPBounds> const& A) { return write_matrix(os,A,16); }

template<class X1, class X2> OutputStream& operator<<(OutputStream& os, Pair<Matrix<X1>,UpperTriangularMatrix<X2>> const& As) {
    ::operator<<(os,As.first); operator<<(os,As.second); return os;
}

namespace Ariadne {
typedef Matrix<FloatDPApproximation> FloatApproximationMatrix;
}


class TestEigenvalues {
    typedef MultiplePrecision PR;
    typedef FloatBounds<PR> X;
    PR pr;
  public:
    TestEigenvalues(PR prec) : pr(prec) { }
    Void test();
  private:
    Void test_misc();
    Void test_approximate();
    Void test_householder();
    Void test_upper_hessenberg();
};

Void
TestEigenvalues::test()
{
    ARIADNE_TEST_CALL(test_misc());
    ARIADNE_TEST_CALL(test_approximate());
    ARIADNE_TEST_CALL(test_householder());
    ARIADNE_TEST_CALL(test_upper_hessenberg());
}

Void
TestEigenvalues::test_misc()
{
    FloatDPApproximation::set_output_places(6);

    {
        Matrix<X> A({{3,1},{2,5}},pr);
        X d(5,6,pr);
//        Vector<X> v({{.3_dec,.4_dec},{0.9_dec,1.0_dec}},pr);
        Vector<X> v({{.3_dec,.4_dec,pr},{0.9_dec,1.0_dec,pr}});

        make_lpair(d,v)=eigenvalue_solve(A,d,v);
        std::cout << "d="<<d<<", v="<<v<<"\nAv-dv="<<A*v-d*v<<"\n";

    }

    // FIXME: Allow construction from Rational
    {
        ARIADNE_TEST_CONSTRUCT(Matrix<FloatMPBounds>,A,({{2,3,4},{3,5,6},{4,6,7}},pr));
        ARIADNE_TEST_PRINT(gram_schmidt(A));
    }


    {
        Matrix<X> A({{4,1,2,3},{1,8,4,5},{2,4,9,6},{3,5,6,10}},pr);
        auto HB = upper_hessenberg_factorisation(A);
        HouseholderProductMatrix<X> H=HB.first;
        UpperHessenbergMatrix<X> B=HB.second;
        std::cout << "A="<<A<<"\nB="<<B<<"\n";
        qr_step(B);
        std::cout << "B="<<B<<"\n";
        qr_step(B);
        std::cout << "B="<<B<<"\n";

    }


    // A0=Q0*R0; A1=R0*Q0; A1=Q1*R1; A2=R1*Q1
    // A2=Q1'*A1*Q1=Q1'*Q0'*A0*Q0*Q1;
    // (Q0*Q1)*A2 = A0*(Q0*Q1);
    {
        Matrix<X> A=Matrix<X>({{1,2,3,4},{5,6,7,8},{9,11,14,15},{10,13,0,12}},pr);
        Matrix<X> B=A;
        Matrix<X> V=Matrix<X>::identity(4,pr);
        Matrix<X> RR(4,4,pr); Matrix<X> QQ(4,4,pr);
        std::cout << "A=" << A << "\n";
        for(SizeType i=0; i!=5; ++i) {
            auto QR=gram_schmidt(B); OrthogonalMatrix<X> Q=std::get<0>(QR); UpperTriangularMatrix<X> R=std::get<1>(QR);
            B = R * Q;
            //std::cout << "Q=" << Q << "\nR=" << R << "\n";
            std::cout << "B=" << B << "\n";
            V = V * Q;
        }
        std::cout << "B=" << B << "\n";
        std::cout << "V=" << V << "\n";
        std::cout << "V'*V=" << transpose(V)*V << "\n";
        std::cout << "V*B-A*V="<< V*B-A*V << "\n";
    }

    {
        Matrix<X> A=Matrix<X>({{1,2,3,4},{5,6,7,8},{9,11,14,15},{10,13,0,12}},pr);
        auto QR=gram_schmidt(A); OrthogonalMatrix<X> Q=QR.first; auto R=std::get<1>(QR);
        std::cout << "A=" << A << "\nQ=" << Q << "\nR=" << R << "\n";
        std::cout << "Q'Q=" << transpose(Q)*Q << "\n";
        std::cout << "QR-A=" << Q*R-A << "\n";
        return;
    }

    {
        Matrix<X> A=Matrix<X>({{1,2,3,4},{5,6,7,8},{9,11,14,15},{10,13,0,12}},pr);
        auto HB = upper_hessenberg_factorisation(A);
        auto H=HB.first; auto B=HB.second;
    }

    {
        //FloatBounds<PR> e=FloatBounds<PR>(-1,+1,pr)/4;
        Matrix<X> A=Matrix<X>({{4,1,0},{1,3,1},{0,1,2}},pr);

        //FloatDPApproximation e=FloatDPApproximation(pr);
        //Vector<X> v({1.25_dyadic+e,0.5_dyadic+e,0+e}); X mu(4.5_dyadic+e);
        FloatBounds<PR> e=FloatBounds<PR>(-1,+1,pr)/16;
        Vector<X> v({0.75_dyadic+e,0.5625_dyadic+e,0.1875_dyadic+e}); X mu(4.75_dyadic+e);
        //v=[0.7887,0.5774,0.2113], mu=4.7320
//        Vector<X> v({0,1,0},pr); X mu(3);
//        Vector<X> v({0,0,1},pr); X mu(2);
        eigenvalue_solve(A,mu,v);
    }
    return;

    Vector<X> v({1,2,2,4},pr);
    SizeType n=v.size();
    Matrix<X> I=Matrix<X>::identity(n,pr);
    HouseholderMatrix<X> h(v);
    std::cout << "v="<<v<<"\n";
    std::cout << "h="<<Matrix<X>(h)<<"\n";
    std::cout << "h*h="<<Matrix<X>(h)*Matrix<X>(h)<<"\n";

    UpperHessenbergMatrix<X> B(n,pr);
    X x(pr); x=x+2; for(SizeType i=0; i!=n; ++i) { for(SizeType j=0; j!=n; ++j) {
        if(i<=j+1u) { B.set(i,j,x); x=x+1; } } }
    std::cout << "B="<<to_matrix(B)<<"\n";
    std::cout << "norm(h*h-I)="<<norm(Matrix<X>(h)*Matrix<X>(h)-I)<<"\n";
}


Void
TestEigenvalues::test_approximate()
{
    typedef FloatApproximation<PR> XA;
    ARIADNE_TEST_CONSTRUCT(Matrix<XA>,A,({{1,2,3,4},{5,6,7,8},{9,11,14,15},{10,13,0,12}},pr));
    Vector<XA> v({1,0,0,0},pr); XA mu(3,pr);
    tie(mu,v)=eigenvalue_solve(A,mu,v);
    ARIADNE_TEST_PRINT(mu);
    ARIADNE_TEST_PRINT(v);
    ARIADNE_TEST_PRINT(A*v-mu*v);

}


Void
TestEigenvalues::test_householder()
{
    Matrix<X> I=Matrix<X>::identity(4,pr);
    ARIADNE_TEST_CONSTRUCT(Matrix<X>,A,({{1,2,3,4},{5,6,7,8},{9,11,14,15},{10,13,0,12}},pr));
    Vector<X> v({2,3,5,7},pr);
    HouseholderMatrix<X> H(v);
    Matrix<X> B=H*A*H;
    ARIADNE_TEST_BINARY_PREDICATE(refines,A,(H*B*H));
    ARIADNE_TEST_BINARY_PREDICATE(refines,I,Matrix<X>(H*H));

}


Void
TestEigenvalues::test_upper_hessenberg()
{
    Matrix<X> I=Matrix<X>::identity(5,pr);
    ARIADNE_TEST_CONSTRUCT(Matrix<X>,A,({{1,2,3,4,5},{5,6,7,8,11},{9,11,14,15,17},{10,13,0,12,14},{15,13,0,17,19}},pr));
    auto HB = upper_hessenberg_factorisation(A);
    HouseholderProductMatrix<X> H=HB.first;
    OrthogonalMatrix<X> Q=H;
    UpperHessenbergMatrix<X> B=HB.second;
    ARIADNE_TEST_PRINT(H);
    ARIADNE_TEST_PRINT(Q);
    ARIADNE_TEST_PRINT(B);
    ARIADNE_TEST_PRINT(inverse(H));
    ARIADNE_TEST_BINARY_PREDICATE(refines,I,(Matrix<X>(H*inverse(H))));
    ARIADNE_TEST_BINARY_PREDICATE(refines,I,(Matrix<X>(inverse(H)*H)));
    ARIADNE_TEST_BINARY_PREDICATE(refines,A,(Matrix<X>(H*B*inverse(H))));
    ARIADNE_TEST_BINARY_PREDICATE(refines,A,(Matrix<X>(Q*B*transpose(Q))));

    ARIADNE_TEST_EXECUTE(qr_step(B));
};

Int main() {
    TestEigenvalues(MultiplePrecision(32)).test();
    return ARIADNE_TEST_FAILURES;
}

