/***************************************************************************
 *            solvers/nonlinear_programming.cpp
 *
 *  Copyright  2010-20  Pieter Collins
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

// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used

#include "function/functional.hpp"
#include "config.hpp"

#include <limits>

#include "utility/macros.hpp"
#include "conclog/logging.hpp"
#include "utility/tuple.hpp"
#include "utility/tribool.hpp"
#include "numeric/numeric.hpp"
#include "algebra/linear_algebra.decl.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/diagonal_matrix.hpp"
#include "algebra/symmetric_matrix.hpp"
#include "algebra/differential.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "function/procedure.hpp"

#include "solvers/nonlinear_programming.hpp"
#include "solvers/solver.hpp"
#include "algebra/multi_index-noaliasing.hpp"
#include "solvers/constraint_solver.hpp"

#include "algebra/expansion.inl.hpp"

using namespace ConcLog;

// CONCLOG_NUM_ARGS extracts the number of arguments of a variadic argument pack
#define _CONCLOG_NUM_ARGS(X8, X7, X6, X5, X4, X3, X2, X1, N, ...) N
#define CONCLOG_NUM_ARGS(...) _CONCLOG_NUM_ARGS(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1)

#define CONCLOG_EXPAND(x)  x
#define CONCLOG_FIRSTARG(x, ...)  (x)
#define CONCLOG_RESTARGS(x, ...)  (__VA_ARGS__)
#define CONCLOG_FOREACH(m, xs)  CONCLOG_FOREACH_(CONCLOG_NUM_ARGS xs, m, xs)
#define CONCLOG_FOREACH_(n, m, xs)  CONCLOG_FOREACH__(n, m, xs)
#define CONCLOG_FOREACH__(n, m, xs)  CONCLOG_FOREACH_##n(m, xs)
#define CONCLOG_FOREACH_1(m, xs)  m xs
#define CONCLOG_FOREACH_2(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_1(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_3(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_2(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_4(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_3(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_5(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_4(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_6(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_5(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_7(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_6(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_8(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_7(m, CONCLOG_RESTARGS xs)

#define CONCLOG_WRITE_VAR(x)    << #x << "=" << (x)
#define CONCLOG_WRITE_VARS(...) CONCLOG_FOREACH(CONCLOG_WRITE_VAR, (__VA_ARGS__))

#define CONCLOG_INITIALISE_AT(level) if (!Logger::instance().is_muted_at(level)) { std::ostringstream logger_stream; logger_stream << std::boolalpha
#define CONCLOG_FINALISE_AT(level) ; Logger::instance().println(level,logger_stream.str()); }

#define CONCLOG_PRINTLN_VARS_AT(level,...) { CONCLOG_INITIALISE_AT(level) CONCLOG_WRITE_VARS(__VA_ARGS__) CONCLOG_FINALISE_AT(level) }
#define CONCLOG_PRINTLN_VARS(...) { CONCLOG_PRINTLN_VARS_AT(0,__VA_ARGS__) }


namespace Ariadne {

namespace {

template<class FLT> using Exact = FLT;

template<class FLT> requires ARawFloat<FLT> inline Scalar<ExactNumber> to_generic(Scalar<Exact<FLT>> x) { return Scalar<ExactNumber>(x); }
template<class FLT> inline Scalar<ValidatedNumber> to_generic(Scalar<Bounds<FLT>> x) { return Scalar<ValidatedNumber>(x); }
template<class FLT> inline Scalar<ApproximateNumber> to_generic(Scalar<Approximation<FLT>> x) { return Scalar<ApproximateNumber>(x); }

template<class FLT> inline Vector<ExactNumber> to_generic(Vector<Exact<FLT>> v) { return Vector<ExactNumber>(v); }
template<class FLT> inline Vector<ValidatedNumber> to_generic(Vector<Bounds<FLT>> v) { return Vector<ValidatedNumber>(v); }
template<class FLT> inline Vector<ApproximateNumber> to_generic(Vector<Approximation<FLT>> v) { return Vector<ApproximateNumber>(v); }

inline Sweeper<FloatDP> default_sweeper() { return Sweeper<FloatDP>(); }

typedef DiagonalMatrix<FloatDPBounds> FloatDPBoundsDiagonalMatrix;

typedef Vector<FloatDPApproximation> ApproximateVectorType;
typedef VectorRange<ApproximateVectorType> ApproximateVectorTypeRange;
typedef Matrix<FloatDPApproximation> FloatDPApproximationMatrix;
typedef DiagonalMatrix<FloatDPApproximation> FloatDPApproximationDiagonalMatrix;
typedef Differential<FloatDPApproximation> FloatDPApproximationDifferential;

typedef Vector<FloatDPBounds> FloatDPBoundsVector;
typedef Matrix<FloatDPBounds> FloatDPBoundsMatrix;
typedef Differential<FloatDPBounds> FloatDPBoundsDifferential;
typedef Vector<FloatDP> ExactFloatDPVectorType;

typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Matrix<UpperIntervalType> UpperIntervalMatrixType;

typedef FloatDPApproximation ApproximateNumericType;

template<class T> struct Pretty;

template<> struct Pretty<Matrix<FloatDPApproximation>> {
    static const uint prec=4; static const uint wdth=prec+6;
    Matrix<FloatDPApproximation> const& A;
    friend OutputStream& operator<<(OutputStream& os, Pretty<Matrix<FloatDPApproximation>> const& pA) {
        auto& A=pA.A;
        os << "\n";
        for(SizeType i=0; i!=A.row_size(); ++i) {
            os << "[";
            for(SizeType j=0; j!=A.column_size(); ++j) {
                if (j!=0) { os << ","; } os << std::setprecision(prec) << std::setw(wdth) << A[i][j].get_d();
            }
            os << "]\n";
        }
        return os;
    }
};
template<class T> Pretty<T> pretty(T const& t) { return Pretty<T>{t};}

class NegatedIdentityMatrix;
class IdentityMatrix {
    SizeType _n;
  public:
    IdentityMatrix(SizeType n) : _n(n) { }
    SizeType size() const { return this->_n; }
    Int at(SizeType i, SizeType j) const { return (i==j); }
    friend NegatedIdentityMatrix operator-(IdentityMatrix);
};
class NegatedIdentityMatrix {
    SizeType _n;
  public:
    NegatedIdentityMatrix(SizeType n) : _n(n) { }
    SizeType size() const { return this->_n; }
    Int at(SizeType i, SizeType j) const { return -(i==j); }
    friend NegatedIdentityMatrix operator-(IdentityMatrix I) { return NegatedIdentityMatrix(I.size()); }
};


inline Vector<Differential<RawFloatDP>>const& cast_raw(Vector<Differential<FloatDPApproximation>>const& v) {
    return reinterpret_cast<Vector<Differential<RawFloatDP>>const&>(v);
}

inline Vector<FloatDPApproximation>& cast_approximate(Vector<RawFloatDP>& v) {
    return reinterpret_cast<Vector<FloatDPApproximation>&>(v);
}

inline Vector<FloatDPApproximation>& cast_approximate(Vector<FloatDPApproximation>& v) {
    return v;
}
inline Vector<Differential<FloatDPApproximation>>const& cast_approximate(Vector<Differential<RawFloatDP>>const& v) {
    return reinterpret_cast<Vector<Differential<FloatDPApproximation>>const&>(v);
}


inline FloatDPApproximation operator*(ApproximateDouble x1, FloatDPApproximation x2) {
    return FloatDPApproximation(x1,x2.precision()) * x2;
}
inline FloatDPApproximation operator*(FloatDPApproximation x1, ApproximateDouble x2) {
    return x1 * FloatDPApproximation(x2,x1.precision());
}
inline ApproximateKleenean operator<(FloatDPApproximation x1, ApproximateDouble x2) {
    return x1 < FloatDPApproximation(x2,x1.precision());
}


inline UpperIntervalType dot(Vector<UpperIntervalType> const& bx1, Vector<ExactIntervalType> const& bx2) {
    return dot(bx1,Vector<UpperIntervalType>(bx2));
}

template<class X1, class IVL2> inline decltype(auto) dot(Vector<X1> const& x1, Box<IVL2> const& bx2) {
    return dot(x1,cast_vector(bx2)); }

template<class X> inline
DiagonalMatrix<X> diagonal_matrix(const Vector<X>& v) {
    return DiagonalMatrix<X>(v.array());
}

template<class X> inline
Bool epos(const Vector<X>& x) {
    for(SizeType i=0; i!=x.size(); ++i) { if(x[i]<=0) { return false; } } return true;
}

template<class X> inline
Bool eneg(const Vector<X>& x) {
    for(SizeType i=0; i!=x.size(); ++i) { if(x[i]>=0) { return false; } } return true;
}

template<class X, class XX> inline
Bool egtr(const Vector<X>& x, const XX& s) {
    for(SizeType i=0; i!=x.size(); ++i) { if(decide(x[i]<=s)) { return false; } } return true;
}

template<class X, class XX> inline
Bool elss(const Vector<X>& x, const XX& s) {
    for(SizeType i=0; i!=x.size(); ++i) { if(decide(x[i]>=s)) { return false; } } return true;
}

template<class X> inline
Vector<X> eadd(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]+y[i]; } return r;
}

template<class X> inline
Vector<X> eadd(const Vector<X>& x, const X& s) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]+s; } return r;
}

template<class X> inline
Vector<X> eadd(const X& s, const Vector<X>& x) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=s+x[i]; } return r;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const X& s) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]-s; } return r;
}

template<class X> inline
Vector<X> esub(const X& s, const Vector<X>& x) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=s-x[i]; } return r;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]-y[i]; } return r;
}

template<class X1, class X2, class X3> using Product3Type = decltype(declval<X1>()*declval<X2>()*declval<X3>());

template<class X1, class X2> inline
Vector<ProductType<X1,X2>> emul(const Vector<X1>& x1, const Vector<X2>& x2) {
    Vector<ProductType<X1,X2>> r(x1.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x1[i]*x2[i]; } return r;
}

template<class X1, class X2, class X3> inline
Vector<Product3Type<X1,X2,X3>> emul(const Vector<X1>& x1, const Vector<X2>& x2, const Vector<X3>& x3) {
    Vector<Product3Type<X1,X2,X3>> r(x1.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x1[i]*x2[i]*x3[i]; } return r;
}

template<class X1, class X2, class X3> requires AScalar<X1> inline
Vector<Product3Type<X1,X2,X3>> emul(const Scalar<X1>& s1, const Vector<X2> x2, const Vector<X3>& x3) {
    Vector<Product3Type<X1,X2,X3>> r(x2.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=s1*x2[i]*x3[i]; } return r;
}

template<class X, class XX> inline
Vector<X> ediv(const Vector<X>& x, const Vector<XX>& z) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]/z[i]; } return r;
}

template<class X> inline
Vector<X> ediv(const X& s, const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=s/z[i]; } return r;
}

template<class X> inline
Vector<X> erec(const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=rec(z[i]); } return r;
}

template<class X> inline
Vector<X> esqr(const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=sqr(z[i]); } return r;
}

inline
ExactIntervalType eivl(const RawFloatDPVector& x) {
    ARIADNE_ASSERT(x.size()>0); ExactIntervalType r=ExactIntervalType(FloatDP(x[0]));
    for(SizeType i=1; i!=x.size(); ++i) { r=hull(r,FloatDP(x[i])); } return r;
}

} // namespace

Matrix<ApproximateNumericType> join(Matrix<ApproximateNumericType> const&, Matrix<ApproximateNumericType> const&, Matrix<ApproximateNumericType> const&);

Matrix<ApproximateNumericType> join(Matrix<ApproximateNumericType> const& A1, Matrix<ApproximateNumericType> const& A2, Matrix<ApproximateNumericType> const& A3) {
    SizeType m=A1.row_size(); SizeType n1=A1.column_size(); SizeType n2=A2.column_size(); SizeType n3=A3.column_size();
    Matrix<ApproximateNumericType> A123(m,n1+n2+n3,(A1.zero_element()+A2.zero_element()+A3.zero_element()));
    project(A123,range(0,m),range(0,n1))=A1;
    project(A123,range(0,m),range(n1,n1+n2))=A2;
    project(A123,range(0,m),range(n1+n2,n1+n2+n3))=A3;
    return A123;
}

template<class X> Matrix<X> cojoin(Matrix<X> const& A1, Matrix<X> const& A2, Matrix<X> const& A3) {
    SizeType n=A1.column_size(); SizeType m1=A1.row_size(); SizeType m2=A2.row_size(); SizeType m3=A3.row_size();
    Matrix<X> A123(m1+m2+m3,n,(A1.zero_element()+A2.zero_element()+A3.zero_element()));
    project(A123,range(0,m1),range(0,n))=A1;
    project(A123,range(m1,m1+m2),range(0,n))=A2;
    project(A123,range(m1+m2,m1+m2+m3),range(0,n))=A3;
    return A123;
}


// Compute S+=ADA^T, where D is diagonal and S is symmetric.
template<class X>
Void adat(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            X ADij=A[i1][j]*D[j];
            for(SizeType i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(SizeType i1=1; i1!=m; ++i1) {
        for(SizeType i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
Void atda(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    ARIADNE_PRECONDITION(S.row_size()==S.column_size());
    ARIADNE_PRECONDITION(S.column_size()==A.column_size());
    ARIADNE_PRECONDITION(D.size()==A.row_size());

    const SizeType m=A.column_size();
    const SizeType n=A.row_size();
    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            X ATDij=A[j][i1]*D[j];
            for(SizeType i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ATDij*A[j][i2];
            }
        }
    }
    for(SizeType i1=1; i1!=m; ++i1) {
        for(SizeType i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
Void atda(Matrix<X>& S, const Matrix<X>& A, const DiagonalMatrix<X>& D)
{
    atda(S,A,D.diagonal());
}

// Compute S=ADA^T, where D is diagonal.
template<class X>
Matrix<X> adat(const Matrix<X>& A, const Vector<X>& D)
{
    const SizeType m=A.row_size();
    Matrix<X> S=Matrix<X>::zero(m,m);
    adat(S,A,D);
    return S;
}

// Compute S+=AA^T
template<class X>
Matrix<X> amulat(const Matrix<X>& A)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    Matrix<X> S(m,m);
    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            for(SizeType i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=A[i1][j]*A[i2][j];
            }
        }
    }
    for(SizeType i1=1; i1!=m; ++i1) {
        for(SizeType i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
    return S;
}

template<class X> inline Bool all_greater(const Vector<X>& x, const X& e) {
    for(SizeType i=0; i!=x.size(); ++i) { if(x[i]<=e) { return false; } } return true;
}






template<class X> Vector< Differential<X> > second_derivative(const ValidatedVectorMultivariateFunction& f, const Vector<X>& x) {
    Vector< Differential<X> > d=Differential<X>::variables(f.result_size(),f.argument_size(),2);
    return f.evaluate(d);
}

template<class Vec, class Diff> Void set_gradient(Vec& g, const Diff& D) {
    SizeType i=0;
    typename Diff::ConstIterator iter=D.begin();
    if(iter!=D.end() && iter->index().degree()==0) { ++iter; }
    while(iter!=D.end() && iter->index().degree()<=2) {
        while(iter->index()[i]==0) { ++i; }
        g[i]=iter->coefficient();
        ++iter;
    }
}

template<class Mx, class Diff> Void set_jacobian_transpose(Mx& A, const Vector<Diff>& D) {
    for(SizeType j=0; j!=A.column_size(); ++j) {
        for(SizeType i=0; i!=A.row_size(); ++i) {
            A[i][j]=D[j][i];
        }
    }
}

template<class Mx, class Diff> Void set_hessian(Mx& H, const Diff& D) {
    typedef typename Diff::ValueType X;
    SizeType i=0; SizeType j=1;
    typename Diff::ConstIterator iter=D.begin();
    while(iter!=D.end() && iter->index().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->index().degree()<=2) {
        UniformConstReference<MultiIndex> a=iter->index();
        UniformConstReference<X> c=iter->coefficient();
        while(a[i]==0) { ++i; j=i+1; }
        if(a[i]==2) { H[i][i]=c; }
        else { while(a[j]==0) { ++j; } H[i][j]=c; H[j][i]=c; }
        ++iter;
    }
}

template<class Mx, class S, class Diff> Void add_hessian(Mx& H, const S& s, const Diff& D) {
    typedef typename Diff::ValueType X;
    typename Diff::ConstIterator iter=D.begin();
    while(iter!=D.end() && iter->index().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->index().degree()==2) {
        UniformConstReference<MultiIndex> a=iter->index();
        UniformConstReference<X> c=iter->coefficient();
        SizeType i=0;
        while(a[i]==0) { ++i; }
        if(a[i]==2) { H[i][i]+=s*c; }
        else { SizeType j=i+1; while(a[j]==0) { ++j; } H[i][j]+=s*c; H[j][i]+=s*c; }
        ++iter;
    }
}

// Compute the product (A -A I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_mul(const Matrix<XX>& A, const Vector<X>& v)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    ARIADNE_ASSERT(v.size()==2*(m+n));
    Vector<X> r(m+1u,v.zero_element());
    for(SizeType i=0; i!=m; ++i) {
        r[i]=v[2*n+i]-v[2*n+m+i];
        for(SizeType j=0; j!=n; ++j) {
            r[i]+=A[i][j]*(v[j]-v[n+j]);
        }
    }
    for(SizeType k=0; k!=2*(m+n); ++k) {
        r[m]+=v[k];
    }
    return r;
}

// Compute the product (AT 1 \\ -AT 1 \\ I 1 \\ -I 1) I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_trmul(const Matrix<XX>& A, const Vector<X>& w)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    ARIADNE_ASSERT(w.size()==m+1);
    Vector<X> r(2*(m+n),w.zero_element());
    for(SizeType j=0; j!=n; ++j) {
        r[j]=0;
        for(SizeType i=0; i!=m; ++i) {
            r[j]+=A[i][j]*w[i];
        }
        r[n+j]=-r[j];
        r[j]+=w[m];
        r[n+j]+=w[m];
    }
    for(SizeType i=0; i!=m; ++i) {
        r[2*n+i]=w[i]+w[m];
        r[2*n+m+i]=-w[i]+w[m];
    }
    return r;
}


// Compute the product \f$\hat{A}^T \hat{D} \hat{A} + \hat{H}\f$ where \f$\hat{A}=\left(\begin{matrix}A&-A&I&-I\\1&1&1&1\end{matrix}\right)\f$ and \f$\hat{D}=D\f$ is diagonal.
template<class X> Matrix<X> feasibility_adat(const Matrix<X>& H, const Matrix<X>& A, const Vector<X>& D)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    ARIADNE_ASSERT(H.row_size()==m);
    ARIADNE_ASSERT(H.column_size()==m);
    ARIADNE_ASSERT(D.size()==2*(m+n));
    Matrix<X> S(m+1,m+1,H.zero_element());

    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=m; ++j) { S[i][j] = H[i][j]; } }
    for(SizeType i=0; i!=m; ++i) { S[i][m]=0; S[m][i]=0; } S[m][m]=0;

    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            X ADij=A[i1][j]*(D[j]+D[n+j]);
            for(SizeType i2=0; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(SizeType i=0; i!=m; ++i) {
        S[i][i]+=(D[2*n+i]+D[2*n+m+i]);
    }
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            S[i][m]+=A[i][j]*(D[j]-D[n+j]);
        }
        S[i][m]+=(D[2*n+i]-D[2*n+m+i]);
        S[m][i]=S[i][m];
    }
    for(SizeType k=0; k!=2*(m+n); ++k) {
        S[m][m]+=D[k];
    }

    return S;
}




template<class P> OutputStream& operator<<(OutputStream& os, FeasibilityProblem<P> const& p) {
    return os << "FeasibilityProblem( D=" << p.D << ", g=" << p.g << ", C=" << p.C << " )";
}
template OutputStream& operator<<(OutputStream&, FeasibilityProblem<ApproximateTag> const&);
template OutputStream& operator<<(OutputStream&, FeasibilityProblem<ValidatedTag> const&);

template<class P> OutputStream& operator<<(OutputStream& os, OptimisationProblem<P> const& p) {
    return os << "OptimisationProblem( f=" << p.f << ", D=" << p.D << ", g=" << p.g << ", C=" << p.C << " )";
}
template OutputStream& operator<<(OutputStream&, OptimisationProblem<ApproximateTag> const&);
template OutputStream& operator<<(OutputStream&, OptimisationProblem<ValidatedTag> const&);


template<class R>
class ConstrainedFeasibilityMatrix {
    ConstrainedFeasibilityMatrix(const Vector<R>& x, const Vector<R>& z, const Matrix<R>& a, const Matrix<R>& h)
        : X(x), Z(z), D(ediv(x,z)), A(a), H(h) { }

    template<class RR> Tuple< Vector<RR>,Vector<RR>,Vector<RR> >
    mul(const Vector<RR>& x, const Vector<RR>& yt, const Vector<RR>& z) const {
        Vector<RR> nx=Z*x+X*x;
        Vector<RR> nyt=H*yt-A*x;
        Vector<RR> nz=H*yt-A*x;
        return make_tuple(nx,nyt,nz);
    }

    template<class RR> Tuple< Vector<RR>,Vector<RR>,Vector<RR> >
    solve(const Vector<RR>& x, const Vector<RR>& yt, const Vector<RR>& z) const {
        Sinv=inverse(feasibility_adat(H,A,D));
        Vector<RR> rx=Z.solve(x);
        Vector<RR> ryt=yt+feasibility_mul(A,rx-D*z);
        Vector<RR> rz=z;
        ryt=Sinv*ryt;
        rz=rz-feasibility_trmul(A,ryt);
        rx=rz-D*rz;
        return make_tuple(rx,ryt,rz);
    }

    DiagonalMatrix<R> X;
    DiagonalMatrix<R> Z;
    DiagonalMatrix<R> D;
    const Matrix<R>& A;
    const Matrix<R>& H;
    Matrix<R> Sinv;
};


inline Box<Interval<FloatDP>> cast_exact_widen(Box<Interval<FloatDP>> const& bx, FloatDP e) {
    Box<Interval<FloatDP>> r(bx);
    for(SizeType i=0; i!=bx.size(); ++i) {
        r[i]=Interval<FloatDP>(sub(down,bx[i].lower_bound(),e),add(up,bx[i].upper_bound(),e));
    }
    return r;
}



//------- FeasibilityChecker -----------------------------------//

FeasibilityChecker* FeasibilityChecker::
clone() const
{
    return new FeasibilityChecker(*this);
}


ApproximateKleenean FeasibilityChecker::
almost_feasible_point(ValidatedFeasibilityProblem p, ApproximateVector ax, ApproximateNumber e) const
{
    Vector<Approximation<FLT>> axf(ax,pr);
    Approximation<FLT> ef(e,pr);

    if (!probably(contains(p.D,axf))) { return false; }
    Vector<Approximation<FLT>> gxf=p.g(axf);
    return contains(cast_exact_widen(p.C,cast_exact(ef)),gxf);
}


ValidatedKleenean FeasibilityChecker::
is_feasible_point(ValidatedFeasibilityProblem p, ExactVector x) const
{
    Vector<Bounds<FLT>> xf(x,pr);

    if (not possibly(contains(p.D,xf))) { return false; }
    Vector<Bounds<FLT>> gxf=p.g(xf);
    return contains(p.C,gxf);
}

ValidatedKleenean FeasibilityChecker::
contains_feasible_point(ValidatedFeasibilityProblem p, UpperBoxType X) const
{
    if (this->validate_feasibility(p,cast_singleton(X))) {
        return true;
    } else {
        p.D=intersection(p.D,cast_exact_box(X));
        if (this->validate_infeasibility(p)) {
            return false;
        } else {
            return indeterminate;
        }
    }
}


ValidatedKleenean FeasibilityChecker::
check_feasibility(ValidatedFeasibilityProblem p,
                  ValidatedVector x, ExactVector y) const
{
    return this->check_feasibility(p, Vector<Bounds<FLT>>(x,dp), Vector<Bounds<FLT>>(y,pr));
}

ValidatedKleenean FeasibilityChecker::
check_feasibility(ValidatedFeasibilityProblem p,
                  Vector<Bounds<FLT>> x, Vector<Bounds<FLT>> y) const
{
    if (this->validate_feasibility(p,x)) {
        return true;
    } else if (this->validate_infeasibility(p,x,y)) {
        return false;
    } else {
        return indeterminate;
    }
}


Bool FeasibilityChecker::
validate_feasibility(ValidatedVectorMultivariateFunction h,
                     ValidatedVector x) const
{
    return this->validate_feasibility(h, Vector<Bounds<FLT>>(x,pr));
}


Bool FeasibilityChecker::
validate_feasibility(ValidatedVectorMultivariateFunction h,
                     Vector<Bounds<FLT>> X) const
{
    // Let h:R^n->R^m with m<n, and X define a box in R^n.
    // Attempt to solve h(c+R*B*w)=0 for w, where c=mid(X), R is a diagonal scaling matrix, and B is an n*m matrix.
    // Let z lie in the unit box Z=[-1:+1]^n, and take x=R*z+c where R=diag(rad(X)) is a scaling matrix.
    // Then the function z->h(c+R*z) has Jacobian Dh(X)*R = [A]*R over Z.
    // Then for an interval Newton step, we have [A]*R*B dw = -h(c+B*w)
    // Take A to be an approximation to [A], such as A=mid([A]) or A=Dh(C), and set B=(A*R)^T=R*A^T
    // Then we solve h(c+R*R*AT*w)=0 using Newton's method centred at w=0, yielding W' = -([A]*R*R*A^T)\h(c)

    Vector<Approximation<FLT>> ca = midpoint(X);
    Vector<FLT> c = cast_exact(ca);
    Matrix<FLT> AT = transpose(cast_exact(h.jacobian(ca)));
    CONCLOG_PRINTLN("A="<<transpose(AT));


    Vector<Bounds<FLT>> W=h(X);
    Matrix<Bounds<FLT>> A = h.jacobian(X);
    DiagonalMatrix<FloatDP> R(Array<FloatDP>(X.size(), [&X](SizeType i){return cast_exact(X[i].error());}));


    // Could take w to be centre of set W
    // Vector<FloatDP> w=cast_exact(W);
    // new_W = w - gs_solve(A*AT,h(x+AT*w));
    // Easier to take w = 0, which should be an element of W
    Vector<Bounds<FLT>> new_W = - gs_solve(A*(R*R*AT),h(c));
    Vector<Bounds<FLT>> new_X = c + (R * R) * (AT * new_W);

    if (refines(new_X,X)) {
        return true;
    } else {
        if (refines(new_W,W)) {
            ARIADNE_WARN("Did not verify feasible point in "<<X<<", but one may exist.");
        }
        return false;
    }
}

Bool FeasibilityChecker::
validate_feasibility(ValidatedFeasibilityProblem p,
                     ValidatedVector x) const
{
    return this->validate_feasibility(p,Vector<Bounds<FLT>>(x,pr));
}

Bool FeasibilityChecker::
validate_feasibility(ValidatedFeasibilityProblem p,
                     Vector<Bounds<FLT>> x) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    CONCLOG_PRINTLN("x="<<x);

    ARIADNE_PRECONDITION(x.size()==D.size());


    for(SizeType i=0; i!=D.size(); ++i) {
        CONCLOG_PRINTLN_AT(1,"x["<<i<<"]="<<x[i]<<", D["<<i<<"]="<<D[i]);
        if (!possibly(contains(D[i],x[i]))) {
            return false;
        } else if (definitely(contains(D[i],x[i]))) {
        } else {
            x[i]=cast_singleton(intersection(UpperIntervalType(x[i]),D[i]));
        }
    }

    Vector<Bounds<FLT>> w=g(x);
    CONCLOG_PRINTLN_AT(1,"w=g(x)="<<w);

    List<SizeType> equalities;
    for(SizeType j=0; j!=C.size(); ++j) {
        CONCLOG_PRINTLN_AT(1,"w["<<j<<"]="<<w[j]<<", C["<<j<<"]="<<C[j]);
        if (!possibly(contains(C[j],w[j]))) {
            return false;
        } else if (definitely(contains(C[j],w[j]))) {
        } else {
            // NOTE: It is safe to try solving a constraint as an equality
            if ( true || decide(C[j].lower_bound()==C[j].upper_bound()) ) {
                equalities.append(j);
            }
        }
    }
    if(equalities.empty()) {
        return true;
    }

    ValidatedVectorMultivariateFunction h(equalities.size(),g.domain());
    for(SizeType i=0; i!=equalities.size(); ++i) {
        SizeType j=equalities[i];
        h[i] = g[j]-static_cast<ExactNumber>(intersection(UpperIntervalType(w[j]),C[j]).midpoint());
    }
    CONCLOG_PRINTLN_AT(1,"h="<<h);

    return this->validate_feasibility(h,x);
}


Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       UpperBoxType B, ExactVector y) const
{
    return validate_infeasibility(p,B,Vector<Bounds<FLT>>(y,pr));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       UpperBoxType B, Vector<Exact<FLT>> y) const
{
    return validate_infeasibility(p,B,Vector<Bounds<FLT>>(y));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       UpperBoxType B, Vector<Bounds<FLT>> y) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;
    ExactBoxType DB = intersection(D,cast_exact_box(B));
    Vector<Approximation<FLT>> xa = midpoint(B);
    return this->validate_infeasibility(ValidatedFeasibilityProblem(DB,g,C),xa,y);
}


Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       ApproximateVector xa, ExactVector y) const
{
    return validate_infeasibility(p, Vector<Approximation<FLT>>(xa,dp), Vector<Bounds<FLT>>(y,pr));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       Vector<Approximation<FLT>> xa, Vector<Exact<FLT>> y) const
{
    return validate_infeasibility(p, xa, Vector<Bounds<FLT>>(y));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       Vector<Approximation<FLT>> xa, Vector<Bounds<FLT>> y) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", xa="<<xa<<", y="<<y);

    ARIADNE_PRECONDITION(xa.size()==D.size());
    ARIADNE_PRECONDITION(y.size()==C.size());

    if(y.size()==0) { return D.is_empty(); }

    UpperIntervalType yC = dot(y,C);

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
        for(SizeType j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j];
    }
    UpperIntervalType ygD = apply(tyg.function(),D);
    // UpperIntervalType ygD = dot(y,apply(g,D));

    if(definitely(disjoint(yC,ygD))) {
        CONCLOG_PRINTLN("Infeasible");
        return true;
    }

    UpperIntervalMatrixType dgD = jacobian_range(g,cast_vector(D));
    UpperIntervalVectorType ydgD = transpose(dgD)*y;

    Vector<Bounds<FLT>> x = cast_exact(xa);
    Bounds<FLT> ygx = tyg(x);
    // ValidatedNumericType ygx = dot(y,g(x));
    UpperIntervalType ygDx = UpperIntervalType(ygx);
    for(SizeType i=0; i!=x.size(); ++i) {
        ygDx += ydgD[i] * (D[i]-x[i]);
    }

    CONCLOG_PRINTLN("yC="<<yC<<", ygD="<<ygD<<", ygx="<<ygx<<", ydgD="<<ydgD<<", ygDx="<<ygDx);

    if(definitely(disjoint(yC,intersection(ygD,ygDx)))) {
        CONCLOG_PRINTLN("Infeasible"); return true; }
    else { return false; }
}


Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p, ExactVector y) const
{
    return this->validate_infeasibility(p,Vector<Bounds<FLT>>(y,pr));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p, Vector<Exact<FLT>> y) const
{
    return this->validate_infeasibility(p,Vector<Bounds<FLT>>(y));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p, Vector<Bounds<FLT>> y) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", y="<<y);

    ARIADNE_PRECONDITION(y.size()==C.size());

    if(y.size()==0) { return D.is_empty(); }

    UpperIntervalType yC = dot(y,C);

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(SizeType j=0; j!=y.size(); ++j) {
        tyg+=y[j]*tg[j];
    }
    UpperIntervalType ygD = apply(tyg.function(),D);
    // UpperIntervalType ygD = dot(y,apply(g,D));

    if(definitely(disjoint(yC,ygD))) {
        CONCLOG_PRINTLN("Infeasible");
        return true;
    } else {
        return false;
    }
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);

    if (D.is_empty()) {
        return true;
    }

    UpperBoxType W=apply(g,D);
    return (definitely(disjoint(W,C)));
}


OutputStream& operator<<(OutputStream& os, FeasibilityChecker const& fc) {
    return os << "FeasibilityChecker()";
}



//------- OptimiserBase -----------------------------------//

const OptimiserBase::FLT OptimiserBase::zero = FLT(0,pr);
const OptimiserBase::FLT OptimiserBase::one = FLT(1,pr);


auto OptimiserBase::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;
    return this->minimise(ValidatedOptimisationProblem{f,D,g,C});
}

auto OptimiserBase::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->minimise(f,D,gh,C);
}

auto OptimiserBase::
minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const -> ApproximateVector
{
    CONCLOG_SCOPE_CREATE;
    return this->minimise(ApproximateOptimisationProblem{f,D,g,C});
}

auto OptimiserBase::
minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateVectorMultivariateFunction h) const -> ApproximateVector
{
    CONCLOG_SCOPE_CREATE;
    ApproximateVectorMultivariateFunction gh=join(g,h);
    ApproximateBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->minimise(f,D,gh,C);
}



auto OptimiserBase::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;
    return this->feasible(ValidatedFeasibilityProblem{D,g,C});
}

auto OptimiserBase::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->feasible(D,gh,C);
}


//------- PenaltyFunctionOptimiser ------------------------------------------//

PenaltyFunctionOptimiser* PenaltyFunctionOptimiser::
clone() const
{
    return new PenaltyFunctionOptimiser(*this);
}

auto PenaltyFunctionOptimiser::
minimise(ValidatedOptimisationProblem f) const -> ValidatedVector
{
    ARIADNE_NOT_IMPLEMENTED;
}

auto PenaltyFunctionOptimiser::
minimise(ApproximateOptimisationProblem p) const ->ApproximateVector
{
    ARIADNE_NOT_IMPLEMENTED;
}

auto PenaltyFunctionOptimiser::
feasible(ValidatedFeasibilityProblem p) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<p.D<<" g="<<p.g<<" C="<<p.C<<" ");

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ApproximateVectorType x=midpoint(D);

    ApproximateVectorType w=midpoint(C);
    for(SizeType i=0; i!=C.size(); ++i) {
        if(C[i].upper_bound()==+infty) { w[i]=C[i].lower_bound()+one; }
        else if(C[i].lower_bound()==-infty) { w[i]=C[i].upper_bound()-one; }
    }

    ApproximateVectorType y(C.size(),zero);

    CONCLOG_PRINTLN("x="<<x<<" w="<<w<<" y="<<y);

    for(SizeType i=0; i!=10; ++i) {
        this->feasibility_step(ApproximateFeasibilityProblem(D,g,C),w,x,y);
    }
    return FeasibilityChecker().check_feasibility(p,cast_exact(x),cast_exact(y));
}

auto PenaltyFunctionOptimiser::
feasible(ApproximateFeasibilityProblem p) const -> ApproximateKleenean
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<p.D<<" g="<<p.g<<" C="<<p.C<<" ");

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ApproximateVectorType x=midpoint(D);

    ApproximateVectorType w=midpoint(C);
    for(SizeType i=0; i!=C.size(); ++i) {
        if (decide(C[i].upper_bound()==+infty)) { w[i]=C[i].lower_bound()+one; }
        else if (decide(C[i].lower_bound()==-infty)) { w[i]=C[i].upper_bound()-one; }
    }

    ApproximateVectorType y(C.size(),zero);

    CONCLOG_PRINTLN("x="<<x<<" w="<<w<<" y="<<y);

    for(SizeType i=0; i!=10; ++i) {
        this->feasibility_step(p,w,x,y);
    }
    return C.contains(g(x));
}

Void PenaltyFunctionOptimiser::
feasibility_step(ApproximateFeasibilityProblem p,
                 ApproximateVectorType& w, ApproximateVectorType& x, ApproximateNumberType& mu) const
{
    CONCLOG_SCOPE_CREATE;
    auto& d=p.D; auto& g=p.g; auto& c=p.C;

    ApproximateVectorMultivariateFunction h(0u,d.dimension());
    const SizeType n=d.size();
    const SizeType m=c.size();
    const SizeType l=h.result_size();

    CONCLOG_PRINTLN("x="<<x);
    CONCLOG_PRINTLN("w="<<w);

    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));

    mu *= 0.5;
    CONCLOG_PRINTLN("mu="<<mu);

    // G is the constraint value vector
    ApproximateVectorType gx = ddgx.value();
    ApproximateVectorType hx = ddhx.value();
    CONCLOG_PRINTLN("g(x)="<<gx);
    CONCLOG_PRINTLN("h(x)="<<hx);

    // A is the transpose derivative matrix aij=dgi/dxj
    FloatDPApproximationMatrix A = transpose(ddgx.jacobian());
    CONCLOG_PRINTLN("A=Dg(x)="<<A);
    FloatDPApproximationMatrix B = transpose(ddhx.jacobian());
    // FIXME: Due to problems with zero-element differential, need to resize matrix if no h
    if(l==0) { B.resize(n,0); }
    CONCLOG_PRINTLN("B=Dh(x)="<<B);

    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j] lambda[j]*dg[j]/dx[i1]dx[i2]
    FloatDPApproximationMatrix H(n,n,dp);
    for(SizeType j=0; j!=m; ++j) { H += (gx[j]-w[j]) * ddgx[j].hessian(); }
    for(SizeType k=0; k!=l; ++k) { H += (hx[k]) * ddhx[k].hessian(); }
    CONCLOG_PRINTLN("H="<<H);

    FloatDPApproximationDiagonalMatrix E(n,dp);
    FloatDPApproximationDiagonalMatrix D(m,dp);
    for(SizeType i=0; i!=n; ++i) { E[i] = rec(sqr(x[i]-d[i].lower_bound())) + rec(sqr(d[i].upper_bound()-x[i])); }
    for(SizeType j=0; j!=m; ++j) { D[j] = rec(sqr(w[j]-c[j].lower_bound())) + rec(sqr(c[j].upper_bound()-w[j])); }
    CONCLOG_PRINTLN("E="<<E);
    CONCLOG_PRINTLN("D="<<D);

    FloatDPApproximationMatrix S = H + B * transpose(B);
    S += E;
    CONCLOG_PRINTLN("S="<<S);

    FloatDPApproximationMatrix R=inverse(S);
    CONCLOG_PRINTLN("inverse(S)="<<R);

    // Compute residuals
    ApproximateVectorType rx = A*gx + B * hx ; // + 1/(x.upper_bound()-x) + 1/x.lower_bound()-x if no regularisation
    ApproximateVectorType rw = w-gx;

    CONCLOG_PRINTLN("rx="<<rx);
    CONCLOG_PRINTLN("rw="<<rw);

    ApproximateVectorType dx = R * (rx + A * rw);
    ApproximateVectorType dw = rw + transpose(A)*dx;
    CONCLOG_PRINTLN("dx="<<dx);
    CONCLOG_PRINTLN("dw="<<dw);


    ApproximateVectorType newx(n,dp);
    ApproximateVectorType neww(m,dp);

    static const FloatDPApproximation ALPHA_SCALE_FACTOR = 0.75_approx;

    FloatDPApproximation alpha = 1.0_approx;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        alpha *= ALPHA_SCALE_FACTOR;
    } while ( ! decide(contains(d,newx)) || ! decide(contains(c,neww)) );
    alpha /= ALPHA_SCALE_FACTOR;

    CONCLOG_PRINTLN("alpha="<<alpha);
    CONCLOG_PRINTLN("newx="<<newx<<", neww="<<neww);

    x=newx;
    w=neww;

    return;
}


Void PenaltyFunctionOptimiser::
feasibility_step(ValidatedFeasibilityProblem p,
                 ValidatedVectorType& w, ValidatedVectorType& x) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


// Use a penalty approach without multipliers on the constraint functions
// Solve g(x)=w, x in D, w in C; Lagrangian y.(g(x)-w)
Void PenaltyFunctionOptimiser::
feasibility_step(ApproximateFeasibilityProblem p,
                 ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y) const
{
    CONCLOG_SCOPE_CREATE;
    auto& D=p.D; auto& g=p.g; auto& C=p.C;
    auto m=y.size(); auto n=x.size();

    Vector<ApproximateNumberType> cl=lower_bounds(C);
    Vector<ApproximateNumberType> cu=upper_bounds(C);
    Vector<ApproximateNumberType> dl=lower_bounds(D);
    Vector<ApproximateNumberType> du=upper_bounds(D);

    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    CONCLOG_PRINTLN("dl="<<dl<<", du="<<du);
    CONCLOG_PRINTLN("cl="<<cl<<", cu="<<cu);
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y);

    ARIADNE_ASSERT_MSG(g.argument_size()==D.size(),"D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_ASSERT_MSG(g.result_size()==C.size(),  "D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_ASSERT(w.size()==m);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    // Solve the problem
    //   minimise Sum -log(w-cl)-log(cu-w)-log(x-dl)-log(du-x)
    //   subject to g(x)-w=0

    // Lagrangian
    //   -log(w-cl)-log(cu-w)-log(x-dl)-log(du-x) - y.(g(x)-w)

    // Conditions for a constrained minimum
    // 1/(cu-w)-1/(w-cu) + y       = 0
    // 1/(du-x)-1/(x-dl) - y.Dg(x) = 0
    //          w - g(x)           = 0

    // Second-order conditions
    // (1/(w-cl)^2 + 1/(cu-w)^2) dw                 +   dy = - ( 1/(cu-w) - 1/(w-cl) + y       )
    //   (1/(x-dl)^2 + 1/(du-x)^2 - y.D^2x) dx - Dg'(x) dy = - ( 1/(du-x) - 1/(x-dl) - y.Dg(x) )
    //                           dw - Dg(x) dx             = - ( w - g(x) )

    Vector<Differential<ApproximateNumberType>> ddgx=g.evaluate(Differential<ApproximateNumberType>::variables(2,x));
    CONCLOG_PRINTLN("ddgx="<<ddgx);

    Vector<ApproximateNumberType> gx = ddgx.value();
    CONCLOG_PRINTLN("g(x)="<<gx);
    Matrix<ApproximateNumberType> A = ddgx.jacobian();
    CONCLOG_PRINTLN("Dg(x)="<<A);

    Vector<ApproximateNumberType> yA=transpose(A)*y;

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<ApproximateNumberType> YH(x.size(),x.size(),dp);
    for(SizeType i=0; i!=y.size(); ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    CONCLOG_PRINTLN("Y.D2g(x)="<<YH);

    Vector<ApproximateNumberType> recwu=cu-w; recwu=erec(recwu);
    Vector<ApproximateNumberType> recwl=w-cl; recwl=erec(recwl);
    Vector<ApproximateNumberType> recxu=du-x; recxu=erec(recxu);
    Vector<ApproximateNumberType> recxl=x-dl; recxl=erec(recxl);

    Vector<ApproximateNumberType> diagDw=esqr(recwu)+esqr(recwl);
    Matrix<ApproximateNumberType> Dw(m,m,dp); for(SizeType i=0; i!=m; ++i) { Dw[i][i]=diagDw[i]; }
    DiagonalMatrix<ApproximateNumberType> Dx(esqr(recxu)+esqr(recxl));


    for(SizeType i=0; i!=n; ++i) { YH[i][i]-=Dx[i]; }

    Matrix<ApproximateNumberType> AT=transpose(A);
    Matrix<ApproximateNumberType> Znm(n,m,dp);
    Matrix<ApproximateNumberType> Zmn(m,n,dp);
    Matrix<ApproximateNumberType> Zmm(m,m,dp);
    Matrix<ApproximateNumberType> Im=Matrix<ApproximateNumberType>::identity(m,dp);


    Matrix<ApproximateNumberType> S=cojoin(join(Dw,Zmn,Im),join(Znm,-YH,-AT),join(Im,-A,Zmm));
    Vector<ApproximateNumberType> r=join(recwu-recwl+y,recxu-recxl-yA,w-gx);

    for(SizeType j=0; j!=m; ++j) {
        if( decide(C[j].lower_bound()==C[j].upper_bound()) ) {
            S[j][j]=1;
            S[j][m+n+j]=0;
            S[m+n+j][j]=0;
            r[j]=0;
        }
    }

    Vector<ApproximateNumberType> swxy = -solve(S,r);

    Vector<ApproximateNumberType> sw(m,dp),sx(n,dp),sy(m,dp);
    sw = project(swxy,range(0,m));
    sx = project(swxy,range(m,m+n));
    sy = project(swxy,range(m+n,m+n+m));

    ApproximateNumberType al=one;
    Vector<ApproximateNumberType> nw=w+al*sw;
    Vector<ApproximateNumberType> nx=x+al*sx;
    Vector<ApproximateNumberType> ny(m,dp);
    CONCLOG_PRINTLN("sx="<<sx);
    CONCLOG_PRINTLN("sw="<<sw);
    while( ! decide(contains(C,nw)) || ! decide(contains(D,nx)) ) {
        al*=0.75;
        nw=w+al*sw;
        nx=x+al*sx;
    }
    CONCLOG_PRINTLN("al="<<sw);
    ny=y+al*sy;

    w=nw; x=nx; y=ny;
}


//------- ApproximateOptimiser ------------------------------------------//

ApproximateOptimiser* ApproximateOptimiser::clone() const {
    return new ApproximateOptimiser(*this);
}

ValidatedKleenean ApproximateOptimiser::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", h="<<h);
    ApproximateVectorType x=midpoint(D);
    ApproximateVectorType y(h.result_size(),zero);

    for(SizeType i=0; i!=8; ++i) {
        this->feasibility_step(D,h,x,y);
    }

    if( decide(norm(h(x))<1e-10) ) { return true; }

    if(!possibly(contains(UpperIntervalType(dot(UpperIntervalVectorType(cast_exact(y)),apply(h,D))),zero))) { return false; }

    return indeterminate;
}

Void ApproximateOptimiser::
feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& h,
                 ApproximateVectorType& x, ApproximateVectorType& y) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("x="<<x<<" y="<<y);
    static const double SCALE_FACTOR = 0.75;
    const SizeType n=x.size();
    const SizeType m=y.size();
    // Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));
    FloatDPApproximationMatrix A = ddhx.jacobian();
    CONCLOG_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPApproximationMatrix H(n,n,dp);
    for(SizeType i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(SizeType j=0; j!=n; ++j) {
        H[j][j] += rec(sqr(x[j]-D[j].lower_bound()));
        H[j][j] += rec(sqr(D[j].upper_bound()-x[j]));
    }

    ApproximateVectorType rx = transpose(A) * y;
    for(SizeType j=0; j!=n; ++j) {
        rx[j] -= rec(x[j]-D[j].lower_bound());
        rx[j] += rec(D[j].upper_bound()-x[j]);
    }
    ApproximateVectorType ry = ddhx.value();
    CONCLOG_PRINTLN("rx="<<rx<<" ry="<<ry);

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPApproximationMatrix Hinv=inverse(H);
    CONCLOG_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPApproximationMatrix S=A*Hinv*transpose(A);
    FloatDPApproximationMatrix Sinv=inverse(S);
    CONCLOG_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    ApproximateVectorType dy = Sinv * ( A*(Hinv*rx) - ry );
    ApproximateVectorType dx = Hinv * ( rx - transpose(A) * dy);
    CONCLOG_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPApproximation ax = one;
    ApproximateVectorType nx = x-ax*dx;
    while(!contains(D,cast_exact(nx))) {
        ax*=SCALE_FACTOR;
        nx = x - ax * dx;
    }
    ApproximateVectorType ny = y-ax*dy;
    CONCLOG_PRINTLN("nx="<<nx<<" ax="<<ax<<" ny="<<ny);
    CONCLOG_PRINTLN_AT(1,"h(x)="<<h(nx));

    x=nx; y=ny;
}





//------- KKT / Central path matrices -----------------------------------//

template<class T1, class T2, class T3> decltype(auto) join(T1 const& t1, T2 const& t2, T3 const& t3) {
    return join(join(t1,t2),t3);
}

template<class T1, class T2, class T3, class T4> decltype(auto) join(T1 const& t1, T2 const& t2, T3 const& t3, T4 const& t4) {
    return join(join(t1,t2),join(t3,t4));
}

template<class T, class M> Void assign_dense(MatrixRange<Matrix<T>> R, M const& A) {
    R = A;
}

template<class T, class DM> Void assign_diagonal(MatrixRange<Matrix<T>> R, DM const& D) {
    for (SizeType i=0; i!=D.size(); ++i) { R.at(i,i)=D.at(i,i); }
}

namespace {


// [  Q+Z/X   A' ]
// [   YA     W  ]

template<class T> auto
primal_dual_matrix(SymmetricMatrix<T>const& Q, Matrix<T>const& A,
                   DiagonalMatrix<T>const& W, DiagonalMatrix<T>const& X,
                   DiagonalMatrix<T>const& Y, DiagonalMatrix<T>const& Z) -> Matrix<T>
{
    SizeType m=A.row_size(); SizeType n=A.column_size();
    IdentityMatrix In(n);
    T z=Q.zero_element();
    Matrix<T> S(m+n,m+n,z);
    assign_dense(S[range(0,n)][range(0,n)],Q+Z/X);
    assign_dense(S[range(0,n)][range(n,m+n)],transpose(A));
    assign_dense(S[range(n,m+n)][range(0,n)],Y*A);
    assign_diagonal(S[range(n,m+n)][range(n,m+n)],W);
    return S;
}

struct PrimalDualData {
    using T = FloatDPApproximation;
    Vector<T> x,y;
    PrimalDualData()
        : PrimalDualData(0u,0u,dp) { }
    PrimalDualData(SizeType m, SizeType n, DP pr)
        : x(n,pr), y(m,pr) { }
    PrimalDualData(SizeType m, SizeType n, Vector<T> r)
        : x(r[range(0,n)]), y(r[range(n,m+n)]) { ARIADNE_DEBUG_PRECONDITION(r.size()==m+n); }
    PrimalDualData(Vector<T> x_, Vector<T> y_)
        : x(x_), y(y_) { }
    Vector<T> assemble() const {
        return join(x,y); }
    friend PrimalDualData operator-(PrimalDualData d) {
        return PrimalDualData(-d.x,-d.y); }
};

// [  Q   A' -I  ]
// [ YA   W   0  ]
// [  Z   0   X  ]

template<class T> auto
primal_dual_complementary_matrix(SymmetricMatrix<T>const& Q, Matrix<T>const& A,
                                 DiagonalMatrix<T>const& W, DiagonalMatrix<T>const& X,
                                 DiagonalMatrix<T>const& Y, DiagonalMatrix<T>const& Z) -> Matrix<T>
{
    SizeType m=A.row_size(); SizeType n=A.column_size();
    IdentityMatrix In(n);
    T z=Q.zero_element();
    Matrix<T> S(m+2*n,m+2*n,z);
    assign_dense(S[range(0,n)][range(0,n)],Q);
    assign_dense(S[range(0,n)][range(n,m+n)],transpose(A));
    assign_diagonal(S[range(0,n)][range(m+n,m+2*n)],-In);
    assign_dense(S[range(n,m+n)][range(0,n)],Y*A);
    assign_diagonal(S[range(n,m+n)][range(n,m+n)],W);
    assign_diagonal(S[range(m+n,m+2*n)][range(0,n)],Z);
    assign_diagonal(S[range(m+n,m+2*n)][range(m+n,m+2*n)],X);
    return S;
}

struct PrimalDualComplementaryData {
    using T = FloatDPApproximation;
    Vector<T> x,y,z;
    PrimalDualComplementaryData()
        : PrimalDualComplementaryData(0u,0u,dp) { }
    PrimalDualComplementaryData(SizeType m, SizeType n, DP pr)
        : x(n,pr), y(m,pr), z(n,pr) { }
    PrimalDualComplementaryData(SizeType m, SizeType n, Vector<T> r)
        : x(r[range(0,n)]), y(r[range(n,m+n)]), z(r[range(m+n,m+2*n)]) { ARIADNE_DEBUG_PRECONDITION(r.size()==m+2*n); }
    PrimalDualComplementaryData(Vector<T> x_, Vector<T> y_, Vector<T> z_)
        : x(x_), y(y_), z(z_) { }
    Vector<T> assemble() const {
        return join(x,join(y,z)); }
    friend PrimalDualComplementaryData operator-(PrimalDualComplementaryData d) {
        return PrimalDualComplementaryData(-d.x,-d.y,-d.z); }
};

struct PrimalDualComplementaryMatrix {
    using T = FloatDPApproximation;

    SymmetricMatrix<T> Q; Matrix<T> A; DiagonalMatrix<T> W,X,Y,Z;
    PrimalDualComplementaryMatrix(SymmetricMatrix<T> Q_, Matrix<T> A_,
                                  DiagonalMatrix<T> W_, DiagonalMatrix<T> X_, DiagonalMatrix<T> Y_, DiagonalMatrix<T> Z_)
        : Q(Q_), A(A_), W(W_), X(X_), Y(Y_), Z(Z_) { }

    Matrix<T> assemble() const {
        return primal_dual_complementary_matrix(Q,A,W,X,Y,Z); }

    Void _validate() const {
        auto m=A.row_size(); auto n=A.column_size(); ARIADNE_ASSERT(Q.size()==n);
        ARIADNE_ASSERT(W.size()==m); ARIADNE_ASSERT(Y.size()==m);
        ARIADNE_ASSERT(X.size()==n); ARIADNE_ASSERT(Z.size()==n); }

    PrimalDualComplementaryData solve(PrimalDualComplementaryData const& d) const {
        Vector<T> dr = lu_solve(this->assemble(),d.assemble());
        return PrimalDualComplementaryData(this->A.row_size(),this->A.column_size(),dr);
    }
};

} // namespace

// [  Y   0   W   0  ]
// [  0   Q   A' -I  ]
// [ -I   A   0   0  ]
// [  0   Z   0   X  ]

template<class T> auto
slack_primal_dual_complementary_matrix(SymmetricMatrix<T>const& Q, Matrix<T>const& A,
                                       DiagonalMatrix<T>const& W, DiagonalMatrix<T>const& X, DiagonalMatrix<T>const& Y, DiagonalMatrix<T>const& Z) -> Matrix<T>
{
    SizeType m=A.row_size(); SizeType n=A.column_size();
    IdentityMatrix Im(m), In(n);
    T z=Q.zero_element();
    Matrix<T> S(2*(m+n),2*(m+n),z);
    assign_diagonal(S[range(0,m)][range(0,m)],Y);
    assign_diagonal(S[range(0,m)][range(m+n,2*m+n)],W);
    assign_dense(S[range(m,m+n)][range(m,m+n)],Q);
    assign_dense(S[range(m,m+n)][range(m+n,2*m+n)],transpose(A));
    assign_diagonal(S[range(m,m+n)][range(2*m+n,2*(m+n))],In);
    assign_diagonal(S[range(m+n,2*m+n)][range(0,m)],-Im);
    assign_dense(S[range(m+n,2*m+n)][range(m,m+n)],A);
    assign_diagonal(S[range(2*m+n,2*(m+n))][range(m,m+n)],Z);
    assign_diagonal(S[range(2*m+n,2*(m+n))][range(2*m+n,2*(m+n))],X);
    return S;
}



template<class T> struct InteriorPointOptimiser::SlackPrimalDualComplementaryData {
    Vector<T> w,x,y,z;
    SlackPrimalDualComplementaryData()
        : SlackPrimalDualComplementaryData(0u,0u,dp) { }
    SlackPrimalDualComplementaryData(SizeType m, SizeType n, DP pr)
        : w(m,pr), x(n,pr), y(m,pr), z(n,pr) { }
    SlackPrimalDualComplementaryData(SizeType m, SizeType n, Vector<T> r)
        : w(r[range(0,m)]), x(r[range(m,m+n)]), y(r[range(m+n,2*m+n)]), z(r[range(2*m+n,2*(m+n))]) { ARIADNE_DEBUG_PRECONDITION(r.size()==2*(m+n)); }
    SlackPrimalDualComplementaryData(Vector<T> w_, Vector<T> x_, Vector<T> y_, Vector<T> z_)
        : w(w_), x(x_), y(y_), z(z_) { }
    template<class TT> requires Convertible<TT,T> SlackPrimalDualComplementaryData(SlackPrimalDualComplementaryData<TT>const& d)
        : w(d.w), x(d.x), y(d.y), z(d.z) { }
    Vector<T> assemble() const {
        return join(join(w,x),join(y,z)); }
    friend SlackPrimalDualComplementaryData<T> operator-(SlackPrimalDualComplementaryData<T> d1, SlackPrimalDualComplementaryData<T> d2) {
        return SlackPrimalDualComplementaryData<T>(d1.w-d2.w,d1.x-d2.x,d1.y-d2.y,d1.z-d2.z); }
    friend SlackPrimalDualComplementaryData<T> operator-(SlackPrimalDualComplementaryData<T> d) {
        return SlackPrimalDualComplementaryData<T>(-d.w,-d.x,-d.y,-d.z); }

    friend OutputStream& operator<<(OutputStream& os, SlackPrimalDualComplementaryData<T>& d) {
        return os << "SlackPrimalDualComplementaryData(w=" << d.w << ", x=" << d.x << ", y=" << d.y << ", z=" << d.z << ")"; }
};


//! \brief A matrix defining the linear equations for a full Newton step for slack, primal, dual and complementary variables.
//! \details Form is \f[ \begin{pmatrix}Y&0&W&0\\0&Q&A^T&I\\-I&A&0&0\\0&Z&0&X\end{pmatrix} . \f]
template<class T> struct InteriorPointOptimiser::SlackPrimalDualComplementaryMatrix {
    SymmetricMatrix<T> Q; Matrix<T> A; DiagonalMatrix<T> W,X,Y,Z;
    SlackPrimalDualComplementaryMatrix(SymmetricMatrix<T> Q_, Matrix<T> A_,
                                       DiagonalMatrix<T> W_, DiagonalMatrix<T> X_, DiagonalMatrix<T> Y_, DiagonalMatrix<T> Z_)
        : Q(Q_), A(A_), W(W_), X(X_), Y(Y_), Z(Z_) { }

    Matrix<T> assemble() const {
        return slack_primal_dual_complementary_matrix(Q,A,W,X,Y,Z); }

    Void _validate() const {
        auto m=A.row_size(); auto n=A.column_size(); ARIADNE_ASSERT(Q.size()==n);
        ARIADNE_ASSERT(W.size()==m); ARIADNE_ASSERT(Y.size()==m);
        ARIADNE_ASSERT(X.size()==n); ARIADNE_ASSERT(Z.size()==n); }

    SlackPrimalDualComplementaryData<T> solve(SlackPrimalDualComplementaryData<T> const& d) const {
        Vector<T> dr = lu_solve(this->assemble(),d.assemble());
        return SlackPrimalDualComplementaryData<T>(this->A.row_size(),this->A.column_size(),dr);
    }
};


//------- ApproximateOptimiserBase -----------------------------------//

//------- InteriorPointOptimiser -----------------------------------//

InteriorPointOptimiser::InteriorPointOptimiser() {
}

auto
InteriorPointOptimiser::clone() const -> InteriorPointOptimiser* {
    return new InteriorPointOptimiser(*this);
}

OutputStream& operator<<(OutputStream& os, InteriorPointOptimiser const& opt) {
    return os << "InteriorPointOptimiser()";
}



auto
InteriorPointOptimiser::minimise(ValidatedOptimisationProblem p) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;

    ARIADNE_PRECONDITION(p.D.has_nonempty_interior());
    ARIADNE_PRECONDITION(p.C.has_nonempty_interior());
    ApproximateVectorType x0=p.D.midpoint();
    ApproximateVectorType y0(p.C.dimension(),zero);
    return std::get<1>(this->minimise_hotstarted(p,x0,y0));
}

auto
InteriorPointOptimiser::feasible(ValidatedFeasibilityProblem p) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;

    ARIADNE_PRECONDITION(p.D.has_nonempty_interior());
    ARIADNE_PRECONDITION(p.C.has_nonempty_interior());

    ApproximateVectorType x0=p.D.midpoint();
    ApproximateVectorType y0(p.C.dimension(),zero);

    return std::get<0>(this->feasible_hotstarted(p,x0,y0));
}


auto
InteriorPointOptimiser::minimise(ApproximateOptimisationProblem p) const -> ApproximateVector
{
    CONCLOG_SCOPE_CREATE;

    auto x0 = midpoint(p.D);
    auto y0 = ApproximateVectorType(p.C.dimension(),zero);
    return std::get<1>(minimise_hotstarted(p,x0,y0));
}

auto
InteriorPointOptimiser::minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> Tuple<ValidatedNumber,ValidatedVector,ValidatedVector>
{
    UpperBoxType B=p.D;

    auto r = this->splitting_minimise_hotstarted(p,B,x0,y0);

    auto ro0 = std::get<OptimalityCertificate>(r)[0];
    return make_tuple(ro0.value(),ro0.primal(),ro0.dual());
}


auto
InteriorPointOptimiser::nonsplitting_minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> Tuple<ValidatedNumber,ValidatedVector,ValidatedVector>
{

    auto r = this->minimise_hotstarted(ApproximateOptimisationProblem(p),x0,y0);
    // FIXME: Find a better way of making the answer exact
    const Scalar<ApproximateNumber>& av=std::get<0>(r);
    const Vector<ApproximateNumber>& ax=std::get<1>(r);
    const Vector<ApproximateNumber>& ay=std::get<2>(r);

    Scalar<ValidatedNumber> v=to_generic(cast_exact(Scalar<Approximation<FLT>>(av,pr)));
    Vector<ValidatedNumber> x=to_generic(cast_exact(Vector<Approximation<FLT>>(ax,pr)));
    Vector<ValidatedNumber> y=to_generic(cast_exact(Vector<Approximation<FLT>>(ay,pr)));

    return make_tuple(v,x,y);
}


auto
InteriorPointOptimiser::splitting_minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    UpperBoxType B, ApproximateVectorType x, ApproximateVectorType y) const
        -> Variant<OptimalityCertificate,InfeasibilityCertificate>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    // Extra variables for use of optimization solver
    ApproximateVectorType w=midpoint(C);
    Approximation<FLT> mu = one;
    ApproximateVectorType z=this->compute_z(p,x,mu);

    for(SizeType i=0; i!=18; ++i) {
        CONCLOG_PRINTLN_AT(1,"w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x))
        this->minimisation_step(p, w,x,y,z, mu);
        if (i>=2 && i<17) { mu=mu/4; }
    }
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", mu="<<mu<<", f(x)="<<f(x)<<", g(x)="<<g(x))

    auto r=this->check_minimality(p,w,x,y,z);
    if (definitely(std::get<0>(r))) {
        CONCLOG_PRINTLN("r="<<r);
        FloatDPBounds v=f(cast_exact(x));
        return OptimalityCertificate({LocalOptimalityCertificate(v,cast_exact(x),cast_exact(y))});
    } else {
        ConstraintSolver cslvr;
        UpperBoxType oldB=B;
        cslvr.reduce(B,g,C);
        CONCLOG_PRINTLN("r="<<r<<", D="<<D<<", B="<<oldB<<"="<<B);
        if (definitely(B.is_empty())) {
            return InfeasibilityCertificate({LocalInfeasibilityCertificate(cast_exact_box(oldB),cast_exact(y))});
        } else {
            if (decide(measure(B)<measure(oldB)/2)) {
                return splitting_minimise_hotstarted(p,B,midpoint(B),y);
            } else {
                auto sB=split(B); auto& B1=sB.first; auto& B2=sB.second;
                auto r1=splitting_minimise_hotstarted(p,B1,midpoint(B),y);
                auto r2=splitting_minimise_hotstarted(p,B2,midpoint(B),y);

                if (std::holds_alternative<OptimalityCertificate>(r1) && std::holds_alternative<OptimalityCertificate>(r2)) {
                    // Both sub-boxes are feasible. Filter out clearly non-optimal points
                    auto o=catenate(std::get<OptimalityCertificate>(r1),std::get<OptimalityCertificate>(r2));
                    ARIADNE_DEBUG_ASSERT(not o.empty());
                    auto v=o[0].v;
                    for (auto lo : o) { v=min(v,lo.v); }
                    OptimalityCertificate ro;
                    for (auto lo : o) { if (possibly(lo.v<=v)) { ro.append(lo); } }
                    return ro;
                } else if (std::holds_alternative<OptimalityCertificate>(r1)) {
                    // B1 is feasible, but B2 is not
                    return r1;
                } else if (std::holds_alternative<OptimalityCertificate>(r2)) {
                    return r2;
                } else {
                    return InfeasibilityCertificate(catenate(std::get<InfeasibilityCertificate>(r1),std::get<InfeasibilityCertificate>(r2)));
                }
            }
        }
    }

}

auto
InteriorPointOptimiser::minimise_hotstarted(
    const ApproximateOptimisationProblem& p,
    const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const
        -> Tuple<ApproximateNumber,ApproximateVector,ApproximateVector>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    CONCLOG_PRINTLN("f="<<f<<", D="<<D<<", g="<<g<<", C="<<C);

    StepData d(g.result_size(),g.argument_size(),zero.precision());
    auto& w=d.w; auto& x=d.x; auto& y=d.y; auto& z=d.z;

    d.mu = one;
    d.x = x0;
    d.y = y0;
    d.w = midpoint(C);
    d.z = this->compute_z(p,x0,d.mu);

    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x));
        this->minimisation_step(p, d);
        if (2<=i && i<=10) { d.mu *= 0.25_exact; }
    }
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z)

    return make_tuple(f(d.x),d.x,d.y);
}


auto
InteriorPointOptimiser::feasible(ApproximateFeasibilityProblem p) const -> ApproximateKleenean
{
    CONCLOG_SCOPE_CREATE;

    ApproximateVectorType x0 = p.D.midpoint();
    ApproximateVectorType y0(p.g.result_size(),zero);
    return std::get<0>(this->feasible_hotstarted(p,x0,y0));
}



auto
InteriorPointOptimiser::splitting_feasible_hotstarted(
    const ValidatedFeasibilityProblem& p,
    UpperBoxType B, ApproximateVectorType x, ApproximateVectorType y) const
        -> Variant<FeasibilityCertificate,InfeasibilityCertificate>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    // Extra variables for use of optimization solver
    ApproximateVectorType w=midpoint(C);
    Approximation<FLT> mu = one;
    ApproximateVectorType z=this->compute_z(p,x,mu);

    // FIXME: Allow more steps
    Approximation<FLT> t=this->compute_t(p,x);

    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        this->feasibility_step(p, w,x,y,z);
        t=this->compute_t(p,x);
        if(probably(LogicalValue(t>0))) {
            break;
        }
    }
    CONCLOG_PRINTLN("t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x))

    auto r=this->check_feasibility(p,x,y);
    if (definitely(std::get<0>(r))) {
        CONCLOG_PRINTLN("r="<<r);
        return FeasibilityCertificate(Pair<ValidatedVector,ExactVector>(cast_exact(x),cast_exact(y)));
    } else {
        ConstraintSolver cslvr;
        UpperBoxType oldB=B;
        cslvr.reduce(B,g,C);
        CONCLOG_PRINTLN("r="<<r<<", D="<<D<<", B="<<oldB<<"="<<B);
        if (definitely(B.is_empty())) {
            return InfeasibilityCertificate({LocalInfeasibilityCertificate(cast_exact_box(oldB),cast_exact(y))});
        } else {
            if (decide(measure(B)<measure(oldB)/2)) {
                return splitting_feasible_hotstarted(p,B,midpoint(B),y);
            } else {
                auto sB=split(B);
                auto r1=splitting_feasible_hotstarted(p,sB.first,midpoint(B),y);
                if (std::holds_alternative<FeasibilityCertificate>(r1)) { return r1; }
                auto r2=splitting_feasible_hotstarted(p,sB.second,midpoint(B),y);
                if (std::holds_alternative<FeasibilityCertificate>(r2)) { return r2; }
                return InfeasibilityCertificate(catenate(std::get<InfeasibilityCertificate>(r1),std::get<InfeasibilityCertificate>(r2)));
            }
        }
    }

}

auto
InteriorPointOptimiser::feasible_hotstarted(
    const ValidatedFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p<<", x0="<<x0<<", y0="<<y0)

    auto r=splitting_feasible_hotstarted(p,p.D,x0,y0);
    if (std::holds_alternative<FeasibilityCertificate>(r)) {
        FeasibilityCertificate fc=std::get<FeasibilityCertificate>(r);
        CONCLOG_PRINTLN("fc="<<fc)
        return make_tuple(true,fc.x,fc.y);
    } else {
        InfeasibilityCertificate ifc=std::get<InfeasibilityCertificate>(r);
        CONCLOG_PRINTLN("ifc="<<ifc)
        return make_tuple(false,cast_singleton(ifc[0].B),ifc[0].y);
    }
}

auto
InteriorPointOptimiser::nonsplitting_feasible_hotstarted(
    const ValidatedFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ARIADNE_ASSERT(x0.size()==g.argument_size());
    ARIADNE_ASSERT(y0.size()==g.result_size());

    Vector<Approximation<FLT>> x=x0;
    Vector<Approximation<FLT>> y=y0;

    // Extra variables for use of optimization solver
    Vector<Approximation<FLT>> w=midpoint(C);
    Approximation<FLT> mu = one;
    Vector<Approximation<FLT>> z=this->compute_z(p,x0,mu);

    // FIXME: Allow more steps
    Approximation<FLT> t=this->compute_t(p,x);

    CONCLOG_PRINTLN_AT(1,"D="<<D<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        CONCLOG_PRINTLN_AT(2,"w="<<w<<", z="<<z);
        this->feasibility_step(p, w,x,y,z);
        t=this->compute_t(p,x);
        if(probably(LogicalValue(t>0))) {
            return this->check_feasibility(p,x,y);
        }
    }
    CONCLOG_PRINTLN("t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    return this->check_feasibility(p,x,y);
}


auto
InteriorPointOptimiser::check_feasibility(
    const ValidatedFeasibilityProblem& p,
    const ApproximateVectorType& xa, const ApproximateVectorType& ya) const
        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    CONCLOG_SCOPE_CREATE;

    auto& g=p.g;
    const Vector<Exact<FLT>> x=cast_exact(xa);
    const Vector<Exact<FLT>> y=cast_exact(ya);

    CONCLOG_PRINTLN("p="<<p)
    CONCLOG_PRINTLN("x="<<x<<", y="<<y<<", g(x)="<<g(x))

    FeasibilityChecker fc;
    ValidatedKleenean k=fc.check_feasibility(p,x,y);
    return make_tuple(k,to_generic(x),to_generic(y));
}


auto
InteriorPointOptimiser::check_minimality(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& wa, const ApproximateVectorType& xa, const ApproximateVectorType& ya, const ApproximateVectorType& za) const
        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    //CONCLOG_SCOPE_CREATE;
    auto& f=p.f; auto& g=p.g;

    CONCLOG_PRINTLN("InteriorPointOptimiser::check_minimality(ValidatedOptimisationProblem p, ApproximateVectorType wa, xa, ya, za)");
    CONCLOG_PRINTLN("wa="<<wa<<", xa="<<xa<<", ya="<<ya<<", za="<<za<<", f(xa)="<<f(xa)<<", g(xa)="<<g(xa));

    ExactDouble se(std::numeric_limits<float>::epsilon());
    se=ExactDouble(1.0/1024);
    se=ExactDouble(1.0/1048576);
    FloatDPBounds e=se*FloatDPBounds(-1,+1,dp);
    auto cast_exact_widen = [&e](auto va){return eadd(cast_exact(va)*(1+e),e);};

    FloatDPBoundsVector w=cast_exact_widen(wa);
    FloatDPBoundsVector x=cast_exact_widen(xa);
    FloatDPBoundsVector y=cast_exact_widen(ya);
    FloatDPBoundsVector z=cast_exact_widen(za);

    // auto restrict_number = [](FloatDPBounds xb, UpperIntervalType ivl) {
    //     return FloatDPBounds(max(xb.lower(),ivl.lower_bound()),min(xb.upper(),ivl.upper_bound())); };
    // auto restrict_vector = [&](Vector<FloatDPBounds> v, ExactBoxType bx) {
    //     return Vector<FloatDPBounds>( v.size(), [&](SizeType i){return restrict_number(v[i],bx[i]);} ); };
    // w=restrict_vector(w,p.C);
    // x=restrict_vector(x,p.D);

    CONCLOG_PRINTLN_AT(0, "w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x));
    return this->check_minimality(p,w,x,y,z);
}

auto
InteriorPointOptimiser::check_minimality(
    const ValidatedOptimisationProblem& p,
    const FloatDPBoundsVector& w_, const FloatDPBoundsVector& x_, const FloatDPBoundsVector& y_, const FloatDPBoundsVector& z_) const
        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    Vector<Bounds<FLT>> w(w_), x(x_), y(y_), z(z_);

    //CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("InteriorPointOptimiser::check_minimality(ValidatedOptimisationProblem p, Vector<Bounds<FLT>> w, x, y, z)");
    CONCLOG_PRINTLN_AT(1, "w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x));

    // Set-up interval Newton contractor
    // Want to solve KKT conditions:
    // (w-cl)(cu-w)y=0
    // df(x)+dg(x)'*y+z=0
    // g(x)-w=0
    // (x-dl)(du-x)z=0

    // [ Y  0  W  0]
    // [ 0  Q AT  I]
    // [-I  A  0  0]
    // [ 0  Z  0  X]

    auto Cl=lower_bounds(C);
    auto Cu=upper_bounds(C);
    auto Dl=lower_bounds(D);
    auto Du=upper_bounds(D);

    auto midpoint = [](ExactIntervalType const& ivl) {
        return FloatDPBounds(hlf(ivl.lower_bound()+ivl.upper_bound())); };
    auto midpoints = [&midpoint](ExactBoxType const& bx) {
        return Vector<FloatDPBounds>( bx.dimension(), [&midpoint,&bx](SizeType i){return midpoint(bx[i]);} ); };

    bool is_refinement = false;
    int step=0;

    while (not is_refinement && step<8) {
        auto ddfx=differential(f,x,2);
        auto ddgx=differential(g,x,2);

        auto A=ddgx.jacobian();
        auto Q=ddfx.hessian();
        for (SizeType i=0; i!=y.size(); ++i) { Q+=y[i]*ddgx[i].hessian(); }

        auto W=DiagonalMatrix<Bounds<FLT>>(emul(w-Cl,Cu-w));
        auto X=DiagonalMatrix<Bounds<FLT>>(emul(x-Dl,Du-x));
        auto Y=DiagonalMatrix<Bounds<FLT>>(emul(2,midpoints(C)-w,y));
        auto Z=DiagonalMatrix<Bounds<FLT>>(emul(2,midpoints(D)-x,z));

        Vector<Exact<FLT>> we=cast_exact(w);
        Vector<Exact<FLT>> xe=cast_exact(x);
        Vector<Exact<FLT>> ye=cast_exact(y);
        Vector<Exact<FLT>> ze=cast_exact(z);
        Covector<Bounds<FLT>> ce=f.gradient(xe);
        Matrix<Bounds<FLT>> Ae=g.jacobian(xe);

        Vector<Bounds<FLT>> rw=emul(we-Cl,Cu-we,ye);
        Vector<Bounds<FLT>> rx=transpose(ce)+transpose(Ae)*ye+ze;
        Vector<Bounds<FLT>> ry=g(xe)-we;
        Vector<Bounds<FLT>> rz=emul(xe-Dl,Du-xe,ze);

        SlackPrimalDualComplementaryMatrix<Bounds<FLT>> S(Q,A,W,X,Y,Z);
        SlackPrimalDualComplementaryData<Bounds<FLT>> r(rw,rx,ry,rz);
        SlackPrimalDualComplementaryData<Exact<FLT>> q(we,xe,ye,ze);

        auto Sa=Matrix<Approximation<FLT>>(S.assemble());
        CONCLOG_PRINTLN("step="<<step);
        CONCLOG_PRINTLN("C="<<C<<", D="<<D);
        CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z);
        CONCLOG_PRINTLN("Sa="<<pretty(Sa))

        SlackPrimalDualComplementaryData<Bounds<FLT>> nd=q-S.solve(r);
        auto& nw=nd.w; auto& nx=nd.x; auto& ny=nd.y; auto& nz=nd.z;
        CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z)
        CONCLOG_PRINTLN("nw="<<nw<<", nx="<<nx<<", ny="<<ny<<", nz="<<nz)
        CONCLOG_PRINTLN("refines(nw,w)="<<refines(nw,w)<<", refines(nx,x)="<<refines(nx,x)<<", refines(ny,y)="<<refines(ny,y)<<", refines(nz,z)="<<refines(nz,z));
        CONCLOG_PRINTLN("inconsistent(nw,w)="<<inconsistent(nw,w)<<", inconsistent(nx,x)="<<inconsistent(nx,x)<<", inconsistent(ny,y)="<<inconsistent(ny,y)<<", inconsistent(nz,z)="<<inconsistent(nz,z));
        if (refines(nw,w) && refines(nx,x) && refines(ny,y) && refines(nz,z)) {
            return make_tuple(true,nx,ny);
        } else if (inconsistent(nw,w) || inconsistent(nx,x) || inconsistent(ny,y) || inconsistent(nz,z)) {
            return make_tuple(false,x,ny);
        } else {
            w=refinement(w,nw); x=refinement(x,nx); y=refinement(y,ny); z=refinement(z,nz);
            ++step;
        }
    }
    return make_tuple(indeterminate,x,y);
}

auto
InteriorPointOptimiser::feasible_hotstarted(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> Tuple<ApproximateKleenean,ApproximateVector,ApproximateVector>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& g=p.g; auto& C=p.C;

    ARIADNE_ASSERT(x0.size()==g.argument_size());
    ARIADNE_ASSERT(y0.size()==g.result_size());

    Vector<Approximation<FLT>> x=x0;
    Vector<Approximation<FLT>> y=y0;

    // Extra variables for use of optimization solver
    Vector<Approximation<FLT>> w=midpoint(C);
    Approximation<FLT> mu = one;
    Vector<Approximation<FLT>> z=this->compute_z(p,x0,mu);

    // FIXME: Allow more steps
    Approximation<FLT> t=this->compute_t(p,x);

    for (SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        this->feasibility_step(p, w,x,y,z);
        t=this->compute_t(p,x);
        if (probably(LogicalValue(t>0))) {
            return make_tuple(true,x,y);
        }
    }
    CONCLOG_PRINTLN("x="<<x<<", y="<<y<<", g(x)="<<g(x));
    return make_tuple(false,x,y);
}





// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used


// min f(x) | x\in D & w\in C | g(x) = w
// Lagrange multipliers y d(g(x)-w); z dx
Void
InteriorPointOptimiser::minimisation_step(const ApproximateOptimisationProblem& p, StepData& d) const
{
    auto& w=d.w; auto& x=d.x; auto& y=d.y; auto& z=d.z; auto& mu=d.mu;
    this->minimisation_step(p, w,x,y,z,mu);
}

auto
InteriorPointOptimiser::minimisation_update(
    const ApproximateOptimisationProblem& p,
    ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z,
    Approximation<FLT>& mu) const -> SlackPrimalDualComplementaryData<Approximation<FLT>>
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    const SizeType m=g.result_size();
    const SizeType n=g.argument_size();

    ARIADNE_DEBUG_PRECONDITION(w.size()==y.size());
    ARIADNE_DEBUG_PRECONDITION(x.size()==z.size());
    ARIADNE_DEBUG_PRECONDITION(w.size()==m);
    ARIADNE_DEBUG_PRECONDITION(x.size()==n);

    ARIADNE_DEBUG_PRECONDITION(decide(contains(D,x)));
    ARIADNE_DEBUG_PRECONDITION(decide(contains(C,w)));
    ARIADNE_DEBUG_PRECONDITION(decide(mu>0));

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("x="<<x);
    CONCLOG_PRINTLN("w="<<w);
    CONCLOG_PRINTLN("y="<<y);
    CONCLOG_PRINTLN("z="<<z);
    CONCLOG_PRINTLN("mu="<<mu);


    Differential<Approximation<FLT>> ddfx=f.evaluate(Differential<Approximation<FLT>>::variables(2,x));
    Vector<Differential<Approximation<FLT>>> ddgx=g.evaluate(Differential<Approximation<FLT>>::variables(2,x));

    // G is the constraint value vector
    Approximation<FLT> fx = ddfx.value();
    Vector<Approximation<FLT>> gx = ddgx.value();
    CONCLOG_PRINTLN("f(x)="<<fx);
    CONCLOG_PRINTLN("g(x)="<<gx);
    CONCLOG_PRINTLN_AT(1,"g(x)-w="<<(gx-w));

    // A, B are the derivative matrices aij=dgi/dxj
    // HACK: Need to explicitly set size of Jacobian if g or h have result_size of zero
    Vector<Approximation<FLT>> c = transpose(ddfx.gradient());
    CONCLOG_PRINTLN_AT(1,"c=df(x)="<<c);
    Matrix<Approximation<FLT>> A = ddgx.jacobian();
    if(m==0) { A=Matrix<Approximation<FLT>>(m,n,dp); }
    CONCLOG_PRINTLN("A=Dg(x)="<<A);

    // Q is the Hessian matrix Q[i1,i2] = df/dx[i1]dx[i2] + Sum_[j]y[j]*dg[j]/dx[i1]dx[i2]
    Matrix<Approximation<FLT>> Q = ddfx.hessian();
    for(SizeType j=0; j!=m; ++j) { Q += y[j] * ddgx[j].hessian(); }
    CONCLOG_PRINTLN("Q="<<Q);

    // Determines the weighting to give to the relaxation parameter mu
    // for equality constraints relative to other constraints
    // Currently unused
    //    static const double EQUALITY_RELAXATION_MULTIPLIER = 1.0;

    // Set-up equations
    //   (w-cl)*(cu-w)*y + (cl+cu-2*w) * mu = 0
    //   Df(x) + sum y_j Dg_j(x) + z = 0
    //   g(x) - w = 0
    //   (x-dl)*(du-x)*z + (dl+du-2*w) * mu = 0

    auto Cl = lower_bounds(C);
    auto Cu = upper_bounds(C);
    auto Dl = lower_bounds(D);
    auto Du = upper_bounds(D);

    // Compute the residuals and contributions from slack in x and w
    //   rx[i] = df/dx[i] + Sum[j] dg[j]/dx[i] * y[j] + z[i]
    Vector<Approximation<FLT>> rw = emul(w-Cl,Cu-w,y) + (Cl+Cu-2*w)*mu;
    Vector<Approximation<FLT>> rx = c + transpose(A) * y + z;
    Vector<Approximation<FLT>> ry = gx - w;
    Vector<Approximation<FLT>> rz = emul(x-Dl,Du-x,z) + (Dl+Du-2*x)*mu;

    DiagonalMatrix<Approximation<FLT>> W(emul(w-Cl,Cu-w));
    DiagonalMatrix<Approximation<FLT>> Y(esub(emul(Cl+Cu-2*w,y),2*mu));
    DiagonalMatrix<Approximation<FLT>> X(emul(x-Dl,Du-x));
    DiagonalMatrix<Approximation<FLT>> Z(esub(emul(Dl+Du-2*x,z),2*mu));

    // FIXME: Make this more consistent
    for (SizeType j=0; j!=m; ++j) {
        if (decide(Cl[j]==-inf)) { rw[j]=(Cu[j]-w[j])*y[j]-mu; W._at(j)=Cu[j]-w[j]; Y._at(j)=y[j]; }
        else if (decide(Cu[j]==+inf)) { rw[j]=(w[j]-Cl[j])*y[j]-mu; W._at(j)=w[j]-Cl[j]; Y._at(j)=y[j]; }
    }

    CONCLOG_PRINTLN("rw="<<rw);
    CONCLOG_PRINTLN("rx="<<rx);
    CONCLOG_PRINTLN("ry="<<ry);
    CONCLOG_PRINTLN("rz="<<rz);

    SlackPrimalDualComplementaryMatrix<Approximation<FLT>> S(Q,A,W,X,Y,Z);
    SlackPrimalDualComplementaryData<Approximation<FLT>> r(rw,rx,ry,rz);

    SlackPrimalDualComplementaryData<Approximation<FLT>> dr=S.solve(r);

    return dr;
}

Void
InteriorPointOptimiser::minimisation_step(
    const ApproximateOptimisationProblem& p,
    ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z,
    Approximation<FLT>& mu) const
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ARIADNE_DEBUG_PRECONDITION(w.size()==g.result_size());
    ARIADNE_DEBUG_PRECONDITION(x.size()==f.argument_size());
    ARIADNE_DEBUG_PRECONDITION(y.size()==w.size());
    ARIADNE_DEBUG_PRECONDITION(z.size()==x.size());

    ARIADNE_DEBUG_PRECONDITION(decide(contains(D,x)));
    ARIADNE_DEBUG_PRECONDITION(decide(contains(C,w)));
    ARIADNE_DEBUG_PRECONDITION(decide(mu>0));

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", mu="<<mu);

    SlackPrimalDualComplementaryData<Approximation<FLT>> dr = this->minimisation_update(p, w,x,y,z, mu);

    static const Approximation<FLT> ALPHA_SCALE_FACTOR = 0.75_approx;
    static const Approximation<FLT> MINIMUM_ALPHA = 1e-16_approx;

    // Compute distance to move variables preserving feasibility
    // FIXME: Current implementation might fail due to getting too close to boundary!
    ApproximateVectorType newx(x.size(),dp);
    ApproximateVectorType neww(w.size(),dp);
    Approximation<FLT> alpha = 1.0_approx;
    Bool success = false;
    do {
        newx = x - alpha * dr.x;
        neww = w - alpha * dr.w;
        if (probably(contains(D,newx)) && probably(contains(C,neww))) { success = true; }
        else { alpha *= ALPHA_SCALE_FACTOR; }
        if (probably(alpha<MINIMUM_ALPHA)) { throw NearBoundaryOfFeasibleDomainException(); }
    } while (!success);
    CONCLOG_PRINTLN("alpha="<<alpha);

    ApproximateVectorType newz = z - alpha * dr.z;
    ApproximateVectorType newy = y - alpha * dr.y;

    CONCLOG_PRINTLN("neww="<<neww<<", newx="<<newx<<", newy="<<newy<<", newz="<<newz);

    x=newx; w=neww; y=newy; z=newz;
}


Void
InteriorPointOptimiser::feasibility_step(
    const ApproximateFeasibilityProblem& p,
    ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z) const
{
    CONCLOG_SCOPE_CREATE;

    EffectiveScalarMultivariateFunction f(p.g.argument_size());
    ApproximateOptimisationProblem optp(f,p.D,p.g,p.C);

    Approximation<FLT> mu=one;
    this->minimisation_step(optp, w,x,y,z, mu);
}



auto
InteriorPointOptimiser::compute_dual(
    const ApproximateBoxType& D,
    const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    ApproximateVectorType z(x.size(),x.zero_element());
    for(SizeType i=0; i!=D.size(); ++i) {
        if (decide(D[i].lower_bound() == D[i].upper_bound())) { }
        else if (decide(D[i].lower_bound()==-infty)) { z[i] = mu / (D[i].upper_bound()-x[i]); }
        else if (decide(D[i].upper_bound()==+infty)) { z[i] = -mu / (x[i]-D[i].lower_bound()); }
        else { z[i] = mu * ( rec(D[i].upper_bound()-x[i]) - rec(x[i]-D[i].lower_bound()) ); }
    }
    return z;
}

auto
InteriorPointOptimiser::initial_step_data(
    const ApproximateFeasibilityProblem& p) const
        -> StepData
{
    ApproximateVectorType x=this->compute_x(p);
    ApproximateVectorType y=this->compute_y(p,x,one);
    return this->initial_step_data_hotstarted(p,x,y);
}

auto
InteriorPointOptimiser::initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x, const ApproximateVectorType& y) const
        -> StepData
{
    Approximation<FLT> mu=one;
    ApproximateVectorType w=this->compute_w(p,x,y,mu);
    ApproximateVectorType z=this->compute_z(p,x,mu);
    return StepData(w,x,y,z,mu);
}


auto
InteriorPointOptimiser::compute_x(const ApproximateFeasibilityProblem& p) const
        -> ApproximateVectorType
{
    return p.D.midpoint();
}

auto
InteriorPointOptimiser::compute_y(const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    return this->compute_dual(p.C,p.g(x),mu);
}

auto
InteriorPointOptimiser::compute_w(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const ApproximateVectorType& y, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    return this->compute_w(p,x,mu);
}

auto
InteriorPointOptimiser::compute_w(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    static const ApproximateDouble GAMMA=0.0009765625; // 1.0/1024;

    auto& C=p.C;
    auto gx=p.g(x);

    ApproximateVectorType w(gx.size(),gx.zero_element());
    for(SizeType j=0; j!=C.size(); ++j) {
        if (decide(C[j].radius()<=mu*GAMMA)) {
            w[j]=C[j].midpoint();
        } else {
            w[j]=min( max( gx[j], C[j].lower_bound()+mu*GAMMA ), C[j].upper_bound()-mu*GAMMA );
        }
    }
    return w;
}

auto
InteriorPointOptimiser::compute_z(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    return this->compute_dual(p.D,x,mu);
}

auto
InteriorPointOptimiser::compute_t(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x) const
        -> Approximation<FLT>
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    // Compute a the minimum constraint satisfaction
    Approximation<FLT> t=zero; t=+infty;
    ApproximateVectorType gx = g(x);

    for(SizeType i=0; i!=D.size(); ++i) {
        t = min( t, min(x[i]-D[i].lower_bound(), D[i].upper_bound()-x[i]) );
    }
    for(SizeType j=0; j!=C.size(); ++j) {
        t = min( t, min(gx[j]-C[j].lower_bound(), C[j].upper_bound()-gx[j]) );
    }
    return t;
}

auto
InteriorPointOptimiser::compute_mu(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const ApproximateVectorType& y) const
        -> ApproximateNumberType
{
    auto& g=p.g; auto& C=p.C;

    // Compute the relaxation parameter mu as the average of the product of the Lyapunov exponents and constraint satisfactions
    Approximation<FLT> mu=zero;
    ApproximateVectorType gx = g(x);

    for(SizeType i=0; i!=C.size(); ++i) {
        if (decide(C[i].lower_bound()==C[i].upper_bound())) { }
        else if (decide(C[i].lower_bound()==-infty)) { mu += y[i] * (gx[i] - C[i].upper_bound()); }
        else if (decide(C[i].upper_bound()==+infty)) { mu += y[i] * (gx[i] - C[i].lower_bound()); }
        else {
            if ( decide(y[i] <=0.0) ) { mu += y[i] * (gx[i] - C[i].upper_bound()); }
            else { mu += y[i] * (gx[i] - C[i].lower_bound()); }
        }
    }
    mu /= C.size();
    return mu;
}


Void InteriorPointOptimiser::compute_tz(
    const ApproximateFeasibilityProblem& p,
    ApproximateVectorType& x, ApproximateNumberType& t, ApproximateVectorType& z) const
{
    t = this->compute_t(p,x);
    z = this->compute_z(p,x,one);
}



//------- PrimalDualOnlyInteriorPointOptimiser -------------------------//

auto PrimalDualOnlyInteriorPointOptimiser::clone() const -> PrimalDualOnlyInteriorPointOptimiser* {
    return new PrimalDualOnlyInteriorPointOptimiser(*this);
}

auto PrimalDualOnlyInteriorPointOptimiser::
feasible_hotstarted(const ApproximateFeasibilityProblem& p,
                    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
                        -> Tuple<ApproximateKleenean,ApproximateVector,ApproximateVector>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& g=p.g;

    ARIADNE_ASSERT(x0.size()==g.argument_size());
    ARIADNE_ASSERT(y0.size()==g.result_size());

    Vector<Approximation<FLT>> x=x0;
    Vector<Approximation<FLT>> y=y0;
    Approximation<FLT> t(dp);

    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        this->feasibility_step(p,x,y,t);
        if(probably(LogicalValue(t>0))) {
            return make_tuple(true,to_generic(x),to_generic(y));
        }
    }
    CONCLOG_PRINTLN("t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    return make_tuple(false,x,y);
}


Void
PrimalDualOnlyInteriorPointOptimiser::feasibility_step(
    const ApproximateFeasibilityProblem& p,
    ApproximateVectorType& x, ApproximateVectorType& y, FloatDPApproximation& t) const
{
    CONCLOG_SCOPE_CREATE;

    auto& d=p.D; auto& g=p.g; auto& c=p.C;

    static const ApproximateDouble GAMMA=0.0009765625; // 1.0/1024;
    static const ApproximateDouble SIGMA=0.125;
    static const ApproximateDouble SCALE=0.75;

    const SizeType m=d.size();
    const SizeType n=c.size();
    CONCLOG_PRINTLN("m="<<m<<" n="<<n);

    ApproximateVectorType z(n,dp);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    CONCLOG_PRINTLN("ddgx="<<ddgx);

    Vector<FloatDPApproximation> gx = ddgx.value();
    CONCLOG_PRINTLN("g(x)="<<gx<<" ");
    Matrix<FloatDPApproximation> A = transpose(ddgx.jacobian());
    CONCLOG_PRINTLN("A="<<A<<" ");

    // Q is the Hessian matrix of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) \lambda_k$
    Matrix<FloatDPApproximation> Q(n,n,dp);
    for(SizeType i=0; i!=m; ++i) {
        Q+=y[i]*ddgx[i].hessian();
    }
    CONCLOG_PRINTLN("Q="<<Q<<" ");


    // Add correction for singleton domain to diagonal elements of Hessian
    for(SizeType i=0; i!=m; ++i) {
    }

    // Compute diagonal entries of KKT Hessian
    Vector<FloatDPApproximation> D(m,dp);
    for(SizeType j=0; j!=m; ++j) {
        if (decide(c[j].lower_bound()==c[j].upper_bound())) {
        } else if (decide(c[j].upper_bound()==+inf)) {
        } else if (decide(c[j].lower_bound()==-inf)) {
        } else {
            ARIADNE_DEBUG_ASSERT(decide(-infty<c[j].lower_bound() && c[j].lower_bound()<c[j].upper_bound() && c[j].upper_bound()<+infty));
        }
    }

    FloatDPApproximation sigma(SIGMA,dp);
    FloatDPApproximation mu=dot(x,z)/n;
    if(!egtr(emul(x,z),GAMMA*mu)) {
        CONCLOG_PRINTLN("WARNING: near-degeneracy in Lyapunov multipliers in interior-point solver:");
        CONCLOG_PRINTLN("x="<<x<<", y="<<y<<", z="<<z);
        x=(1-sigma)*x+ApproximateVectorType(n,sigma/n);
        mu=dot(x,z)/n;
    }

    ApproximateVectorType yt=join(y,t);
    CONCLOG_PRINTLN("x="<<x<<" yt="<<yt<<" z="<<z);


    // Construct diagonal matrices
    ApproximateVectorType DE=ediv(x,z);
    CONCLOG_PRINTLN("D="<<DE);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    ApproximateVectorType gye(2*(m+n),dp);
    //for(SizeType j=0; j!=n; ++j) { gxe[j]=gy[j]-c[j].upper_bound()+t; gye[n+j]=c[j].lower_bound()-gy[j]+t; }
    //for(SizeType i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper_bound()+t; gye[2*n+m+i]=d[i].lower_bound()-y[i]+t; }
    CONCLOG_PRINTLN("GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    FloatDPApproximationMatrix AE(m+1,2*(m+n),dp);
    //for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    //for(SizeType i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    //for(SizeType k=0; k!=o; ++k) { AE[m][k]=1; }
    FloatDPApproximationMatrix AET=transpose(AE);

    // Construct the symmetric matrix and its inverse
    //FloatDPMatrix S(m+1,m+1); adat(S,AE,DE);
    //CONCLOG_PRINTLN("S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //CONCLOG_PRINTLN("S="<<S);
    FloatDPApproximationMatrix S=feasibility_adat(Q,A,DE);
    CONCLOG_PRINTLN("S="<<S);
    FloatDPApproximationMatrix Sinv=inverse(S);
    CONCLOG_PRINTLN("Sinv="<<Sinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    ApproximateVectorType rx=esub(emul(x,z),mu*sigma);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    ApproximateVectorType ryt=-feasibility_mul(A,x); ryt[m]+=1; // FIXME: Need hessian
    ApproximateVectorType rz=gye+z;
    CONCLOG_PRINTLN("rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;
    ApproximateVectorType rr=ryt + AE*ediv(ApproximateVectorType(rx-emul(x,rz)),z) - ryt;


    // Compute the differences
    ApproximateVectorType dyt=Sinv*rr;
    //RawFloatDPVector dz=-rz-prod(AET,dyt);
    ApproximateVectorType dz=-rz-feasibility_trmul(A,dyt);
    ApproximateVectorType dx=-ediv(ApproximateVectorType(rx+emul(x,dz)),z);
    CONCLOG_PRINTLN("dx="<<dx<<" dyt="<<dyt<<" dz="<<dz);

    ApproximateVectorType nx,ny,nyt,nz; FloatDPApproximation nt(dp);

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly rather than assuming the
    // linear update of y and z are good enough.
    Bool allpositive=false;
    FloatDPApproximation alpha=1/FloatDPApproximation(SCALE,dp);
    if(!egtr(emul(x,z) , GAMMA*mu/16)) {
        CONCLOG_PRINTLN("WARNING: x="<<x<<", z="<<z<< ", x.z="<<emul(x,z)<<"<"<<GAMMA*mu / 16);
        throw NearBoundaryOfFeasibleDomainException();
    }
    while(!allpositive) {
        alpha=alpha*SCALE;
        nx=x+alpha*dx;
        nyt=yt+alpha*dyt;
        ny=project(nyt,range(0,m));
        nt=nyt[m];
        //InteriorPointOptimiser::compute_z(d,g,c,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),GAMMA*mu);
    }
    CONCLOG_PRINTLN("alpha="<<alpha);
    CONCLOG_PRINTLN("nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<emul(nx,nz));

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];
}


//------- PrimalSplitDualInteriorPointOptimiser -------------------------//


auto
PrimalSplitDualInteriorPointOptimiser::minimise_hotstarted(
    const ApproximateOptimisationProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> Tuple<ApproximateNumber,ApproximateVector,ApproximateVector>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    CONCLOG_PRINTLN("f="<<f<<", D="<<D<<", g="<<g<<", C="<<C)

    StepData d(g.result_size(),g.argument_size(),zero.precision());
    auto& x=d.x;

    d.mu = one;
    d.x = x0;

    d.w = this->compute_w(p,d.x,d.yu-d.yl,d.mu);

    d.yl = ediv(d.mu,d.w-lower_bounds(C));
    d.yu = ediv(d.mu,upper_bounds(C)-d.w);
    d.zl = ediv(d.mu,d.x-lower_bounds(D));
    d.zu = ediv(d.mu,upper_bounds(D)-d.x);

    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"w="<<d.w<<", d.x="<<d.x<<", d.yl="<<d.yl<<", d.yu="<<d.yu<<", d.zl="<<d.zl<<", d.zu="<<d.zu)
        CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", g(x)="<<g(x));
        this->minimisation_step(p, d.w,d.x,d.yl,d.yu,d.zl,d.zu, d.mu);
        if(i%3==0 && i<=10) { d.mu *= 0.25_exact; }
    }
    CONCLOG_PRINTLN_AT(0,"w="<<d.w<<", d.x="<<d.x<<", d.yl="<<d.yl<<", d.yu="<<d.yu<<", d.zl="<<d.zl<<", d.zu="<<d.zu)
    CONCLOG_PRINTLN_AT(0,"f(x)="<<f(x)<<", g(x)="<<g(x));

    return make_tuple(f(d.x),d.x,d.yu-d.yl);
}


namespace {
template<class X> Vector<X> operator-(X const& s1, Vector<X> v2) { return esub(s1,v2); }
template<class X> Vector<X> operator*(Vector<X> v1, Vector<X> const& v2) { return emul(v1,v2); }
template<class X> Vector<X> operator/(Vector<X> v1, Vector<X> const& v2) { return ediv(v1,v2); }
}

Void
PrimalSplitDualInteriorPointOptimiser::minimisation_step(
    const ApproximateOptimisationProblem& p,
    ApproximateVectorType& w, ApproximateVectorType& x,
    ApproximateVectorType& yl, ApproximateVectorType& yu, ApproximateVectorType& zl, ApproximateVectorType& zu,
    ApproximateNumberType mu) const
{
    using XA = ApproximateNumberType;

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    CONCLOG_PRINTLN("f="<<f<<", D="<<D<<", g="<<g<<", C="<<C)
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", yl="<<yl<<", yu="<<yu<<", zl="<<zl<<", zu="<<zu<<", mu="<<mu)

    auto Cl = lower_bounds(p.C);
    auto Cu = upper_bounds(p.C);
    auto Dl = lower_bounds(p.D);
    auto Du = upper_bounds(p.D);

    auto y=yu-yl;
    auto z=zu-zl;

    auto ddfx = differential(f,x,2);
    auto ddgx = differential(g,x,2);

    auto c=ddfx.gradient();
    auto A=ddgx.jacobian();
    auto Q=ddfx.hessian();
    for (SizeType i=0; i!=y.size(); ++i) {
        Q+=y[i]*ddgx[i].hessian();
    }

// Equations (w-Cl)*yl=mu, (Cu-w)*yu=mu
// y=yu-yl
// (w-Cl)*dyl + yl * dw = (mu - (w-Cl)*yl)
// (Cu-w)*dyu - yu * dw = (mu - (Cu-w)*yu)
// (w-Cl)*(Cu-w)*dyl + (Cu-w)*yl * dw = (Cu-w)*(mu - (w-Cl)*yl)
// (w-Cl)*(Cu-w)*dyu - (w-Cl)*yu * dw = (w-Cl)*(mu - (Cu-w)*yu)
// (w-Cl)*(Cu-w)*dy  - ((w-Cl)*yu+(Cu-w)*yl) * dw = (w-Cl)*(mu - (Cu-w)*yu)-(Cu-w)*(mu - (w-Cl)*yl)
//      = (2*w-Cl-Cu)*mu - (w-Cl)*(Cu-w)*y
//
// If Cu=+inf, yu=0, so (w-Cl)*(-dyl) - yl*dw = - (mu-(w-Cl)*yl)
// If Cl=-inf, yl=0, so (Cu-w)*(+dyu) - yu*dw = + (mu-(Cu-w)*yu)
// If Cl=Cu=0, -w^2*dy - (w*yu - w*yl) * dw = w*(mu + w*yu)+w*(mu-w*yl)
//    -w*dy - (yu-yl) * dw = 2*mu + w * (yu-yl)


    auto W = DiagonalMatrix<XA>((w-Cl)*(Cu-w));
    auto X = DiagonalMatrix<XA>((x-Dl)*(Du-x));
    auto Y = DiagonalMatrix<XA>((w-Cl)*yu+(Cu-w)*yl);
    auto Z = DiagonalMatrix<XA>((x-Dl)*zu+(Du-x)*zl);

    auto rwl = mu-(w-Cl)*yl;
    auto rwu = mu-(Cu-w)*yu;
    auto rzl = mu-(x-Dl)*zl;
    auto rzu = mu-(Du-x)*zu;

    auto rw = (w-Cl)*rwu-(Cu-w)*rwl;
    auto rx = transpose(c)+transpose(A)*y+z;
    auto ry = g(x)-w;
    auto rz = (x-Dl)*rzu-(Du-x)*rzl;

    SlackPrimalDualComplementaryMatrix<XA> S(Q,A,W,X,Y,Z);
    SlackPrimalDualComplementaryData<XA> r(rw,rx,ry,rz);
    SlackPrimalDualComplementaryData<XA> d=S.solve(r);

    auto& dw=d.w; auto& dx=d.x; auto& dy=d.y; auto& dz=d.z;


    // Use relations for dw,dy (similary for dx,dz)
    // Note must have dy = dyu-dyl
    // (cu-w) * dyu - yu * dw = (mu - (cu-w)*yu)
    // (w-cl) * dyl + yl * dw = (mu - (w-cl)*yl)

    auto dyl = (rwl-yl*dw) / (w-Cl);
    auto dyu = (rwu+yu*dw) / (Cu-w);
    auto dzl = (rzl-zl*dx) / (x-Dl);
    auto dzu = (rzu+zu*dx) / (Du-x);

    ARIADNE_DEBUG_ASSERT(decide(norm(dyu-dyl-dy)<1e-8));
    ARIADNE_DEBUG_ASSERT(decide(norm(dzu-dzl-dz)<1e-8));

    XA alpha=one;

    static const ExactDouble GAMMA(1.0/1024);

    auto max_vector_step = [&](Vector<XA> const& v, Vector<XA> const& dv) {
        XA m(1,v.zero_element().precision()); XA e = mu*GAMMA; for (SizeType i=0; i!=v.size(); ++i) { if (decide(dv[i]>0)) { m=min((v[i]-e)/dv[i],m); } } return m; };

    alpha = min( one, min( min(max_vector_step(yl,dyl),max_vector_step(yu,dyu)), min(max_vector_step(zl,dzl),max_vector_step(zu,dzu)) ) );

    w += alpha * dw;
    x += alpha * dx;
    yl += alpha * dyl;
    yu += alpha * dyu;
    zl += alpha * dzl;
    zu += alpha * dzu;
}


//------- InfeasibleInteriorPointOptimiser -------------------------//

InfeasibleInteriorPointOptimiser::
InfeasibleInteriorPointOptimiser() {
}

auto InfeasibleInteriorPointOptimiser::
clone() const -> InfeasibleInteriorPointOptimiser* {
    return new InfeasibleInteriorPointOptimiser(*this);
}

OutputStream& operator<<(OutputStream& os, InfeasibleInteriorPointOptimiser const& opt) {
    return os << "InfeasibleInteriorPointOptimiser()";
}



auto InfeasibleInteriorPointOptimiser::
minimise(ApproximateOptimisationProblem p) const -> ApproximateVector
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    static const ApproximateDouble VALUE_TOLERANCE=1e-8;
    static const ApproximateDouble STATE_TOLERANCE=1e-8;
    static const CounterType MAXIMUM_STEPS=24;

    StepData v(g.result_size(),g.argument_size(),dp);
    ApproximateVectorType& x=v.x;
    ApproximateVectorType& y=v.y;

    this->setup_feasibility(p,v);

    ApproximateVectorType oldx=x;

    static const ApproximateDouble MU_MIN = 1e-12;

    // FIXME: Allow more steps
    for(SizeType i=0; i!=MAXIMUM_STEPS; ++i) {
        CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        oldx=x;
        FloatDPApproximation oldfx=f(oldx);
        this->minimisation_step(p,v);
        FloatDPApproximation fx=f(x);
        if(probably(mag(fx-oldfx)<VALUE_TOLERANCE) && probably(norm(oldx-x)<STATE_TOLERANCE)) {
            break;
        }
        if(probably(v.mu<MU_MIN)) {
            break;
        }
    }
    CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));

    if (probably(D.contains(x)) && probably(C.contains(g(x)))) {
        CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        return x;
    }
    CONCLOG_PRINTLN("indeterminate_feasibility")
    throw IndeterminateFeasibilityException();
}

auto InfeasibleInteriorPointOptimiser::
minimise(ValidatedOptimisationProblem p) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    static const ApproximateDouble VALUE_TOLERANCE=1e-8;
    static const ApproximateDouble STATE_TOLERANCE=1e-8;
    static const CounterType MAXIMUM_STEPS=24;

    ARIADNE_ASSERT(f.argument_size()==D.size());
    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());
    StepData v;
    ApproximateVectorType& x=cast_approximate(v.x);
    ApproximateVectorType& y=cast_approximate(v.y);
    C=intersection(C,cast_exact_box(apply(g,D)+UpperIntervalVectorType(C.size(),UpperIntervalType(-1,+1))));
    this->setup_feasibility(p,v);
    ApproximateVectorType oldx=x;

    static const ExactDouble MU_MIN = 1e-12_pr;

    // FIXME: Allow more steps
    for(SizeType i=0; i!=MAXIMUM_STEPS; ++i) {
        CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        oldx=x;
        FloatDPApproximation oldfx=f(oldx);
        this->minimisation_step(p,v);
        if(FeasibilityChecker().validate_infeasibility(p,cast_exact(y))) {
            CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            CONCLOG_PRINTLN_AT(1,"Infeasible");
            std::cerr<<"EXCEPTION: "<<InfeasibleProblemException().what()<<"\n";
            throw InfeasibleProblemException();
        }
        FloatDPApproximation fx=f(x);
        if(probably(mag(fx-oldfx)<VALUE_TOLERANCE) && probably(norm(oldx-x)<STATE_TOLERANCE)) {
            break;
        }
        if(v.mu.raw()<MU_MIN) {
            break;
        }
    }
    CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));

    if(FeasibilityChecker().validate_feasibility(p,cast_exact(x))) {
        CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        return to_generic(cast_exact(x));
    }
    CONCLOG_PRINTLN("indeterminate_feasibility");
    throw IndeterminateFeasibilityException();
}

auto InfeasibleInteriorPointOptimiser::
feasible(ValidatedFeasibilityProblem p) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE
    CONCLOG_PRINTLN("p="<<p);

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());

    StepData v;
    ApproximateVectorType& x=cast_approximate(v.x);
    ApproximateVectorType& y=cast_approximate(v.y);

    ApproximateScalarMultivariateFunction f(EuclideanDomain(D.dimension()));
    auto R=intersection(cast_exact_box(widen(apply(g,D),1)),C);

    ApproximateOptimisationProblem optp(f,D,g,R);
    this->setup_feasibility(p,v);

    static const ExactDouble MU_MIN = 1e-12_pr;

    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        this->minimisation_step(optp,v);
        if(FeasibilityChecker().validate_feasibility(p,cast_exact(x))) {
            CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            CONCLOG_PRINTLN("Feasible");
            return true;
        }
        if(FeasibilityChecker().validate_infeasibility(p,cast_exact(y))) {
            CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            CONCLOG_PRINTLN("Infeasible");
            return false;
        }
        if(v.mu.raw()<MU_MIN) {
            break;
        }
    }
    CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    CONCLOG_PRINTLN("Indeterminate");
    return indeterminate;
}


auto InfeasibleInteriorPointOptimiser::
feasible(ApproximateFeasibilityProblem p) const -> ApproximateKleenean
{
    CONCLOG_SCOPE_CREATE
    CONCLOG_PRINTLN("p="<<p);

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());

    StepData v;
    ApproximateVectorType& x=cast_approximate(v.x);
    ApproximateVectorType& y=cast_approximate(v.y);

    ApproximateScalarMultivariateFunction f(EuclideanDomain(D.dimension()));

    ApproximateOptimisationProblem optp(f,D,g,C);
    this->setup_feasibility(p,v);

    static const ExactDouble MU_MIN = 1e-12_pr;

    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        this->minimisation_step(optp,v);
        if(v.mu.raw()<MU_MIN) {
            break;
        }
    }
    CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    return C.contains(g(x));
}


auto InfeasibleInteriorPointOptimiser::
feasible_hotstarted(ValidatedFeasibilityProblem p,
                    const PrimalDualData& wxy0) const
                        -> Pair<ValidatedKleenean,ApproximateVectorType>
{
    ARIADNE_NOT_IMPLEMENTED;
}


Void InfeasibleInteriorPointOptimiser::
setup_feasibility(const ApproximateFeasibilityProblem& p,
                  StepData& v) const
{
    ExactIntervalType I(-1,+1);
    SizeType m=p.C.size(); SizeType n=p.D.size();

    v.x=midpoint(p.D);
    v.y=ApproximateVectorType(m,zero);
    v.w=midpoint(p.C);

    //stp.xl=lower_bound(D)-x;
    v.wl=Vector(m,-one);
    v.wu=Vector(m,+one);
    v.xl=lower_bounds(p.D)-v.x;
    v.xu=upper_bounds(p.D)-v.x;
    v.vl=Vector(m,-one);
    v.vu=Vector(m,+one);
    v.zl=Vector(n,-one);
    v.zu=Vector(n,+one);
    // FIXME: What should relaxation parameter be?
    v.mu=1.0_x;
}


Void
InfeasibleInteriorPointOptimiser::minimisation_step(
    const ApproximateOptimisationProblem& p,
    StepData& v) const
{
    auto& f=p.f; auto& d=p.D; auto& g=p.g; auto& c=p.C;

    ApproximateVectorType& w=v.w; ApproximateVectorType& x=v.x; ApproximateVectorType& y=v.y;
    FloatDPApproximation& mu=v.mu;
    ApproximateVectorType& wl=v.wl; ApproximateVectorType& wu=v.wu;
    ApproximateVectorType& xl=v.xl; ApproximateVectorType& xu=v.xu;
    ApproximateVectorType& vl=v.vl; ApproximateVectorType& vu=v.vu;
    ApproximateVectorType& zl=v.zl; ApproximateVectorType& zu=v.zu;
    ApproximateVectorType cl=lower_bounds(c); ApproximateVectorType cu=upper_bounds(c);
    ApproximateVectorType dl=lower_bounds(d); ApproximateVectorType du=upper_bounds(d);

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    CONCLOG_PRINTLN("w ="<<w<<",  x ="<<x<<", y ="<<y<<" mu="<<mu);
    CONCLOG_PRINTLN("wl="<<wl<<", wu="<<wu<<", xl="<<xl<<", xu="<<xu);
    CONCLOG_PRINTLN("vl="<<vl<<", vu="<<vu<<", zl="<<zl<<", zu="<<zu);
    CONCLOG_PRINTLN("cl-wl="<<cl-wl<<", dl-xl="<<dl-xl);
    CONCLOG_PRINTLN("  w  ="<<w<<",   x  ="<<x);
    CONCLOG_PRINTLN("cu-wu="<<cu-wu<<", du-xu="<<du-xu);
    static const ExactDouble gamma=0.0009765625_x;
    static const ExactDouble sigma=0.125_x;
    static const ExactDouble scale=0.75_x;

    const SizeType n=d.size();
    const SizeType m=c.size();

    ARIADNE_ASSERT_MSG(f.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.result_size()==c.size(),  "f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT(w.size()==m);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    mu = mu * sigma;

    ApproximateVectorType ax(x);
    FloatDPApproximationDifferential ddfx=f.differential(ax,2u);
    CONCLOG_PRINTLN("ddfx="<<ddfx);
    Vector<FloatDPApproximationDifferential> ddgx=g.differential(ax,2u);
    CONCLOG_PRINTLN("ddgx="<<ddgx);

    FloatDPApproximation fx = ddfx.value();
    Vector<FloatDPApproximation> gx = ddgx.value();
    CONCLOG_PRINTLN("f(x)="<<fx);
    CONCLOG_PRINTLN("g(x)="<<gx);
    Vector<FloatDPApproximation> Jfx = transpose(ddfx.gradient());
    Matrix<FloatDPApproximation> A = ddgx.jacobian();
    Matrix<FloatDPApproximation>& Jgx = A;
    CONCLOG_PRINTLN("Df(x)="<<Jfx);
    CONCLOG_PRINTLN("Dg(x)="<<Jgx);

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<FloatDPApproximation> YH = ddfx.hessian();
    for(SizeType i=0; i!=m; ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    CONCLOG_PRINTLN("D2f(x)="<<ddfx.hessian());
    CONCLOG_PRINTLN("D2f(x)+Y.D2g(x)="<<YH);

    // Set up the system of equations
    // (A^TDA + E - Y.H) dx = A^T(r_w-Dr_y)+r_x
    // dw = A \delta x + r_y
    // dy = r_w - D dw

    FloatDPApproximationDiagonalMatrix const& Vl=diagonal_matrix(vl);
    FloatDPApproximationDiagonalMatrix const& Vu=diagonal_matrix(vu);
    FloatDPApproximationDiagonalMatrix const& Wl=diagonal_matrix(wl);
    FloatDPApproximationDiagonalMatrix const& Wu=diagonal_matrix(wu);
    FloatDPApproximationDiagonalMatrix const& Xl=diagonal_matrix(xl);
    FloatDPApproximationDiagonalMatrix const& Xu=diagonal_matrix(xu);
    FloatDPApproximationDiagonalMatrix const& Zl=diagonal_matrix(zl);
    FloatDPApproximationDiagonalMatrix const& Zu=diagonal_matrix(zu);

    // Compute the diagonal matrices
    //   D=XL/ZL+XU/ZU  E=WL/VL+WU/VU
    FloatDPApproximationDiagonalMatrix Dl=Vl/Wl;
    FloatDPApproximationDiagonalMatrix Du=Vu/Wu;
    FloatDPApproximationDiagonalMatrix D=Dl+Du;
    CONCLOG_PRINTLN("D="<<D);
    FloatDPApproximationDiagonalMatrix El=Zl/Xl;
    FloatDPApproximationDiagonalMatrix Eu=Zu/Xu;
    FloatDPApproximationDiagonalMatrix E=El+Eu;
    CONCLOG_PRINTLN("E="<<E);

    // Construct the residuals
    // The residual for the slack variable xl is given by the duality condition xl.zl=mu as mu/xl-zl
    // The residual for the dual variable zl is given by the slackness condition x-xl-cl
    // The residual for the auxiliary variable w is given by y-(vu-vl)
    // The residual for the dual variable y is given by g(x)-w
    ApproximateVectorType ew=(vl+vu)-y;
    ApproximateVectorType ex=Jfx+transpose(Jgx)*y+(zl+zu);
    ApproximateVectorType ey=gx-w;
    ApproximateVectorType ewl=esub(vl,ediv(mu,wl));
    ApproximateVectorType ewu=esub(vu,ediv(mu,wu));
    ApproximateVectorType exl=esub(zl,ediv(mu,xl));
    ApproximateVectorType exu=esub(zu,ediv(mu,xu));
    ApproximateVectorType evl=w+wl-cl;
    ApproximateVectorType evu=w+wu-cu;
    ApproximateVectorType ezl=x+xl-dl;
    ApproximateVectorType ezu=x+xu-du;

    CONCLOG_PRINTLN("ew="<<ew<<", ex="<<ex<<", ey="<<ey);
    CONCLOG_PRINTLN("ewl="<<ewl<<", ewu="<<ewu<<", exl="<<exl<<" exu="<<exu);
    CONCLOG_PRINTLN("evl="<<evl<<", evu="<<evu<<", ezl="<<ezl<<" ezu="<<ezu);

    ApproximateVectorType rw = ew - (ewl+ewu) + Dl*evl + Du*evu;
    ApproximateVectorType rx = ex - (exl+exu) + El*ezl + Eu*ezu;
    ApproximateVectorType& ry = ey;

    // Solve linear system
    // ( D   0  -I ) (dw)   (rw)
    // ( 0  H+E A^T) (dx) = (rx)
    // (-I   A   0 ) (dy) = (ry)

    // Eliminate dw=A*dx-ry and dy=D*dw-rw
    // Reduce to (H+E)*dx+A'*(D*(A*dx-ry)-rw)=rx
    // Simplify to (H+E+A'*D*A)*dx=rx+A'*(D*ry+rw)

    // normal equation matrix
    FloatDPApproximationMatrix S=YH;
    atda(S,A,D);
    S+=E;

    CONCLOG_PRINTLN("S="<<S);
    ARIADNE_DEBUG_ASSERT(decide(norm(FloatDPApproximationMatrix(S-(YH+E+transpose(A)*(D*A))))/norm(S)<1e-8));
    CONCLOG_PRINTLN("Sinv="<<inverse(S));

    ApproximateVectorType r = transpose(A)*(rw+D*ry)+rx;
    CONCLOG_PRINTLN("rw="<<rw<<" rx="<<rx<<" ry="<<ry);
    CONCLOG_PRINTLN("r="<<r);

    // Compute the differences
    ApproximateVectorType dx = solve(S,r);
    // Apply correction to improve accuracy
    dx = dx + solve(S,r-S*dx);
    CONCLOG_PRINTLN("S*dx="<<S*dx<<" r="<<r);
    CONCLOG_PRINTLN("S*inverse(S)-I="<<S*inverse(S)-FloatDPApproximationMatrix::identity(n,dp));
    ARIADNE_DEBUG_ASSERT(decide(norm(S*dx - r)/max(1.0_x,norm(r))<1e-4));

    ApproximateVectorType dw = A*dx-ry;
    ApproximateVectorType dy = D*dw-rw;
    CONCLOG_PRINTLN("dw="<<dw<<" dx="<<dx<<" dy="<<dy);

    CONCLOG_PRINTLN("YH*dx+E*dx+dy*A="<<(YH*dx+E*dx+transpose(A)*dy)<<", rx="<<rx);

    // Check solution of linear system for residuals
    ARIADNE_DEBUG_ASSERT(decide(norm(D*dw-dy-rw)/max(one,norm(rw))<1e-4));
    ARIADNE_DEBUG_ASSERT(decide(norm(YH*dx+E*dx+transpose(A)*dy-rx)/max(one,norm(rx))<1e-2));
    ARIADNE_DEBUG_ASSERT(decide(norm(-dw+A*dx-ry)/max(one,norm(ry))<1e-4));


    ApproximateVectorType dwl = evl-dw;
    ApproximateVectorType dwu = evu-dw;
    ApproximateVectorType dxl = ezl-dx;
    ApproximateVectorType dxu = ezu-dx;
    ApproximateVectorType dvl = ewl-Dl*dwl;
    ApproximateVectorType dvu = ewu-Du*dwu;
    ApproximateVectorType dzl = exl-El*dxl;
    ApproximateVectorType dzu = exu-Eu*dxu;

    CONCLOG_PRINTLN("dwl="<<dwl<<", dwu="<<dwu<<", dxl="<<dxl<<" dxu="<<dxu);
    CONCLOG_PRINTLN("dvl="<<dvl<<", dvu="<<dvu<<", dzl="<<dzl<<" dzu="<<dzu);

    CONCLOG_PRINTLN("YH*dx+dy*A+dzl+dzu="<<(YH*dx+transpose(A)*dy+dzl+dzu)<<", ex="<<ex);
    // Check solution of linear system
/*
    ARIADNE_DEBUG_ASSERT(norm(-dy+dvl+dvu - ew)/max(1.0,norm(ew))<1e-4);
    ARIADNE_DEBUG_ASSERT(norm(YH*dx+dy*A+dzl+dzu - ex)/max(1.0,norm(ex))<1e-2);
    ARIADNE_DEBUG_ASSERT(norm(-dw+A*dx - ey)/max(1.0,norm(ey))<1e-4);
    ARIADNE_DEBUG_ASSERT(norm(Dl*dwl+dvl - ewl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(Du*dwu+dvu - ewu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(El*dxl+dzl - exl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(Eu*dxu+dzu - exu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dw+dwl - evl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dw+dwu - evu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dx+dxl - ezl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dx+dxu - ezu)<1e-12);
*/

    ApproximateVectorType nw; ApproximateVectorType nx; ApproximateVectorType ny;
    ApproximateVectorType nwl; ApproximateVectorType nwu; ApproximateVectorType nxl; ApproximateVectorType nxu;
    ApproximateVectorType nvl; ApproximateVectorType nvu; ApproximateVectorType nzl; ApproximateVectorType nzu;


    FloatDPApproximation alpha=one;
    nx = x-alpha*dx;
    // Pick an update value which minimises the objective function
    FloatDPApproximation fxmin=f(nx);
    FloatDPApproximation alphamin=one;
    static const CounterType REDUCTION_STEPS=4;
    for(SizeType i=0; i!=REDUCTION_STEPS; ++i) {
        alpha*=scale;
        nx = x-alpha*dx;
        FloatDPApproximation fnx=f(nx);
        if(decide(fnx<fxmin*scale)) {
            fxmin=fnx;
            alphamin=alpha;
        }
    }
    //alpha=alphamin;
    alpha=one;

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly.
    static const double MINIMUM_ALPHA=1e-16;
    Bool allfeasible=false;
    while(decide(alpha>MINIMUM_ALPHA) && !allfeasible) {
        nwl=wl-alpha*dwl;
        nwu=wu-alpha*dwu;
        nxl=xl-alpha*dxl;
        nxu=xu-alpha*dxu;
        nvl=vl-alpha*dvl;
        nvu=vu-alpha*dvu;
        nzl=zl-alpha*dzl;
        nzu=zu-alpha*dzu;
        allfeasible = elss(nwl,mu*gamma) && egtr(nwu,mu*gamma) && elss(nxl,mu*gamma) && egtr(nxu,mu*gamma)
                           && elss(nvl,mu*gamma)&& egtr(nvu,mu*gamma) && elss(nzl,mu*gamma) && egtr(nzu,mu*gamma);
        //allfeasible = eneg(nwl) && epos(nwu) && eneg(nxl) && epos(nxu) && elss(nvl,mu*gamma) && egtr(nvu,mu*gamma) && elss(nzl,mu*gamma) && egtr(nzu,mu*gamma);
        if(!allfeasible) { alpha*=scale; }
    }
    nw=w-alpha*dw;
    nx=x-alpha*dx;
    ny=y-alpha*dy;
    if(decide(alpha<=MINIMUM_ALPHA)) {
        CONCLOG_PRINTLN_AT(1,"w="<<w<<"  x="<<x<<"  y="<<y);
        CONCLOG_PRINTLN_AT(1,"nw="<<nw<<"  nx="<<nx<<"  ny="<<ny);
        throw NearBoundaryOfFeasibleDomainException(); }
    CONCLOG_PRINTLN("alpha="<<alpha);
    CONCLOG_PRINTLN("nw="<<nw<<" nx="<<nx<<" ny="<<ny);
    CONCLOG_PRINTLN("nwl="<<nwl<<", nwu="<<nwu<<", nxl="<<nxl<<" nxu="<<nxu);
    CONCLOG_PRINTLN("nvl="<<nvl<<", nvu="<<nvu<<", nzl="<<nzl<<" nzu="<<nzu);

    w=nw; x=nx; y=ny;
    wl=nwl; wu=nwu; xl=nxl; xu=nxu;
    vl=nvl; vu=nvu; zl=nzl; zu=nzu;

    FloatDPApproximation nmu = zero;
    for(SizeType i=0; i!=m; ++i) {
        nmu = nmu + wl[i]*vl[i] + wu[i]*vu[i];
    }
    for(SizeType j=0; j!=n; ++j) {
        nmu = nmu + xl[j]*zl[j] + xu[j]*zu[j];
    }
    nmu /= (2*(m+n));
    mu = nmu;

    CONCLOG_PRINTLN("nmu="<<nmu);

}





/*

// Solve max log(x-xl) + log(xu-x) + log(zu-z) + log(z-zl) such that g(x)=z
//   if zl[i]=zu[i] then z=zc is hard constraint
//   alternatively, use the penalty (z[j]-zc[j])^2/2 instead

// KKT conditions
//     1/(x-xl) - 1/(xu-x) + y Dg(x) = 0
//     1/(z-zl) - 1/(zu-z) - y = 0
//     g(x) - z = 0
//   If zu[j]=inf, then 1/(z[j-zl[j]) - y[j] = 0, so y[j]>=0
//   If zl[j]=zu[j], then use the equation z[j]=zc[j] instead
//
// Re-write KKT conditions for x,z as
//     (xu-xl) + y g(x) (x-xl)(xu-x) = 0
//     (zu-zl) - y(z-zl)(zu-z) = 0
//  Or 1 - y(z-zl)(zu-z)/(zu-zl) = 0
//   If zu[j] = inf, then 1 - y(z-zl) = 0
//      zl[j]=zu[j]=zc[j], then y(z-zc)^2 = 0
//
// Derivative matrix
//     - 1/(x-xl)^2 - 1/(xu-x)^2 + y D^2g = 0
//
// PROBLEM:
//   Dg can be singular, even at intermediate points
//
// FJ conditions
//     mu/(x-xl) - mu/(xu-x) + y D g(x) = 0
//     mu/(z-zl) - mu/(zu-z) - y = 0
//     g(x) - z = 0
//     sum y^2 - mu = 0
//   If zu[j]=inf, then mu/(z[j-zl[j]) - y[j] = 0, so y[j]>=0
//   If zl[j]=zu[j], then -(z-zc)/mu - y = 0 instead
//
// Re-write FJ conditions for x,z as
// Derivative matrix
//     - mu/(x-xl)^2 dx - mu/(xu-x)^2 dx + y D^2g dx + Dg^T dy + (1/(x-xl) - (1/xu-x)) dmu
//     - mu/(z-zl)^2 dz - mu/(zu-z)^2 dz + I dy + (1/(z-zl) - (1/zu-z)) dmu
//     Dg dx - I dz
//    2y . dy - dmu
Void PenaltyFunctionOptimiser::
feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                 RawFloatDPVector& x, RawFloatDPVector& y, RawFloatDPVector& z) const
{
    CONCLOG_PRINTLN_AT(2,"feasibility_step");
    RawFloatDPVector xl=lower_bounds(D); RawFloatDPVector xu=upper_bounds(D);
    RawFloatDPVector zl=lower_bounds(C); RawFloatDPVector zu=upper_bounds(C);

    const SizeType n=x.size();
    const SizeType m=y.size();

    CONCLOG_PRINTLN_AT(4,"x="<<x<<" y="<<y<<" z="<<z);
    Vector<FloatDPDifferential> ddx = FloatDPDifferential::variables(2,x);
    Vector<FloatDPDifferential> ddgx = g.evaluate(ddx);

    FloatDPMatrix A=ddgx.jacobian();
    CONCLOG_PRINTLN_AT(6,"A="<<A);
    RawFloatDPVector v = join(join(x,z),y);

    RawFloatDPVector r(n+2*m,n+2*m);
    project(r,range(0,n)) = y * A;
    for(SizeType i=0; i!=n; ++i) {
        r[i] += ( rec(x[i]-xl[i]) - rec(xu[i]-x[i]) );
    }
    for(SizeType j=0; j!=m; ++j) {
        if(zl[j]==zu[j]) { assert(zu[j]==zl[j]); r[n+j] = z[j]-zl[j]; }
        else { r[n+j] = ( rec(z[j]-zl[j]) - rec(zu[j]-z[j]) - y[j] ); }
    }
    project(r,range(n+m,n+2*m)) = ddgx.value() - z;
    r[n+2*m]=0.0;
    CONCLOG_PRINTLN_AT(5,"r="<<r);

    FloatDPMatrix S(n+2*m+1,n+2*m+1);
    for(SizeType j=0; j!=m; ++j) {
        FloatDPMatrix H=ddgx[j].hessian();
        for(SizeType i1=0; i1!=n; ++i1) {
            for(SizeType i2=0; i2!=n; ++i2) {
                S[i1][i2]+=y[j]*H[i1][i2];
            }
        }
    }
    for(SizeType j=0; j!=m; ++j) {
        for(SizeType i=0; i!=n; ++i) {
            S[i][j+m+n]=A[j][i];
            S[j+m+n][i]=A[j][i];
        }
    }
    for(SizeType j=0; j!=m; ++j) {
        S[n+j][n+m+j] = -1.0;
        S[n+m+j][n+j] = -1.0;
        //if(zl[j]==zu[j]) { S[n+j][n+j] = -1.0; S[n+j][n+m+j] = 0.0; }
        if(zl[j]==zu[j]) { S[n+j][n+j] = +inf; }
        else { S[n+j][n+j] = - rec(sqr(z[j]-zl[j])) - rec(sqr(zu[j]-z[j])); }
    }
    for(SizeType i=0; i!=n; ++i) {
        S[i][i]-= rec(sqr(xu[i]-x[i]));
        S[i][i]-= rec(sqr(x[i]-xl[i]));
    }

    for(SizeType i=0; i!=n; ++i) {
        S[i][n+2*m] -= rec(xu[i]-x[i]);
        S[i][n+2*m] += rec(x[i]-xl[i]);
    }

    for(SizeType j=0; j!=n; ++j) {
        //S[n+m+j][n+m+j] = -1.0/1024;
    }

    CONCLOG_PRINTLN_AT(5,"S="<<S);
    //CONCLOG_PRINTLN_AT(5,"S="<<std::fixed<<pretty(S));

    FloatDPMatrix Sinv = inverse(S);
    //CONCLOG_PRINTLN_AT(9,"Sinv="<<Sinv);
    //CONCLOG_PRINTLN_AT(5,"Sinv="<<std::fixed<<pretty(Sinv));

    RawFloatDPVector dv = Sinv * r;
    CONCLOG_PRINTLN_AT(5,"dv="<<dv);

    FloatDP alpha = 1.0;
    RawFloatDPVector nv = v-dv;
    while(!contains(D,RawFloatDPVector(project(nv,range(0,n)))) || !contains(C,RawFloatDPVector(project(nv,range(n,n+m)))) ) {
        alpha *= 0.75;
        nv = v-alpha*dv;
    }

    CONCLOG_PRINTLN_AT(4,"nv="<<nv<<" a="<<alpha);

    x=project(nv,range(0,n));
    z=project(nv,range(n,n+m));
    y=project(nv,range(n+m,n+2*m));
    CONCLOG_PRINTLN_AT(4,"g(x)-z="<<g(x)-z);

}
*/

//------- IntervalOptimiser -----------------------------------//

IntervalOptimiser* IntervalOptimiser::clone() const {
    return new IntervalOptimiser(*this);
}


// Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
ValidatedKleenean IntervalOptimiser::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", h="<<h);

    const SizeType n=D.size();

    FloatDPBoundsVector zl(n,dp), zu(n,dp);
    ExactFloatDPVectorType xl = Ariadne::lower_bounds(D);
    ExactFloatDPVectorType xu = Ariadne::upper_bounds(D);

    FloatDPBoundsVector x=cast_singleton(D);
    FloatDPBoundsVector y(h.result_size(),FloatDPBounds(-1,+1,dp));
    FloatDPBounds mu(0,1,dp);

    for(SizeType i=0; i!=8; ++i) {
        this->feasibility_step(xl,xu,h,x,y,zl,zu,mu);
    }

    return indeterminate;
}

Void IntervalOptimiser::
feasibility_step(const ExactFloatDPVectorType& xl, const ExactFloatDPVectorType& xu, const ValidatedVectorMultivariateFunction& h,
                 FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("[x]="<<x<<" [lambda]="<<y<<", [zl]="<<zl<<", [zu]="<<zu<<" [mu]="<<mu);

    const SizeType n=x.size();
    const SizeType m=y.size();

    FloatDPBoundsVector mx=midpoint(x);
    FloatDPBoundsVector my=midpoint(y);
    FloatDPBoundsVector mzl=midpoint(zl);
    FloatDPBoundsVector mzu=midpoint(zu);
    FloatDPBounds mmu(midpoint(mu));
    CONCLOG_PRINTLN_AT(1,"x~"<<x<<" lambda~="<<y<<", mu~"<<mu);

    // Solve equations y Dh(x) - zl + zu = 0; h(x) = 0; (x-xl).zl - mu = 0;  (xu-x).zu - mu = 0; Sum_j y_j^2 - mu = 0
    Vector<FloatDPBoundsDifferential> ddhx=h.evaluate(FloatDPBoundsDifferential::variables(2,x));
    Vector<FloatDPBoundsDifferential> dhmx=h.evaluate(FloatDPBoundsDifferential::variables(1,mx));
    FloatDPBoundsMatrix A = ddhx.jacobian();
    FloatDPBoundsMatrix mA = dhmx.jacobian();
    CONCLOG_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPBoundsVector rx = transpose(mA) * my;
    for(SizeType j=0; j!=n; ++j) {
        rx[j] -= mmu*rec(mx[j]-xl[j]);
        rx[j] += mmu*rec(xu[j]-mx[j]);
    }
    FloatDPBoundsVector ry = dhmx.value();
    FloatDPBoundsVector rzl = esub(emul(FloatDPBoundsVector(mx-xl),mzl),mmu);
    FloatDPBoundsVector rzu = esub(emul(FloatDPBoundsVector(xu-mx),mzu),mmu);
    CONCLOG_PRINTLN("rx="<<rx<<" ry="<<ry<<" rzl="<<rzl<<" rzu="<<rzu);

    FloatDPBoundsMatrix H(n,n,dp);
    for(SizeType i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(SizeType j=0; j!=n; ++j) {
        H[j][j] += mu*rec(sqr(x[j]-xl[j]));
        H[j][j] += mu*rec(sqr(xu[j]-x[j]));
    }

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPBoundsMatrix Hinv=inverse(H);
    CONCLOG_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPBoundsMatrix S=A*Hinv*transpose(A);
    FloatDPBoundsMatrix Sinv=inverse(S);
    CONCLOG_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    FloatDPBoundsVector dy = Sinv * ( A*(Hinv*rx) - ry );
    FloatDPBoundsVector dx = Hinv * ( rx - transpose(A) * dy);
    CONCLOG_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPBoundsVector nx = x-dx;
    FloatDPBoundsVector ny = y-dy;
    CONCLOG_PRINTLN("nx="<<nx<<" ny="<<ny);
    CONCLOG_PRINTLN_AT(1,"h(x)="<<h(nx));

    x = refinement(x,nx); y=refinement(y,ny);
    FloatDPBounds nmu = zero;
    for(SizeType i=0; i!=m; ++i) { nmu += sqr(y[i]); }
    mu=refinement(mu,nmu);
}


//------- Optimality condition functions -----------------------------------//

/*

struct KuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<KuhnTuckerFunctionBody,ExactIntervalType>
{
    ValidatedScalarMultivariateFunction f;
    Array<ValidatedScalarMultivariateFunction> g;
    Array<ValidatedScalarMultivariateFunction> df;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    KuhnTuckerFunctionBody(ValidatedScalarMultivariateFunction _f, ValidatedVectorMultivariateFunction _g) {
        ARIADNE_ASSERT(_f.argument_size()==_g.argument_size());
        const SizeType m=_g.argument_size();
        const SizeType n=_g.result_size();
        g.resize(n); df.resize(m); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        f=_f;
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; }
        for(SizeType i=0; i!=m; ++i) { df[i]=f.derivative(i); }
        for(SizeType j=0; j!=n; ++j) { for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return g.size()*2+f.argument_size(); }
    SizeType argument_size() const { return g.size()*2+f.argument_size(); }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const SizeType m=f.argument_size();
        const SizeType n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        Vector<X> rx(m), rz(n), rs(n);
        for(SizeType i=0; i!=m; ++i) { rx[i]=df[i].evaluate(y); for(SizeType j=0; j!=n; ++j) { rx[i]=rx[i]-x[j]*dg[j][i].evaluate(y); } }
        for(SizeType j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + z[j]; }
        for(SizeType j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
    }
};

struct FeasibilityKuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<FeasibilityKuhnTuckerFunctionBody,ExactIntervalType>
{
    Array<ValidatedScalarMultivariateFunction> g;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    FeasibilityKuhnTuckerFunctionBody(ValidatedVectorMultivariateFunction _g) {
        const SizeType m=_g.argument_size();
        const SizeType n=_g.result_size();
        g.resize(n); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return g.size()*2+g[0].argument_size()+1; }
    SizeType argument_size() const { return g.size()*2+g[0].argument_size()+1; }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const SizeType m=g[0].argument_size();
        const SizeType n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        X t(arg[n+m+n]);
        Vector<X> rx(m), rz(n), rs(n); X rt;
        for(SizeType i=0; i!=m; ++i) { rx[i]=x[0]*dg[0][i].evaluate(y); for(SizeType j=1; j!=n; ++j) { rx[i]=rx[i]+x[j]*dg[j][i].evaluate(y); } }
        for(SizeType j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + t + z[j]; }
        for(SizeType j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        rt=1-x[0]; for(SizeType j=1; j!=n; ++j) { rt=rt-x[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
        res[n+m+n]=rt;
    }
};



struct ConstrainedFeasibilityKuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<FeasibilityKuhnTuckerFunctionBody,ExactIntervalType>
{
    SizeType m;
    SizeType n;
    ExactIntervalVectorType d;
    Array<ValidatedScalarMultivariateFunction> g;
    ExactIntervalVectorType c;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    ConstrainedFeasibilityKuhnTuckerFunctionBody(ExactBoxType D, ValidatedVectorMultivariateFunction _g, ExactBoxType C) {
        m=_g.argument_size();
        n=_g.result_size();
        d=D; c=C;
        g.resize(n); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return 5*m+4*n+1u; }
    SizeType argument_size() const { return 5*m+4*n+1u; }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream& os) const { return os << "KuhnTuckerFunctionBody"; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const X zero=arg[0].zero_element();
        const SizeType l=2*(m+n);
        assert(arg.size()==l+m+l+1);
        Vector<X> x(project(arg,range(0u,l)));
        Vector<X> y(project(arg,range(l,l+m)));
        Vector<X> z(project(arg,range(l+m,l+m+l)));
        X t(arg[l+m+l]);
        Vector<X> rx(m,zero), rz(l,zero), rs(l,zero); X rt(zero);
        Vector<X> gy(n);
        for(SizeType j=0; j!=n; ++j) { gy[j]=g[j].evaluate(y); }
        Matrix<X> dgy(n,m);
        for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { dgy[j][i]=dg[j][i].evaluate(y); } }

        for(SizeType i=0; i!=m; ++i) {
            for(SizeType j=0; j!=n; ++j) { rx[i]+=x[j]*(dgy[j][i]-c[j].upper_bound()); rx[i]+=x[n+j]*(c[j].lower_bound()-dgy[j][i]); }
            rx[i]+=x[2*n+i]*(y[i]-d[i].upper_bound())-x[2*n+m+i]*(d[i].lower_bound()-y[i]);
        }
        for(SizeType j=0; j!=n; ++j) { rz[j]=gy[j] + t + z[j]; rz[n+j]=t+z[n+j]-gy[j]; }
        for(SizeType i=0; i!=m; ++i) { rz[2*n+i]=y[i]+t+z[2*n+i]; rz[2*n+m+i]=y[i]+t+z[2*n+m+i]; }
        for(SizeType k=0; k!=l; ++k) { rs[k]=x[k]*z[k]; }
        rt+=1.0; for(SizeType j=0; j!=2*n; ++j) { rt=rt-x[j]; }
        project(res,range(0,l))=rz;
        project(res,range(l,l+m))=rx;
        project(res,range(l+m,l+m+l))=rs;
        res[l+m+l]=rt;
    }
};

*/


//------- KrawczykOptimiser -----------------------------------//

/*

ValidatedVector KrawczykOptimiser::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


ValidatedKleenean KrawczykOptimiser::
feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    CONCLOG_PRINTLN("KrawczykOptimiser::feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c)");
    CONCLOG_PRINTLN("d="<<d<<", g="<<g<<", c="<<c);

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());

    ExactIntervalType t; ExactIntervalVectorType x,y,z;
    setup_feasibility(d,g,c,x,y,z,t);


    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z);
        try {
            this->feasibility_step(d,g,c,x,y,z,t);
        }
        catch(const SingularMatrixException& e) {
            return indeterminate;
        }
        if(t.lower_bound()>t.upper_bound()) {
            CONCLOG_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return indeterminate;
        }
        if(t.lower_bound()>0.0) {
            CONCLOG_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            if(this->is_feasible_point(d,g,c,midpoint(y))) {
                return true;
            }
        }
        if(t.upper_bound()<0.0) {
            CONCLOG_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return false;
        }
    }
    CONCLOG_PRINTLN_AT(1,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
    if(this->is_infeasibility(d,g,c,midpoint(x))) {
        return false;
    }
    return indeterminate;
}



Void KrawczykOptimiser::setup_feasibility(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                          ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    const SizeType m=g.argument_size();
    const SizeType n=g.result_size();
    const SizeType l=2*(m+n);
    x=ExactIntervalVectorType(l, ExactIntervalType(0,1)/l);
    y=d;
    z.resize(2*(m+n));
    compute_tz(d,g,c,y,t,z);
}


Void KrawczykOptimiser::compute_tz(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                   const ExactIntervalVectorType& y, ExactIntervalType& t, ExactIntervalVectorType& z) const
{
    ARIADNE_ASSERT(d.size()>0u);
    //static const double EPS=1.0/8;
    static const float min_float=std::numeric_limits<float>::min();

    const SizeType m=g.argument_size();
    const SizeType n=g.result_size();

    // Compute the image of y under the constraint function
    ExactIntervalVectorType gy=g(y);
    gy+=ExactIntervalVectorType(gy.size(),ExactIntervalType(-min_float,+min_float));
    ExactIntervalVectorType my=midpoint(y);
    ExactIntervalVectorType mgy=g(my);

    // Find the range of possible values of the optimal t
    // This range is too pessimistic
    t=ExactIntervalType(+inf,+inf);
    for(SizeType j=0; j!=n; ++j) {
        t=min(t,c[j]-gy[j]);
        t=min(t,gy[j]-c[j]);
    }
    for(SizeType i=0; i!=m; ++i) {
        t=min(t,d[i]-y[i]);
        t=min(t,y[i]-d[i]);
    }

    // Find the range of possible values of the optimal t
    FloatDP tmin=+inf;
    FloatDP tmax=+inf;
    for(SizeType j=0; j!=n; ++j) {
        tmax=min(tmax,sub(up,c[j].upper_bound(),gy[j].lower_bound()));
        tmax=min(tmax,sub(up,gy[j].upper_bound(),c[j].lower_bound()));
        tmin=min(tmin,sub(down,c[j].upper_bound(),mgy[j].upper_bound()));
        tmin=min(tmin,sub(down,mgy[j].lower_bound(),c[j].lower_bound()));
    }
    for(SizeType i=0; i!=m; ++i) {
        tmin=min(tmin,sub(up,d[i].upper_bound(),y[i].lower_bound()));
        tmax=min(tmax,sub(up,y[i].upper_bound(),d[i].lower_bound()));
        tmin=min(tmin,sub(down,d[i].upper_bound(),my[i].upper_bound()));
        tmax=min(tmax,sub(down,my[i].lower_bound(),d[i].lower_bound()));
    }
    tmin-=0.0625;
    t=ExactIntervalType(tmin,tmax);


    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(SizeType j=0; j!=n; ++j) {
        z[j]=max(c[j].upper_bound()-gy[j]-t,0.0);
        z[n+j]=max(gy[j]-c[j].lower_bound()-t,0.0);
    }
    for(SizeType i=0; i!=m; ++i) {
        z[2*n+i]=max(d[i].upper_bound()-y[i]-t,0.0);
        z[2*n+m+i]=max(y[i]-d[i].lower_bound()-t,0.0);
    }

    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(SizeType j=0; j!=n; ++j) {
        z[j]=ExactIntervalType(0.0,c[j].upper_bound()-mgy[j].lower_bound()-tmin);
        z[n+j]=ExactIntervalType(0.0,mgy[j].upper_bound()-c[j].lower_bound()-tmin);
    }
    for(SizeType i=0; i!=m; ++i) {
        z[2*n+i]=ExactIntervalType(0.0,d[i].upper_bound()-my[i].lower_bound()-tmin);
        z[2*n+m+i]=ExactIntervalType(0.0,my[i].upper_bound()-d[i].lower_bound()-tmin);
    }

    CONCLOG_PRINTLN_AT(9,"  d="<<d<<", c="<<c<<", y="<<y<<", g(y)="<<gy<<", t="<<t<<", z="<<z);

}


Void KrawczykOptimiser::
minimisation_step(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g,
                  ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z) const
{
    const SizeType m=f.argument_size();
    const SizeType n=g.result_size();

    Differential<UpperIntervalType> ddf=f.evaluate(Differential<UpperIntervalType>::variables(2,y));
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    ExactIntervalMatrixType H(m,m);
    set_hessian(H,ddf);
    for(SizeType j=0; j!=n; ++j) { add_hessian(H,-x[j],ddg[j]); }

    ExactIntervalMatrixType A(m,n);
    set_jacobian_transpose(A,ddg);

    CONCLOG_PRINTLN_AT(9,"f="<<f<<",g="<<g);
    CONCLOG_PRINTLN_AT(9,"x="<<x<<" y="<<y<<" z="<<z);
    CONCLOG_PRINTLN_AT(9,"A="<<A);
    CONCLOG_PRINTLN_AT(9,"H="<<H);

    ARIADNE_NOT_IMPLEMENTED;

}



Void KrawczykOptimiser::feasibility_step(const ValidatedVectorMultivariateFunction& g,
                                         ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    ARIADNE_NOT_IMPLEMENTED;
    const SizeType m=y.size();
    const SizeType n=x.size();

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    CONCLOG_PRINTLN_AT(9,"A="<<A);

    // H is the Hessian matrix Hik = xj*dgj/dyidyk
    ExactIntervalMatrixType H(m,m);
    for(SizeType j=0; j!=n; ++j) {
        add_hessian(H,x[j],ddg[j]);
    }
    CONCLOG_PRINTLN_AT(9," H="<<H);

    FloatDPMatrix mA=midpoint(A);
    CONCLOG_PRINTLN_AT(9," mA="<<mA);
    FloatDPMatrix mH=midpoint(H);
    CONCLOG_PRINTLN_AT(9," mH="<<mH);

    RawFloatDPVector mD(n);
    for(SizeType j=0; j!=n; ++j) { mD[j]=midpoint(x[j])/midpoint(z[j]); }
    CONCLOG_PRINTLN_AT(9," mD="<<mD);

    FloatDPMatrix& mS=mH;
    adat(mS,mA,mD);
    CONCLOG_PRINTLN_AT(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    CONCLOG_PRINTLN_AT(9,"mSinv="<<mSinv);
}

// Feasibility step for dual (inequality constrained) problem without using slack variables
// FIXME: Do we need a slackness parameter mu? Probably not; hopefully the infinities are kept in check...
// This method has the advantage of not needing to update the primal variables
Void KrawczykOptimiser::feasibility_step(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                         ExactIntervalVectorType& y, ExactIntervalType& t) const
{
    const SizeType m=d.size();
    const SizeType n=c.size();

    // Compute function values
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // gy is the vector of values of g(y)
    ExactIntervalVectorType gy(n);
    for(SizeType j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }

    // z is the vector of slack variables z[k]=cu[k]-gy[k]-t or z[k]=gy[k]-cl[k]-t
    ExactIntervalVectorType z(2*(m+n));
    for(SizeType j=0; j!=n; ++j) { z[j]=d[j].upper_bound()-gy[j]-t; z[n+j]=gy[j]-d[j].lower_bound()-t; }
    for(SizeType i=0; i!=m; ++i) { z[i]=c[2*n+i].upper_bound()-y[i]-t; z[2*n+m+i]=y[i]-c[i].lower_bound()-t; }

    ExactIntervalVectorType zr(2*(m+n));
    for(SizeType k=0; k!=2*(m+n); ++k) { zr[k]=1.0/z[k]; }

    ExactIntervalVectorType D(2*(m+n));
    for(SizeType k=0; k!=2*(m+n); ++k) { D[k]=zr[k]*zr[k]; }

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { A[i][j]=ddg[j][i]; } }

    // A is the sum of scaled Hessian matrices hi1i2=zj*ddgj/dyi1yi2
    ExactIntervalMatrixType H(m,m);

    ExactIntervalMatrixType SE(m+1,m+1);
    // SE[0:m][0:m] is the matrix H+/-A(D1+D2)AT+(D3+D4) where D=z.z
    for(SizeType i1=0; i1!=m; ++i1) { for(SizeType i2=0; i2!=m; ++i2) { SE[i1][i2]=H[i1][i2];
        for(SizeType j=0; j!=n; ++j) { SE[i1][i2]+=A[i1][j]*(D[j]+D[n+j])*A[i2][j]; }
    } }
    for(SizeType i=0; i!=m; ++i) { SE[i][i]+=(D[2*n+i]+D[2*n+m+i]); }
    // SE[m][0:m]=SE[0:m][m] is the vector A(D1-D2)e+(D3-D4)e
    for(SizeType i=0; i!=m; ++i) { SE[i][m]=D[2*n+i]-D[2*n+m+i];
        for(SizeType j=0; j!=n; ++j) { SE[i][m]+=A[i][j]*(D[j]-D[n+j]); }
        SE[m][i]=SE[i][m];
    }
    // SE[m][m] is the scalar eT(D1+D2)e+eT(D3+D4)e
    for(SizeType k=0; k!=2*(m+n); ++k) { SE[m][m]+=D[k]; }

    // Vector of residuals
    ExactIntervalVectorType re(m+1);
    for(SizeType i=0; i!=m; ++i) { re[i]+=(zr[2*n+i]-zr[2*n+m+i]);
        for(SizeType j=0; j!=n; ++j) { re[i]+=A[i][j]*(zr[j]-zr[n+j]); }
    }
    for(SizeType j=0; j!=n; ++j) { re[m]+=(zr[j]+zr[n+n]); }
    for(SizeType i=0; i!=m; ++i) { re[m]+=(zr[2*n+i]+zr[2*n+m+i]); }

    // Compute inverse Jacobian matrix
    ExactIntervalMatrixType JE;
    try {
        JE=inverse(midpoint(SE));
    }
    catch(const SingularMatrixException& e) {
        ARIADNE_WARN("Matrix S="<<midpoint(SE)<<" is not invertible");
        CONCLOG_PRINTLN_AT(1,"WARNING: Matrix S="<<midpoint(SE)<<" is not invertible");
        throw e;
    }

    // Krawczyk step
    ExactIntervalVectorType dyt=prod(JE,ExactIntervalVectorType(midpoint(re)))+prod(ExactIntervalMatrixType::identity(m+1)-prod(JE,SE),re-midpoint(re));

    // Extract y and t
    ExactIntervalVectorType yt=join(y,t);
    ExactIntervalVectorType nyt=yt+dyt;

    yt=intersection(yt,nyt);
    y=project(yt,range(0,m));
    t=yt[m];
}




Void KrawczykOptimiser::feasibility_step(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                         ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    const SizeType m=d.size();
    const SizeType n=c.size();
    const SizeType o=2*(m+n);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==o);
    ARIADNE_ASSERT(y.size()==m);
    ARIADNE_ASSERT(z.size()==o);

    ExactIntervalVectorType yt=join(y,t);
    CONCLOG_PRINTLN_AT(9,"m="<<m<<" n="<<n);
    CONCLOG_PRINTLN_AT(9,"x="<<x<<" yt="<<yt<<" z="<<z);

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));
    CONCLOG_PRINTLN_AT(9,"  ddg="<<ddg);

    // gy is the vector of values of g(y)
    UpperIntervalVectorType gy(n); for(SizeType j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }
    CONCLOG_PRINTLN_AT(9,"  g(y)="<<gy<<" ");

    // A is the transpose derivative matrix aij=dgj/dyi, extended with a column of ones
    UpperIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    CONCLOG_PRINTLN_AT(9," A="<<A<<" ");

    // H is the Hessian matrix Hik = (xcuj-xclj)*dgj/dyidyk
    UpperIntervalMatrixType H(m,m);
    for(SizeType j=0; j!=n; ++j) {
        add_hessian(H,x[j]-x[n+j],ddg[j]);
    }
    CONCLOG_PRINTLN_AT(9," H="<<H);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    UpperIntervalVectorType gye(o);
    for(SizeType j=0; j!=n; ++j) { gye[j]=gy[j]-c[j].upper_bound()+t; gye[n+j]=c[j].lower_bound()-gy[j]+t; }
    for(SizeType i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper_bound()+t; gye[2*n+m+i]=d[i].lower_bound()-y[i]+t; }
    CONCLOG_PRINTLN_AT(9,"  GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    UpperIntervalMatrixType AE(m+1,o);
    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    for(SizeType i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    for(SizeType k=0; k!=o; ++k) { AE[m][k]=1; }
    UpperIntervalMatrixType AET=transpose(AE);

    FloatDPMatrix mA=midpoint(A);
    FloatDPMatrix mAE=midpoint(AE);
    FloatDPMatrix mAET=midpoint(AET);
    FloatDPMatrix mH=midpoint(H);
    RawFloatDPVector mx=midpoint(x);
    RawFloatDPVector myt=midpoint(yt);
    RawFloatDPVector mz=midpoint(z);
    RawFloatDPVector mDE=ediv(mx,mz);


    // Construct the symmetric matrix and its inverse
    //FloatDPMatrix S(m+1,m+1); adat(S,AE,DE);
    //CONCLOG_PRINTLN_AT(9,"S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //CONCLOG_PRINTLN_AT(9,"S="<<S);
    FloatDPMatrix mS=feasibility_adat(mH,mA,mDE);
    CONCLOG_PRINTLN_AT(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    CONCLOG_PRINTLN_AT(9,"mSinv="<<mSinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    UpperIntervalVectorType rx=emul(mx,mz);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType ryt=-feasibility_mul(mA,mx); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType rz=midpoint(gye)+mz;
    CONCLOG_PRINTLN_AT(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    // Construct the errors on the residuals ([M]-M)([x]-x)
    UpperIntervalVectorType ex=x-mx;
    UpperIntervalVectorType eyt=yt-myt;
    UpperIntervalVectorType ez=z-mz;
    UpperIntervalMatrixType eA=A-mA;
    UpperIntervalMatrixType eH=H-mH;

    UpperIntervalVectorType erx=2.0*emul(ex,ez);
    UpperIntervalVectorType eryt=UpperIntervalMatrixType(AE-mAE)*ex;
    UpperIntervalVectorType erz=UpperIntervalMatrixType(AET-mAET)*eyt;
    CONCLOG_PRINTLN_AT(9,"erx="<<erx<<" eryt="<<eryt<<" erz="<<erz);

    rx+=2.0*emul(ex,ez);
    ryt+=UpperIntervalMatrixType(AE-mAE)*ex;
    rz+=UpperIntervalMatrixType(AET-mAET)*eyt;
    CONCLOG_PRINTLN_AT(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;

    // Compute the error differences
    UpperIntervalVectorType erxdz=ediv(erx,mz);
    UpperIntervalVectorType edyt=(mSinv*mAE)*erxdz + mSinv*eyt - (mSinv*(mAE*DiagonalMatrix<FloatDP>(mDE))) * ez;
    UpperIntervalVectorType edz=-erz-feasibility_trmul(mA,edyt);
    UpperIntervalVectorType edx=-ediv(UpperIntervalVectorType(erx+emul(mx,edz)),mz);
    CONCLOG_PRINTLN_AT(9,"edx="<<edx<<" edyt="<<edyt<<" edz="<<edz);

    // Compute the error differences
    UpperIntervalVectorType eerr=prod(mAE,ediv(esub(erx,emul(mx,erz)),mz))-eryt;
    CONCLOG_PRINTLN_AT(9,"  eerr="<<eerr);
    UpperIntervalVectorType eedyt=prod(mSinv,eerr);
    UpperIntervalVectorType eedz=-erz-feasibility_trmul(mA,eedyt);
    UpperIntervalVectorType eedx=-ediv(UpperIntervalVectorType(erx+emul(mx,eedz)),mz);
    CONCLOG_PRINTLN_AT(9,"eedx="<<eedx<<" eedyt="<<eedyt<<" eedz="<<eedz);


    // Compute the differences
    UpperIntervalVectorType rr=prod(mAE,ediv(esub(rx,emul(mx,rz)),mz))-ryt;
    UpperIntervalVectorType dyt=prod(mSinv,rr);
    UpperIntervalVectorType dz=-rz-feasibility_trmul(mA,dyt);
    UpperIntervalVectorType dx=-ediv(UpperIntervalVectorType(rx+emul(mx,dz)),mz);
    CONCLOG_PRINTLN_AT(9,"dx="<<dx<<" dyt="<<dyt<<" dz="<<dz);

    UpperIntervalVectorType nx,ny,nyt,nz; FloatDP nt;
    nx=mx+dx;
    nyt=myt+dyt;
    nz=mz+dz;

    x=intersection(x,nx);
    yt=intersection(yt,nyt);
    z=intersection(z,nz);

    y=project(yt,range(0,m));
    t=yt[m];
}

*/


} // namespace Ariadne
