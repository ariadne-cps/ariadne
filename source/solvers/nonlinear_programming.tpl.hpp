/***************************************************************************
 *            solvers/nonlinear_programming.tpl.hpp
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

#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/diagonal_matrix.hpp"
#include "algebra/symmetric_matrix.hpp"

#include "solvers/nonlinear_programming.hpp"

namespace Ariadne {

namespace {

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

} // namespace

template<class X1, class IVL2> inline decltype(auto) dot(Vector<X1> const& x1, Box<IVL2> const& bx2) {
    return dot(x1,cast_vector(bx2)); }

template<> inline UpperIntervalType dot<UpperIntervalType,ExactIntervalType>(Vector<UpperIntervalType> const& bx1, Vector<ExactIntervalType> const& bx2) {
    return dot(bx1,Vector<UpperIntervalType>(bx2));
}

template<class X> Matrix<X> join(Matrix<X> const& A1, Matrix<X> const& A2, Matrix<X> const& A3) {
    SizeType m=A1.row_size(); SizeType n1=A1.column_size(); SizeType n2=A2.column_size(); SizeType n3=A3.column_size();
    Matrix<X> A123(m,n1+n2+n3,(A1.zero_element()+A2.zero_element()+A3.zero_element()));
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
Void iqadat(Matrix<X>& S, const Matrix<X>& A, const DiagonalMatrix<X>& D)
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
Void iqatda(Matrix<X>& S, const Matrix<X>& A, const DiagonalMatrix<X>& D)
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

template<class X>
Matrix<X> qatda(Matrix<X> Q, const Matrix<X>& A, const DiagonalMatrix<X>& D)
{
    iqatda(Q,A,D);
    return Q;
}

// Compute S=ADA^T, where D is diagonal.
template<class X>
Matrix<X> adat(const Matrix<X>& A, const DiagonalMatrix<X>& D)
{
    const SizeType m=A.row_size();
    Matrix<X> S=Matrix<X>::zero(m,m);
    iqadat(S,A,D);
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
template<class X> SymmetricMatrix<X> feasibility_adat(const SymmetricMatrix<X>& H, const Matrix<X>& A, const DiagonalMatrix<X>& D)
{
    ARIADNE_NOT_IMPLEMENTED;

    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    ARIADNE_ASSERT(H.row_size()==m);
    ARIADNE_ASSERT(H.column_size()==m);
    ARIADNE_ASSERT(D.size()==2*(m+n));
    SymmetricMatrix<X> S(m+1,H.zero_element());

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

template<class T>
struct PrimalDualComplementaryMatrix {
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

    PrimalDualComplementaryData<T> solve(PrimalDualComplementaryData<T> const& d) const {
        Vector<T> dr = lu_solve(this->assemble(),d.assemble());
        return PrimalDualComplementaryData(this->A.row_size(),this->A.column_size(),dr);
    }
};

} // namespace

// [ -I   A   0   0  ]
// [  0   Q   A'  I  ]
// [  Y   0   W   0  ]
// [  0   Z   0   X  ]

template<class T> auto
slack_primal_dual_complementary_matrix(SymmetricMatrix<T>const& Q, Matrix<T>const& A,
                                       DiagonalMatrix<T>const& W, DiagonalMatrix<T>const& X, DiagonalMatrix<T>const& Y, DiagonalMatrix<T>const& Z) -> Matrix<T>
{
    SizeType m=A.row_size(); SizeType n=A.column_size();
    IdentityMatrix Im(m), In(n);
    T z=Q.zero_element();
    Matrix<T> S(2*(m+n),2*(m+n),z);
    assign_diagonal(S[range(0,m)][range(0,m)],-Im);
    assign_dense(S[range(0,m)][range(m,m+n)],A);
    assign_dense(S[range(m,m+n)][range(m,m+n)],Q);
    assign_dense(S[range(m,m+n)][range(m+n,2*m+n)],transpose(A));
    assign_diagonal(S[range(m,m+n)][range(2*m+n,2*(m+n))],In);
    assign_diagonal(S[range(m+n,2*m+n)][range(0,m)],Y);
    assign_diagonal(S[range(m+n,2*m+n)][range(m+n,2*m+n)],W);
    assign_diagonal(S[range(2*m+n,2*(m+n))][range(m,m+n)],Z);
    assign_diagonal(S[range(2*m+n,2*(m+n))][range(2*m+n,2*(m+n))],X);
    return S;
}



//! \brief A matrix defining the linear equations for a full Newton step for slack, primal, dual and complementary variables.
//! \details Form is \f[ \begin{pmatrix}Y&0&W&0\\0&Q&A^T&I\\-I&A&0&0\\0&Z&0&X\end{pmatrix} . \f]
template<class T> struct SlackPrimalDualComplementaryMatrix {
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
        return SlackPrimalDualComplementaryData<T>::disassemble(this->A.row_size(),this->A.column_size(),dr);
    }
};

template<class T> auto PrimalDualComplementaryData<T>::disassemble(SizeType m, SizeType n, Vector<T> r) -> PrimalDualComplementaryData<T> {
    ARIADNE_DEBUG_PRECONDITION(r.size()==m+2*n);
    auto x_=r[range(0,n)]; auto y_=r[range(n,m+n)]; auto z_=r[range(m+n,m+2*n)];
    return PrimalDualComplementaryData<T>(x_,y_,z_);
}

template<class T> auto PrimalDualComplementaryData<T>::assemble() const -> Vector<T> {
return join(join(this->x,this->y),this->z); }

template<class T> auto SlackPrimalDualComplementaryData<T>::disassemble(SizeType m, SizeType n, Vector<T> r) -> SlackPrimalDualComplementaryData<T> {
    ARIADNE_DEBUG_PRECONDITION(r.size()==2*(m+n));
    auto w_=r[range(0,m)]; auto x_=r[range(m,m+n)]; auto y_=r[range(m+n,2*m+n)]; auto z_=r[range(2*m+n,2*(m+n))];
    return SlackPrimalDualComplementaryData<T>(w_,x_,y_,z_);
}

template<class T> auto SlackPrimalDualComplementaryData<T>::assemble() const -> Vector<T> {
return join(join(this->w,this->x),join(this->y,this->z)); }


} // namespace Ariadne
