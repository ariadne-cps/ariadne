/***************************************************************************
 *            algebra/eigenvalues.cpp
 *
 *  Copyright 2015-16  Pieter Collins
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

#include "utility/macros.hpp"
#include "utility/array.hpp"
    #include "numeric/logical.hpp"
    #include "numeric/dyadic.hpp"
    #include "numeric/floats.hpp"
#include "numeric/float_bounds.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/eigenvalues.hpp"

namespace Ariadne {

template<class X> UpperTriangularMatrix<X>::UpperTriangularMatrix(UpperHessenbergMatrix<X> const& A)
    : UpperTriangularMatrix<X>(A.row_size(),A.zero_element().precision())
{
    SizeType n=A.column_size();
    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=i; j!=n; ++j) {
            this->set(i,j,A.get(i,j));
        }
    }
}

template<class X> auto UpperTriangularMatrix<X>::_mul(Matrix<X> const& A, UpperTriangularMatrix<X> const& U) -> Matrix<X> {
    SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element();
    Matrix<X> R(m,n,z);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            for(SizeType k=0; k<=j; ++k) {
                R[i][j]+=A[i][k]*U.get(k,j);
            }
        }
    }
    return R;
}

template<class X> auto UpperTriangularMatrix<X>::_mul(UpperTriangularMatrix<X> const& U, Matrix<X> const& A) -> Matrix<X> {
    SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element();
    Matrix<X> R(m,n,z);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            for(SizeType k=i; k!=n; ++k) {
                R[i][j]+=U.get(i,k)*A[k][j];
            }
        }
    }
    return R;
}

template<class X> auto UpperTriangularMatrix<X>::_write(OutputStream& os, SizeType wdth) const -> OutputStream& {
    UpperTriangularMatrix<X>const& A=*this; SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element();
    for(SizeType i=0; i!=m; ++i) {
        os << "[";
        for(SizeType j=0; j!=i; ++j) {
            os << std::setw(wdth) << 0 << ",";
        }
        for(SizeType j=i; j!=n; ++j) {
            os << std::setw(wdth) << A.get(i,j) << (j+1u==n?"]\n":",");
        }
    }
    return os;
}



template<class X> UpperHessenbergMatrix<X>::UpperHessenbergMatrix(UpperTriangularMatrix<X> const& A)
    : UpperHessenbergMatrix<X>(A.row_size(),A.zero_element())
{
    SizeType m=A.row_size(); SizeType n=A.column_size();
    for(SizeType i=0; i!=m; ++i) {
        const SizeType j0=std::max(i,(SizeType)1u)-1u;
        for(SizeType j=j0; j!=n; ++j) {
            this->set(i,j,A.get(i,j));
        }
    }
}

template<class X> UpperHessenbergMatrix<X>::UpperHessenbergMatrix(Matrix<X> const& A)
    : UpperHessenbergMatrix<X>(A.row_size(),A.zero_element())
{
    SizeType m=A.row_size(); SizeType n=A.column_size();
    for(SizeType i=0; i!=m; ++i) {
        const SizeType j0=std::max(i,(SizeType)1u)-1u;
        for(SizeType j=j0; j!=n; ++j) {
            this->set(i,j,A.get(i,j));
        }
    }
}

template<class X> auto UpperHessenbergMatrix<X>::_mul(Matrix<X> const& A, UpperHessenbergMatrix<X> const& UH) -> Matrix<X> {
    ARIADNE_PRECONDITION(A.column_size()==UH.row_size());
    ARIADNE_PRECONDITION(UH.row_size()==UH.column_size());
    SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element();
    Matrix<X> R(m,n,z);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            for(SizeType k=0; k!=std::min(j+2u,n); ++k) {
                R[i][j]+=A[i][k]*UH.get(k,j);
            }
        }
    }
    return R;
}

template<class X> auto UpperHessenbergMatrix<X>::_write(OutputStream& os, SizeType wdth) const -> OutputStream& {
    UpperHessenbergMatrix<X>const& A=*this; SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element();
    for(SizeType i=0; i!=m; ++i) {
        os << "[";
        const SizeType j0=std::max(i,(SizeType)1u)-1u;
        for(SizeType j=0; j!=j0; ++j) {
            os << std::setw(wdth) << 0 << ",";
        }
        for(SizeType j=j0; j!=n; ++j) {
            os << std::setw(wdth) << A.get(i,j) << (j+1u==n?"]\n":",");
        }
    }
    return os;
}


//Q0=H0; A1 = H0*A0*H0 = Q0*A0*Q0'
//Q1=Q0*H1=H0*H1; A2=H1*A1*H1=H1*H0*A0*H0*H1=Q1T*A0*Q1
//Q2=Q1*H2=H0*H1*H2; A3=H2*A2*H2=H2*H1*H0*A0*H0*H1*H2=Q2T*A0*Q2
//An=QnT * A * Qn
//A = Qn * An * QnT
template<class X> Pair<HouseholderProductMatrix<X>,UpperHessenbergMatrix<X>> upper_hessenberg_factorisation(Matrix<X> A) {
    // Set alpha=nrm2(A(1:n,0));
    Matrix<X> old_A=A;
    SizeType n=A.row_size(); X z=A.zero_element();
    HouseholderProductMatrix<X> Q(n);
    for(SizeType j=0; j!=n-2; ++j) {
        //std::cerr<<"  B="<<A<<"\n";
        X alpha=z; for(SizeType i=j+1; i!=n; ++i) { alpha+=sqr(A[i][j]); } alpha=sqrt(alpha);
        Vector<X> v=column(A,j); for(SizeType i=0; i<=j; ++i) { v[i]=z; } v[j+1]-=alpha;
        //std::cerr<<"    v="<<v<<"\n";
        HouseholderMatrix<X> H(v);
        Q*=H;
        Matrix<X> new_A=H*A*H;
        assert(refines(A,H*new_A*H));
        assert(refines(old_A,Q*new_A*inverse(Q)));
        A=new_A;
    }
    //std::cerr<<"  B="<<A<<"\n";
    UpperHessenbergMatrix<X> H(A);
    assert(refines(to_matrix(H),A));
    return std::make_pair(Q,H);
}

template<class X> Pair<X,Vector<X>> eigenvalue_solve(Matrix<X> const& A, X mu, Vector<X> v, ApproximateTag) {
    // Solve equations Av-muv=0; v'v=1
    SizeType n=v.size(); X z=v.zero_element();

    Vector<X> r(n+1,z);
    Matrix<X> Dr(n+1,n+1,z);
    project(Dr,range(0,n),range(0,n))=A;

    //std::cerr<<"v="<<v<<", mu="<<mu<<"\n";
    SizeType step=0;
    while(step<5) {
        ++step;

        for(SizeType i=0; i!=n; ++i) { Dr[i][n]=-v[i]; Dr[n][i]=-v[i]; Dr[i][i]=A[i][i]-mu; }
        r=join(A*v-mu*v,(1-dot(v,v))/2);

        Vector<X> dx=lu_solve(Dr,r);
        assert(not is_nan(dx[n].raw()));
        Vector<X> new_v = v - Vector<X>(project(dx,range(0,n)));
        X new_mu = mu - dx[n];
        //std::cerr<<"  new_v="<<new_v<<", new_mu="<<new_mu<<"\n";
        v=new_v; mu=new_mu;
        //std::cerr<<"v="<<v<<", mu="<<mu<<"\n";
    }
    return std::make_pair(mu,v);
}

template<class X> Pair<X,Vector<X>> eigenvalue_solve(Matrix<X> const& A, X mu, Vector<X> v, ValidatedTag) {
    // Solve equations Av-muv=0; v'v=1
    SizeType n=v.size(); X z=v.zero_element();

    bool is_refinement=false;
    bool is_inconsistent=false;

    SizeType step=0;
    while(step<3) {
        ++step;

        Vector<X> mid_v=midpoint(v); X mid_mu=midpoint(mu);
        Vector<X> mid_r=join(A*mid_v-mid_mu*mid_v,(1-dot(mid_v,mid_v))/2);
        Matrix<X> Dr(n+1,n+1,z);
        project(Dr,range(0,n),range(0,n))=A-mu*Matrix<X>::identity(n,z);
        for(SizeType i=0; i!=n; ++i) { Dr[i][n]=-v[i]; Dr[n][i]=-v[i]; }

        Vector<X> dx=lu_solve(Dr,mid_r);
        Vector<X> new_v = mid_v - Vector<X>(project(dx,range(0,n)));
        X new_mu = mid_mu - dx[n];

        if(not is_refinement) {
           is_refinement = refines(new_mu,mu) and refines(new_v,v);
            if(is_refinement) { std::cerr << "Refinement!\n"; }
        }

        //std::cerr<<"  new_v="<<new_v<<", new_mu="<<new_mu<<"\n";
        is_inconsistent = inconsistent(new_mu,mu) or inconsistent(new_v,v);
        assert(not is_inconsistent);

        //std::cerr<<"v="<<v<<", mu="<<mu<<"\n";
        v=new_v; mu=new_mu;

        //v = refinement(v,new_v);
        //mu = refinement(mu,new_mu);
    }
    return std::make_pair(mu,v);

}

template<class X> Pair<X,Vector<X>> eigenvalue_solve(Matrix<X> const& A, X mu, Vector<X> v) {
    return eigenvalue_solve(A,mu,v, Paradigm<X>());
}


template<class X> Pair<OrthogonalMatrix<X>,UpperTriangularMatrix<X>> gram_schmidt(Matrix<X> A) {
    SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element(); auto pr=z.precision();
    ARIADNE_ASSERT(m==n);
    OrthogonalMatrix<X> Q(n,z.precision()); UpperTriangularMatrix<X> R(n,pr);
    X r=z; X d=z;
    for(SizeType j=0; j!=n; ++j) {
        Vector<X> v=A.column(j);
        for(SizeType i=0; i!=j; ++i) {
            r = dot(v,static_cast<Vector<X>>(Q.column(i)));
            R.set(i,j,r);
            v-=r*static_cast<Vector<X>>(Q.column(i));
        }
        d=two_norm(v);
        R.set(j,j,d);
        v/=d;
        Q.column(j)=v;
    }
    return std::make_pair(std::move(Q),std::move(R));
}

template<class X> Matrix<X> GivensProductMatrix<X>::rmul(Matrix<X> A) const {
    SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element(); auto pr=z.precision();
    for(SizeType j=0; j+1u<n; ++j) {
        auto& alpha=this->_rotations[j].alpha;
        auto& beta=this->_rotations[j].beta;
        for(SizeType i=0; i<m; ++i) {
            X& p=A[i][j]; X& q=A[i][j+1u];
            X t=p*alpha+q*beta; q=q*alpha-p*beta; p=t;
        }
    }
    return A;
}

template<class X> UpperHessenbergMatrix<X> GivensProductMatrix<X>::rmul(UpperTriangularMatrix<X> const& U) const {
    UpperHessenbergMatrix<X> A(U);
    SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element(); auto pr=z.precision();
    for(SizeType j=0; j+1u<n; ++j) {
        auto& alpha=this->_rotations[j].alpha;
        auto& beta=this->_rotations[j].beta;
        for(SizeType i=0; i<std::max(m,j+2u); ++i) {
            X& p=A.at(i,j); X& q=A.at(i,j+1);
            X np=p*alpha+q*beta; X nq=q*alpha-p*beta; p=np; q=nq;
        }
    }
    return A;
}

template<class X> Pair<GivensProductMatrix<X>,UpperTriangularMatrix<X>> gram_schmidt(UpperHessenbergMatrix<X> A) {
    SizeType m=A.row_size(); SizeType n=A.column_size(); X z=A.zero_element(); auto pr=z.precision();
    ARIADNE_ASSERT(m==n);
    GivensProductMatrix<X> Q(n,z.precision());
    for(SizeType i=0; i+1u<n; ++i) {
        X const& a=A.get(i,i);
        X const& b=A.get(i+1,i);
        X rho=sqrt(sqr(a)+sqr(b));
        X alpha = a/rho;
        X beta = b/rho;
        // Would prefer use a/sqrt(a^2+b^2) = 1/(sqrt(1+b^2/a^2)) for accuracy, but problems with sign
        // X as=sqr(a); X bs=sqr(b); X alpha=rec(sqrt(1+bs/as)); X beta=rec(sqrt(1+as/bs));
        Q.set(i,alpha,beta);
        A.set(i,i,rho);
        A.set(i+1,i,0);
        for(SizeType j=i+1; j!=n; ++j) {
            X& p=A.at(i,j); X& q=A.at(i+1,j);
            X np=p*alpha+q*beta; X nq=q*alpha-p*beta; p=np; q=nq;
        }
    }
    UpperTriangularMatrix<X> R(A);
    return std::make_pair(std::move(Q),std::move(R));
}

template<class X> Void qr_step(UpperHessenbergMatrix<X>& H) {
    Pair<GivensProductMatrix<X>,UpperTriangularMatrix<X>> GT=gram_schmidt(H);
    GivensProductMatrix<X> const& G=GT.first;
    UpperTriangularMatrix<X> const& T=GT.second;
//    Matrix<X> A=H; Matrix<X> Q=G; Matrix<X> R=T; Matrix<X> B=R*Q; Matrix<X> QBQT=(Q*B*transpose(Q)); ARIADNE_ASSERT(refines(A,QBQT));
    H=T*G;
}


/*
[ 1, 0, 0, 0, 0] [t00,t01, 0 , 0 , 0 ]   [t00,t01         , 0          , 0   , 0 ]
[ 0, a, b, 0, 0] [ 0 ,t11,t12, 0 , 0 ] = [ 0 , a*t11+b*t12, a*t12+b*t22,b*t23, 0 ]
[ 0,-b, a, 0, 0] [ 0 ,t12,t22,t23, 0 ]   [ 0 ,-b*t11+a*t12,-b*t12+a*t22,a*t23, 0 ]
[ 0, 0, 0, 1, 0] [ 0 , 0 ,t23,t33,t34]   [ 0 , 0          , t23        ,t33  ,t34]
[ 0, 0, 0, 0, 1] [ 0 , 0 , 0 ,t34,t44]   [ 0 , 0          , 0          ,t34  ,t44  ]
*/
#warning
// Algorithm in progress
template<class X> Void qr_step(SymmetricTridiagonalMatrix<X>& T) {
    SizeType n=T.size(); X z=T.zero_element();
    X lt=z;
    GivensProductMatrix<X> G(n,z);
    Matrix<X> A=T;
    for (SizeType i=0; i!=n-1; ++i) {
        X al=A.get(i,i); X be=A.get(i+1,i); X ga=sqrt(sqr(al)+sqrt(be)); al/=ga; be/=ga;
        G.set(i,al,be);
        X ti0i0= al*T.get(i,i)+be*T.get(i,i+1);
        // X ti1i0=-be*T.get(i,i)+al*T.get(i,i+1); // al,be defined so this is 0.
        X ti0i1= al*T.get(i,i+1)+be*T.get(i+1,i+1);
        X ti1i1=-be*T.get(i,i+1)+al*T.get(i+1,i+1);
        A.set(i,i,ti0i0);
        // Todo
    }
    assert(false);
}

template Pair<FloatMPApproximation,Vector<FloatMPApproximation>> eigenvalue_solve(Matrix<FloatMPApproximation> const& A, FloatMPApproximation mu, Vector<FloatMPApproximation> v);

template class UpperTriangularMatrix<FloatMPBounds>;
template class UpperHessenbergMatrix<FloatMPBounds>;
template class GivensProductMatrix<FloatMPBounds>;
template Pair<HouseholderProductMatrix<FloatMPBounds>,UpperHessenbergMatrix<FloatMPBounds>> upper_hessenberg_factorisation(Matrix<FloatMPBounds>);
template Pair<OrthogonalMatrix<FloatMPBounds>,UpperTriangularMatrix<FloatMPBounds>> gram_schmidt(Matrix<FloatMPBounds>);
template Pair<GivensProductMatrix<FloatMPBounds>,UpperTriangularMatrix<FloatMPBounds>> gram_schmidt(UpperHessenbergMatrix<FloatMPBounds>);
template Void qr_step(UpperHessenbergMatrix<FloatMPBounds>& H);

//template<class X> Pair<X,Vector<X>> eigenvalue_solve(Matrix<X> const& A, X mu, Vector<X> v) {
template Pair<FloatMPBounds,Vector<FloatMPBounds>> eigenvalue_solve(Matrix<FloatMPBounds> const& A, FloatMPBounds mu, Vector<FloatMPBounds> v);
//template Pair<FloatMPBounds,Vector<FloatMPBounds>> eigenvalue_solve(Matrix<FloatMPBounds> const& A, FloatMPBounds, Vector<FloatMPBounds>);

} // namespace Ariadne
