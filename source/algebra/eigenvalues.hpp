/***************************************************************************
 *            algebra/eigenvalues.hpp
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

/*! \file algebra/eigenvalues.hpp
 *  \brief
 */



#ifndef ARIADNE_EIGENVALUES_HPP
#define ARIADNE_EIGENVALUES_HPP

#include "utility/array.hpp"
#include "vector.hpp"
#include "covector.hpp"
#include "matrix.hpp"

namespace Ariadne {

template<class T> List<T> reverse(std::vector<T> const& lst) {
    return List<T>(std::vector<T>(lst.rbegin(),lst.rend()));
};


template<class M> decltype(auto) to_matrix(M const& A) {
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    auto z=A.zero_element();
    Matrix<decltype(z)> R(m,n,z);
    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { R[i][j]=A.get(i,j); } }
    return R;
}


template<class X> struct RankOneMatrix {
    Vector<X> _v; Covector<X> _w;
    friend Matrix<X>& operator+=(Matrix<X>& A, RankOneMatrix<X> const& B) {
        for(SizeType i=0; i!=B._v.size(); ++i) { for(SizeType j=0; j!=B._w.size(); ++j) {
            A[i][j]+=B._v[i]*B._w[j];
        } }
        return A;
    }
    friend Matrix<X>& operator-=(Matrix<X>& A, RankOneMatrix<X> const& B) {
        for(SizeType i=0; i!=B._v.size(); ++i) { for(SizeType j=0; j!=B._w.size(); ++j) {
            A[i][j]-=B._v[i]*B._w[j];
        } }
        return A;
    }
};
template<class X> RankOneMatrix<X> operator*(Vector<X> v, Covector<X> w) { return RankOneMatrix<X>{std::move(v),std::move(w)}; }

template<class X> class OrthogonalMatrix : public Matrix<X> {
  public:
    OrthogonalMatrix(SizeType n, PrecisionType<X> pr) : Matrix<X>(n,n,X(pr)) { }
    explicit OrthogonalMatrix(Matrix<X> const& A) : Matrix<X>(A) { }
    friend OrthogonalMatrix<X> operator*(OrthogonalMatrix<X> const& Q1, OrthogonalMatrix<X> const& Q2) {
        return static_cast<Matrix<X>const&>(Q1)*static_cast<Matrix<X>const&>(Q2); }
    friend OrthogonalMatrix<X> operator*(Transpose<OrthogonalMatrix<X>> const& Q1, OrthogonalMatrix<X> const& Q2) {
        return static_cast<Matrix<X>>(Q1)*static_cast<Matrix<X>const&>(Q2); }
    friend Matrix<X> operator*(Transpose<Matrix<X>> const& A1, OrthogonalMatrix<X> const& Q2) {
        return Matrix<X>(A1)*static_cast<Matrix<X>const&>(Q2); }
    friend Matrix<X> operator*(Matrix<X> const& A1, OrthogonalMatrix<X> const& Q2) {
        return A1*static_cast<Matrix<X>const&>(Q2); }
};


template<class X> class UpperHessenbergMatrix;

template<class X> class UpperTriangularMatrix {
    SizeType _n;
    X _z;
    Array<X> _ary;
    constexpr SizeType _pos(SizeType i, SizeType j) const { return j*(j+1)/2+i; }
    OutputStream& _write(OutputStream& os, SizeType wdth) const;
    static Matrix<X> _mul(UpperTriangularMatrix<X> const& R, Matrix<X> const& A);
    static Matrix<X> _mul(Matrix<X> const& A, UpperTriangularMatrix<X> const& R);
  public:
    UpperTriangularMatrix(SizeType n, PrecisionType<X> pr) : _n(n), _z(pr), _ary(n*(n+1)/2,_z) { }
    explicit UpperTriangularMatrix(UpperHessenbergMatrix<X> const& UH);
    SizeType row_size() const { return _n; }
    SizeType column_size() const { return _n; }
    X const& zero_element() const { return _z; }
    X& at(SizeType i, SizeType j) { ARIADNE_DEBUG_ASSERT(i<_n&&j<_n); ARIADNE_ASSERT(i<=j); return _ary[_pos(i,j)]; }
    X const& get(SizeType i, SizeType j) const { ARIADNE_DEBUG_ASSERT(i<_n&&j<_n); return (i<=j) ? _ary[_pos(i,j)] : _z; }
    Void set(SizeType i, SizeType j, X const& x) { assert(i<=j); _ary[_pos(i,j)]=x; }
    operator Matrix<X>() const { return to_matrix(*this); }
    friend Matrix<X> operator*(UpperTriangularMatrix<X> const& R, OrthogonalMatrix<X> const& Q) {
        return R*static_cast<Matrix<X>const&>(Q); }
    friend Matrix<X> operator*(Matrix<X> const& A, UpperTriangularMatrix<X> const& U) {
        return UpperTriangularMatrix<X>::_mul(A,U); }
    friend Matrix<X> operator*(UpperTriangularMatrix<X> const& U, Matrix<X> const& A) {
        return UpperTriangularMatrix<X>::_mul(U,A); }
    friend Matrix<X> operator*(Matrix<X> const&, Matrix<X> const&);
    friend OutputStream& operator<<(OutputStream& os, UpperTriangularMatrix<X> const& U) {
        return U._write(os,20); }
};


template<class X> class UpperHessenbergMatrix {
    SizeType _n;
    X _z;
    Array<X> _ary;
    constexpr SizeType _pos(SizeType i, SizeType j) const { return (j+1)*(j+2)/2-1+i; }
    static Matrix<X> _mul(Matrix<X> const& A, UpperHessenbergMatrix<X> const& UH);
    OutputStream& _write(OutputStream& os, SizeType wdth) const;
  public:
    typedef X ScalarType;
    template<class PR> requires Constructible<X,PR> UpperHessenbergMatrix(SizeType n, PR pr) : _n(n), _z(pr), _ary((n+1)*(n+2)/2-2,_z) { }
    explicit UpperHessenbergMatrix<X>(UpperTriangularMatrix<X> const& A);
    explicit UpperHessenbergMatrix<X>(Matrix<X> const& A);
    operator Matrix<X> () const { return to_matrix(*this); }
    SizeType row_size() const { return _n; }
    SizeType column_size() const { return _n; }
    X const& zero_element() const { return _z; }
    X& at(SizeType i, SizeType j) { ARIADNE_DEBUG_ASSERT(i<_n&&j<_n); ARIADNE_ASSERT(i<=j+1u); return _ary[_pos(i,j)]; }
    X const& get(SizeType i, SizeType j) const { ARIADNE_DEBUG_ASSERT(i<_n&&j<_n); return (i<=j+1u) ? _ary[_pos(i,j)] : _z; }
    Void set(SizeType i, SizeType j, X const& x) { ARIADNE_ASSERT(i<=j+1u); _ary[_pos(i,j)]=x; }
    template<class XX> requires Assignable<X,XX> Void set(SizeType i, SizeType j, XX const& x) { ARIADNE_ASSERT(i<=j+1u); _ary[_pos(i,j)]=x; }
    friend Matrix<X> operator*(Matrix<X> const& A, UpperHessenbergMatrix<X> const& UH) {
        return UpperHessenbergMatrix<X>::_mul(A,UH); }
    friend OutputStream& operator<<(OutputStream& os, UpperHessenbergMatrix<X> const& UH) {
        return UH._write(os,20); }
};


template<class X> class SymmetricTridiagonalMatrix {
    SizeType _n;
    X _z;
    Array<X> _ary;
    constexpr SizeType _pos(SizeType i, SizeType j) const { return (i==j) ? i : (i<j) ? j+_n : i+_n; }
    OutputStream& _write(OutputStream& os, SizeType wdth) const;
//    static Matrix<X> _mul(TridiagonalMatrix<X> const& R, Matrix<X> const& A);
//    static Matrix<X> _mul(Matrix<X> const& A, TridiagonalMatrix<X> const& R);
  public:
    SymmetricTridiagonalMatrix(SizeType n, PrecisionType<X> pr) : _n(n), _z(pr), _ary(3*n-1,_z) { }
    explicit SymmetricTridiagonalMatrix(UpperHessenbergMatrix<X> const& UH);
    SizeType row_size() const { return _n; }
    SizeType column_size() const { return _n; }
    X const& zero_element() const { return _z; }
    X& at(SizeType i, SizeType j) { ARIADNE_DEBUG_ASSERT(i<_n&&j<_n); ARIADNE_ASSERT(i<=j+1&&j<=i+1); return _ary[_pos(i,j)]; }
    X const& get(SizeType i, SizeType j) const { ARIADNE_DEBUG_ASSERT(i<_n&&j<_n); return (i<=j+1&&j<=i+1) ? _ary[_pos(i,j)] : _z; }
    Void set(SizeType i, SizeType j, X const& x) { ARIADNE_DEBUG_ASSERT(i<_n&&j<_n); ARIADNE_ASSERT(i<=j+1&&j<=i+1); _ary[_pos(i,j)]=x; }
    operator Matrix<X>() const { return to_matrix(*this); }
    friend Matrix<X> operator*(Matrix<X> const&, Matrix<X> const&);
    friend OutputStream& operator<<(OutputStream& os, SymmetricTridiagonalMatrix<X> const& T) {
        return Matrix<X>(T)._write(os,20); }
};


template<class X> class GivensProductMatrix {
    struct Rotation { X alpha; X beta; }; // Rotation matrix [a,-b;b,a]
    friend OutputStream& operator<<(OutputStream& os, Rotation const& G) { return os << "G{" << G.alpha << "," << G.beta << "}"; }
    Array<Rotation> _rotations;
  public:
    SizeType size() const { return _rotations.size()+1u; }
    UpperHessenbergMatrix<X> rmul(UpperTriangularMatrix<X> const& A) const;
    Matrix<X> rmul(Matrix<X> A) const;
    template<class PR> requires Constructible<X,Int,PR>
        GivensProductMatrix(SizeType n, PR pr) : _rotations(n-1,Rotation{X(1,pr),X(0,pr)}) { }
    Void set(SizeType i, X alpha, X beta) { ARIADNE_DEBUG_ASSERT(i<_rotations.size()+1u); _rotations[i]=Rotation{alpha,beta}; }
    operator OrthogonalMatrix<X> () const {
        assert(_rotations.size()!=0); Matrix<X> I=Matrix<X>::identity(this->size(),this->_rotations[0].alpha.precision()); return OrthogonalMatrix<X>(I*(*this)); }

    friend UpperHessenbergMatrix<X> operator*(UpperTriangularMatrix<X> const& R, GivensProductMatrix<X> const& Q) {
        return Q.rmul(R); }
    friend OrthogonalMatrix<X> operator*(OrthogonalMatrix<X> const& A, GivensProductMatrix<X> const& Q) {
        return OrthogonalMatrix<X>(static_cast<Matrix<X>const&>(A)*Q); }
    friend Matrix<X> operator*(Matrix<X> const& A, GivensProductMatrix<X> const& Q) {
        return Q.rmul(A); }
    friend OutputStream& operator<<(OutputStream& os, GivensProductMatrix<X> const& G) {
        return os << G._rotations; }
};

template<class X> class HouseholderProductMatrix;

template<class X> class HouseholderMatrix : public MatrixExpression<HouseholderMatrix<X>> {
  public:
    Vector<X> _u;
  public:
    typedef X ScalarType;
    HouseholderMatrix(Vector<X> v) : _u(v/sqrt(dot(v,v))) { }
    SizeType row_size() const { return _u.size(); }
    SizeType column_size() const { return _u.size(); }
    X zero_element() const { return _u.zero_element(); }
    X get(SizeType i, SizeType j) const { X res=((i==j)?1:0)+(-2*_u[i]*_u[j]); return res; }
    operator Matrix<X>() const {
        return to_matrix(*this);
        SizeType n=this->_u.size(); Matrix<X> R(n,n,this->zero_element());
        for(SizeType i=0; i!=n; ++i) { for(SizeType j=0; j!=n; ++j) { R[i][j]=this->get(i,j); } } return R; }
    friend Matrix<X> operator*(HouseholderMatrix<X> const& H, Matrix<X> A) {
        A -= (2*H._u) * (transpose(H._u)*A); return A; }
    friend Matrix<X> operator*(Matrix<X> A, HouseholderMatrix<X> const& H) {
        A -= (A*H._u) * (2*transpose(H._u)); return A; }
    friend HouseholderProductMatrix<X> operator*(HouseholderMatrix<X> const& H1, HouseholderMatrix<X> const& H2);
    friend OutputStream& operator<<(OutputStream& os, HouseholderMatrix<X> const& H) {
        return os << static_cast<Matrix<X>>(H); }
};

template<class X> class HouseholderProductMatrix {
    List<Vector<X>> _us;
    explicit HouseholderProductMatrix(List<Vector<X>> us) : _us(us) { }

    friend Matrix<X> operator*(HouseholderProductMatrix<X> const& H, Matrix<X> A) {
        auto iter=H._us.rbegin(); while(iter!=H._us.rend()) { Vector<X> const& u=*iter; A -= (2*u) * (transpose(u)*A); ++iter; } return A; }
    friend Matrix<X> operator*(Matrix<X> A, HouseholderProductMatrix<X> const& H) {
        auto iter=H._us.begin(); while(iter!=H._us.end()) { Vector<X> const& u=*iter; A -= (A*u) * (2*transpose(u)); ++iter; } return A; }

    friend Matrix<X> operator*(HouseholderProductMatrix<X> const& H, UpperHessenbergMatrix<X> UH) {
        return H*Matrix<X>(UH); }
    friend Matrix<X> operator*(UpperHessenbergMatrix<X> UH, HouseholderProductMatrix<X> const& H) {
        return Matrix<X>(UH)*H; }
    friend OrthogonalMatrix<X> operator*(HouseholderProductMatrix<X> const& H, OrthogonalMatrix<X> Q) {
        return OrthogonalMatrix<X>(H*static_cast<Matrix<X>&>(Q)); }
    friend OrthogonalMatrix<X> operator*(OrthogonalMatrix<X> Q, HouseholderProductMatrix<X> const& H) {
        return OrthogonalMatrix<X>(static_cast<Matrix<X>&>(Q)*H); }

    friend HouseholderProductMatrix<X> operator*(HouseholderProductMatrix<X> H1, HouseholderProductMatrix<X> const& H2) {
        return HouseholderProductMatrix<X>(catenate(H1._us,H2._us)); }
    friend HouseholderProductMatrix<X>& operator*=(HouseholderProductMatrix<X>& H1, HouseholderMatrix<X> const& H2) {
        H1._us.append(H2._u); return H1; }
    friend HouseholderProductMatrix<X> operator*(HouseholderProductMatrix<X> H1, HouseholderMatrix<X> const& H2) {
        H1*=H2; return H1; }
    friend HouseholderProductMatrix<X> operator*(HouseholderMatrix<X> const& H1, HouseholderMatrix<X> const& H2) {
        return HouseholderProductMatrix<X>({H1._u,H2._u}); }
  public:
    HouseholderProductMatrix(SizeType i) { };
    SizeType size() const { return _us[0].size(); }
    SizeType row_size() const { return _us[0].size(); }
    SizeType column_size() const { return _us[0].size(); }
    X zero_element() const { return _us[0].zero_element(); }
    Void append(HouseholderMatrix<X> const& H) { _us.append(H._u); }
    HouseholderProductMatrix<X> inverse() const { return HouseholderProductMatrix<X>(reverse(this->_us)); }
    friend HouseholderProductMatrix<X> inverse(HouseholderProductMatrix<X> const& H) { return H.inverse(); }
    operator OrthogonalMatrix<X>() const {
        Matrix<X> I=Matrix<X>::identity(this->size(),this->zero_element()); return OrthogonalMatrix<X>((*this)*I); }
    friend OutputStream& operator<<(OutputStream& os, HouseholderProductMatrix<X> const& H) {
        return os << "HouseholderProductMatrix("<< H._us << ")"; }
        //return os << H.operator Matrix<X>(); }
};

template<class X> class BlockDiagonalMatrix;

template<class X> Pair<OrthogonalMatrix<X>,UpperTriangularMatrix<X>> gram_schmidt(Matrix<X> A);
template<class X> Pair<X,Vector<X>> eigenvalue_solve(Matrix<X> const& A, X mu, Vector<X> v);
template<class X> Pair<HouseholderProductMatrix<X>,UpperHessenbergMatrix<X>> upper_hessenberg_factorisation(Matrix<X> A);
template<class X> Pair<GivensProductMatrix<X>,UpperTriangularMatrix<X>> hessenberg_qr_factorisation(Matrix<X> A);
template<class X> Void qr_step(UpperHessenbergMatrix<X>& H);

template<class X> Void qr_step(SymmetricTridiagonalMatrix<X>& T);


} // namespace Ariadne

#endif
