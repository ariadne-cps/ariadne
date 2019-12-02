/***************************************************************************
 *            zonotope.cpp
 *
 *  Copyright 2008--17  Alberto Casagrande, Pieter Collins
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

#include "../config.hpp"

#include "../geometry/zonotope.hpp"

#include <iostream>
#include <vector>
#include <algorithm>

#include "../utility/macros.hpp"
#include "../utility/array.hpp"
#include "../numeric/floatdp.hpp"
#include "../numeric/float_bounds.hpp"
#include "../numeric/float_error.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/diagonal_matrix.hpp"
#include "../algebra/pivot_matrix.tpl.hpp"
#include "../function/function.hpp"
    #include "../algebra/algebra.hpp"
    #include "../function/formula.hpp"
    #include "../function/taylor_model.hpp"
#include "../solvers/linear_programming.hpp"
#include "../geometry/point.hpp"
#include "../geometry/box.hpp"
#include "../geometry/list_set.hpp"
#include "../function/function.hpp"
#include "../output/geometry2d.hpp"


namespace Ariadne {

Vector<FloatDP> const& cast_raw(Vector<FloatDPValue> const& v) { return reinterpret_cast<Vector<FloatDP>const&>(v); }
Vector<FloatDP> const& cast_raw(Vector<FloatDPError> const& v) { return reinterpret_cast<Vector<FloatDP>const&>(v); }
Matrix<FloatDP> const& cast_raw(Matrix<FloatDPValue> const& v) { return reinterpret_cast<Matrix<FloatDP>const&>(v); }
Vector<FloatMP> const& cast_raw(Vector<FloatMPValue> const& v) { return reinterpret_cast<Vector<FloatMP>const&>(v); }
Vector<FloatMP> const& cast_raw(Vector<FloatMPError> const& v) { return reinterpret_cast<Vector<FloatMP>const&>(v); }
Matrix<FloatMP> const& cast_raw(Matrix<FloatMPValue> const& v) { return reinterpret_cast<Matrix<FloatMP>const&>(v); }

template<class PRE, class PR, EnableIf<IsDefaultConstructible<PRE>>> inline
    PRE _error_precision(PR const&);
template<class PRE, class PR, DisableIf<IsDefaultConstructible<PRE>>, EnableIf<IsConstructible<PRE,PR>>> inline
    PRE _error_precision(PR const& pr);

template<class X> Vector<Bounds<X>> pm(Vector<Error<X>> const& ve) { return Vector<Bounds<X>>(ve.size(),[&ve](SizeType i){return pm(ve[i]);}); }

Matrix<FloatDPApproximation>& cast_approximate(Matrix<FloatDP>& A) { return reinterpret_cast<Matrix<FloatDPApproximation>&>(A); }
Matrix<FloatDPApproximation> const& cast_approximate(Matrix<FloatDP> const& A) { return reinterpret_cast<Matrix<FloatDPApproximation>const&>(A); }

Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& e) { return reinterpret_cast<Vector<FloatDPValue>const&>(e); }
DiagonalMatrix<FloatDPValue> const& cast_exact(DiagonalMatrix<FloatDPError> const& D) { return reinterpret_cast<DiagonalMatrix<FloatDPValue>const&>(D); }

Vector<FloatMPValue> const& cast_exact(Vector<FloatMPError> const& e) { return reinterpret_cast<Vector<FloatMPValue>const&>(e); }
DiagonalMatrix<FloatMPValue> const& cast_exact(DiagonalMatrix<FloatMPError> const& D) { return reinterpret_cast<DiagonalMatrix<FloatMPValue>const&>(D); }

inline Vector<FloatDP> add(RoundApproximately,const Vector<FloatDP>& v1, const Vector<FloatDP>& v2) { return v1+v2; }
inline Vector<FloatDP> sub(RoundApproximately,const Vector<FloatDP>& v1, const Vector<FloatDP>& v2) { return v1-v2; }

Vector<FloatDP> add(RoundUpward upw, const Vector<FloatDP>& v1, const Vector<FloatDP>& v2) {
    Vector<FloatDP> result;
    for(SizeType i=0; i!=v1.size(); ++i) {
        result[i]=add(upw,v1[i],v2[i]);
    }
    return result;
}

Vector<Bounds<FloatMP>> operator+(Vector<Bounds<FloatMP>> v1, Vector<Bounds<FloatDP>> const& v2) {
    return Vector<Bounds<FloatMP>>(v1.size(),[&](SizeType i){return v1[i]+Bounds<FloatMP>(v2[i],v1[i].precision());}); }


template<class X> class LinearProgram {
  public:
    LinearProgram(const Matrix<X>& A) { ARIADNE_NOT_IMPLEMENTED; }
    Bool is_feasible() { ARIADNE_NOT_IMPLEMENTED; }
};


template<class BS>
ListSet<BS>
subdivide(const BS& bs, const FloatDP& r)
{
    ListSet<BS> result;
    ListSet<BS> working(bs);
    while(!working.size()==0) {
        BS set=working.pop();
        if(set.radius()<r) {
            result.adjoin(set);
        } else {
            working.adjoin(split(set));
        }
    }
    return result;
}

FloatDP add(FloatDP::RoundingModeType, FloatDP, FloatDP);
FloatDP sub(FloatDP::RoundingModeType, FloatDP, FloatDP);
FloatDP mul(FloatDP::RoundingModeType, FloatDP, FloatDP);


template<class X> Vector<Error<X>>
errors(const Vector<Bounds<X>>& pt)
{
    return Vector<Error<X>>(pt.size(),[&pt](SizeType i){return pt[i].error();});
}

template<class PRE, class X> Vector<FloatError<PRE>> make_errors(const Vector<Bounds<X>>& pt) {
    if constexpr (IsSame<PRE,typename X::PrecisionType>::value) { return Vector<FloatError<PRE>>(pt.size(),[&pt](SizeType i){return pt[i].error();}); }
    else { PRE pre; return Vector<FloatError<PRE>>(pt.size(),[&](SizeType i){return FloatError<PRE>(pt[i].error(),pre);}); }
}

template<class PRE, class X> FloatError<PRE> convert_error(const Error<X>& e) {
    if constexpr (IsSame<PRE,typename X::PrecisionType>::value) { return e; }
    else { PRE pre; return FloatError<PRE>(e,pre); }
}

template<class PRE, class X> Vector<FloatError<PRE>> convert_errors(const Vector<Error<X>>& pt) {
    if constexpr (IsSame<PRE,typename X::PrecisionType>::value) { return pt; }
    else { PRE pre; return Vector<FloatError<PRE>>(pt.size(),[&](SizeType i){return FloatError<PRE>(pt[i],pre);}); }
}


template<class X>
Vector<PositiveUpperBound<X>>
row_norms(const Matrix<Bounds<X>>& A)
{
    SizeType const m=A.row_size();
    SizeType const n=A.column_size();
    Vector<PositiveUpperBound<X>> e(m,mag(A.zero_element()));
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            e[i]+=mag(A[i][j]);
        }
    }
    return e;
}

template<class X> Vector<Error<X>>
row_norms(const DiagonalMatrix<Bounds<X>>& A)
{
    SizeType const n=A.row_size();
    Vector<Error<X>> e(n);
    for(SizeType i=0; i!=n; ++i) {
        e[i]=mag(A[i]);
    }
    return e;
}

template<class X> Vector<Error<X>>
row_errors(const Matrix<Bounds<X>>& A)
{
    SizeType const m=A.row_size();
    SizeType const n=A.column_size();
    Vector<Error<X>> e(m);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            e[i]+=A[i][j].error();
        }
    }
    return e;
}

template<class X> Vector<Error<X>>
row_errors(const Vector<Bounds<X>>& b, const Matrix<Bounds<X>>& A)
{
    assert(b.size()==A.row_size());
    Vector<Error<X>> result(b.size());
    for(SizeType i=0; i!=A.row_size(); ++i) {
        result[i]=b[i].error();
        for(SizeType j=0; j!=A.column_size(); ++j) {
            result[i]+=A[i][j].error();
        }
    }
    return result;
}

ValidatedKleenean
norm_grtr(const Vector<FloatDP>& v1, const Vector<FloatDP>& v2)
{
    return norm(v1)>norm(v2);
}







template<class X, class XE> Zonotope<ValidatedTag,X,XE>::~Zonotope()
{
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope()
    : _centre(), _generators(), _error()
{
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(DimensionType d)
    : _centre(d), _generators(d,0), _error(d)
{
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(DimensionType d, SizeType m)
    : _centre(d), _generators(d,m), _error(d)
{
}

template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(InitializerList< Tuple<X,InitializerList<X>,XE> > lst)
    : _centre(lst.size()), _generators(lst.size(),lst.size()==0?0u:std::get<1>(*lst.begin()).size()), _error(lst.size())
{
    for(typename InitializerList< Tuple<X,InitializerList<X>,XE> >::const_iterator aff_iter=lst.begin();
        aff_iter!=lst.end(); ++aff_iter)
    {
        SizeType i=static_cast<SizeType>(aff_iter-lst.begin()); // Cast needed to change signed to unsigned
        this->_centre[i]=static_cast<Value<X>>(std::get<0>(*aff_iter));
        for(SizeType j=0; j!=_generators.column_size(); ++j) {
            this->_generators[i][j]=static_cast<Value<X>>(*(std::get<1>(*aff_iter).begin()+j));
        }
        this->_error[i]=static_cast<Error<XE>>(std::get<2>(*aff_iter));
        ARIADNE_ASSERT(this->_error[i].raw()>=0);
    }
}



template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<X>& c, const Matrix<X>& G)
    : _centre(c), _generators(G), _error(c.size(),_error_precision<typename XE::PrecisionType>(c.zero_element().precision()))
{
    assert(c.size()==G.row_size());
}

template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<X>& c, const Matrix<X>& G, const Vector<XE>& e)
    : _centre(c), _generators(G), _error(e)
{
    assert(c.size()==G.row_size());
    assert(c.size()==e.size());
}

template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<Value<X>>& c, const Matrix<Value<X>>& G, const Vector<Error<XE>>& e)
    : _centre(cast_raw(c)), _generators(cast_raw(G)), _error(cast_raw(e))
{
    assert(c.size()==G.row_size());
    assert(c.size()==e.size());
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<Value<X>>& c, const Matrix<Value<X>>& G)
    : _centre(cast_raw(c)), _generators(cast_raw(G)), _error(c.size(),_error_precision<PRE>(c.zero_element().precision()))
{
    assert(c.size()==G.row_size());
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<Bounds<X>>& c, const Matrix<Value<X>>& G)
    : _centre(midpoint(c)), _generators(G), _error(make_errors<PRE>(c))
{
    assert(c.size()==G.row_size());
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<Value<X>>& c, const Matrix<Bounds<X>>& G)
    : _centre(c), _generators(cast_raw(midpoint(G))), _error(convert_errors<PRE>(row_errors(G)))
{
    assert(c.size()==G.row_size());
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<Bounds<X>>& c, const Matrix<Bounds<X>>& G)
    : _centre(midpoint(c)), _generators(midpoint(G)), _error(convert_errors<PRE>(row_errors(c,G)))
{
    assert(c.size()==G.row_size());
}

template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Vector<Bounds<X>>& c, const Matrix<Value<X>>& G, const Vector<Error<XE>>& e)
    : _centre(midpoint(c)), _generators(G), _error(e+make_errors<PRE>(c))
{
    assert(c.size()==G.row_size());
}



template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Point<Bounds<X>>& pt)
    : _centre(pt.dimension(),[&pt](SizeType i){return pt[i].value();}), _generators(pt.dimension(),0), _error(pt.dimension(),[&pt](SizeType i){return convert_error<PRE>(pt[i].error());})
{
}

template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Box<Interval<UpperBound<X>>>& bx)
    : _centre(bx.dimension(),[&bx](SizeType i){return cast_exact(bx[i].centre());})
    , _generators(bx.dimension(),bx.dimension())
    , _error(bx.dimension())
{
    for(SizeType i=0; i!=this->dimension(); ++i) {
        this->_generators[i][i]=cast_exact(bx[i].radius());
    }
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>::Zonotope(const Zonotope& z)
    : _centre(z._centre), _generators(z._generators), _error(z._error)
{
}



template<class X, class XE> Zonotope<ValidatedTag,X,XE>&
Zonotope<ValidatedTag,X,XE>::operator=(const Zonotope& z)
{
    if(this!=&z) {
        this->_centre=z._centre;
        this->_generators=z._generators;
        this->_error=z._error;
    }
    return *this;
}


template<class X, class XE> Zonotope<ValidatedTag,X,XE>*
Zonotope<ValidatedTag,X,XE>::clone() const
{
    return new Zonotope(*this);
}

template<class X, class XE>
Bool
operator==(const Zonotope<ValidatedTag,X,XE>& z1, const Zonotope<ValidatedTag,X,XE>& z2)
{
    return (z1._centre==z2._centre) && (z1._generators==z2._generators)
        && (cast_raw(z1._error)==cast_raw(z2._error));
}

template<class X, class XE>
DimensionType
Zonotope<ValidatedTag,X,XE>::dimension() const
{
    return this->_centre.size();
}


template<class X, class XE>
SizeType
Zonotope<ValidatedTag,X,XE>::number_of_generators() const
{
    return this->_generators.column_size();
}


template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::centre() const -> const Vector<Value<X>>&
{
    return this->_centre;
}

template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::generators() const -> const Matrix<Value<X>>&
{
    return this->_generators;
}

template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::error() const -> const Vector<Error<XE>>&
{
    return this->_error;
}


template<class X, class XE>
Box<ExactIntervalType>
Zonotope<ValidatedTag,X,XE>::domain() const
{
    return Box<ExactIntervalType>(this->number_of_generators(),ExactIntervalType(-1,1));
}


template<class X, class XE>
UpperBoxType
Zonotope<ValidatedTag,X,XE>::bounding_box() const
{
    const Zonotope& z=*this;
    Vector<Bounds<X>> d(cast_singleton(z.domain()),z.centre().zero_element().precision());
    UpperBoxType b=UpperBoxType(Vector<Bounds<FloatDP>>(Vector<Bounds<X>>(z.centre()+z.generators()*d+pm(z.error())),double_precision));
    return b;
}


template<class X, class XE>
PositiveFloatDPUpperBound
Zonotope<ValidatedTag,X,XE>::radius() const
{
    return Ariadne::radius(this->bounding_box());
}

template<class X, class XE>
UpperBoxType
bounding_box(const Zonotope<ValidatedTag,X,XE>& z)
{
    return z.bounding_box();
}













template<class X> decltype(auto) properties(X const& x) { return x.precision(); }

template<class X> Matrix<X> tensor_product(const Matrix<X>& A1, const Matrix<X>& A2) {
    SizeType m1=A1.row_size(); SizeType m2=A2.row_size(); SizeType n1=A1.column_size(); SizeType n2=A2.column_size();
    X z=max(A1.zero_element(),A2.zero_element());
    Matrix<X> R(m1+m2,n1+n2,z);
    R[range(0,m1)][range(0,n1)]=A1;
    R[range(m1,m1+m2)][range(n1,n1+n2)]=A2;
    return R;
}


template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::_product(const Zonotope<P,X,XE>& z1, const Zonotope<P,X,XE>& z2) -> Zonotope<P,X,XE>
{
    return Zonotope(join(z1.centre(),z2.centre()),tensor_product(z1.generators(),z2.generators()),join(z1.error(),z2.error()));
}

template<class X, class XE> Zonotope<ValidatedTag,X,XE>
Zonotope<ValidatedTag,X,XE>::_product(const Zonotope& z1, const Interval<UpperBound<X>>& ivl2)
{
    return _product(z1,Box<Interval<UpperBound<X>>>(1u,ivl2));
}

template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::_split(const Zonotope<P,X,XE>& z) -> ListSet< Zonotope<P,X,XE> >
{
    ListSet< Zonotope<P,X,XE>  > result;

    auto pr=z.centre().zero_element().precision();
    DimensionType d=z.dimension();
    SizeType m=z.number_of_generators();
    Vector<Value<X>> const& c=z.centre();
    Matrix<Value<X>> const& G=z.generators();
    Vector<Error<XE>> const& e=z.error();

    Array<PositiveUpperBound<X>> norms(m,PositiveUpperBound<X>(pr));
    for(SizeType j=0; j!=m; ++j) {
        norms[j]=norm(Vector<Value<X>>(column(G,j)));
    }

    PositiveUpperBound<X> max_norm(pr);
    SizeType longest_generator=0;
    for(SizeType j=0; j<m; ++j) {
        if(norms[j].raw()>max_norm.raw()) {
            max_norm=norms[j];
            longest_generator=j;
        }
    }
    for(SizeType k=0; k<d; ++k) {
        if(e[k].raw()>max_norm.raw()) {
            max_norm=Error<X>(e[k],pr);
            longest_generator=m+k;
        }
    }

    if(longest_generator<m) {
        Matrix<Value<X>> new_generators=z.generators();
        SizeType j=longest_generator;
        for(SizeType i=0; i!=d; ++i) {
            new_generators[i][j]=hlf(new_generators[i][j]);
        }

        Vector<Value<X>> v=column(new_generators,j);
        Vector<Bounds<X>> new_centre=c-v;
        result.adjoin(Zonotope(new_centre,new_generators,e));
        new_centre=c+v;
        result.adjoin(Zonotope(new_centre,new_generators,e));
    } else {
        SizeType k=longest_generator-m;
        Vector<Bounds<X>> new_centre = z.centre();
        const Matrix<Value<X>>& new_generators = z.generators();
        Vector<Error<XE>> new_error=e;
        new_error[k]=hlf(new_error[k]);
        Vector<Value<X>> new_exact_error=cast_exact(Vector<Error<X>>(new_error,pr));
        new_centre[k]=z.centre()[k]+new_exact_error[k];
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
        new_centre[k]=z.centre()[k]-new_exact_error[k];
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
    }
    return result;
}











template<class X, class XE> auto
approximation(const ValidatedZonotope<X,XE>& z) -> ValidatedZonotope<X,XE>
{
    return Zonotope(z.centre(),z.generators(),z.error());
}


template<class X, class XE> auto
over_approximation(const ValidatedZonotope<X,XE>& z) -> ValidatedZonotope<X,XE>
{
    return z;
}



template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::_error_free_over_approximation(const Zonotope<P,X,XE>& z) -> Zonotope<P,X,XE>
{
    SizeType d=z.dimension();
    SizeType m=z.number_of_generators();

    // Count number of nonzero error values
    SizeType e=0;
    for(SizeType i=0; i!=d; ++i) {
        if(possibly(z.error()[i]>0)) { ++e; }
    }
    Matrix<Value<X>> nG(d,m+e);
    nG[range(0,d)][range(0,m)]=z.generators();

    SizeType j=m;
    for(SizeType i=0; i!=d; ++i) {
        if(possibly(z.error()[i]>0)) {
            nG[i][j]=cast_exact(Error<X>(z.error()[i],z.centre().zero_element().precision()));
            ++j;
        }
    }
    return Zonotope(z.centre(),nG);
}





template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::_orthogonal_over_approximation(const Zonotope<P,X,XE>& z) -> Zonotope<P,X,XE>
{
    //assert(iz.size()==iz.number_of_generators());
    Zonotope ez=error_free_over_approximation(z);

    const Vector<Value<X>>& c=ez.centre();
    const Matrix<Value<X>>& G=ez.generators();

    const Matrix<Approximation<X>> aG=G;
    PivotMatrix P;
    Matrix<Approximation<X>> aQ,aR;
    make_ltuple(aQ,aR,P)=orthogonal_decomposition(aG,true);
    Matrix<Value<X>> Q=cast_exact(aQ);
    Matrix<Bounds<X>> Qinv=inverse(Q); // Don't use transpose as Q is not exactly orthogonal!
    Matrix<Bounds<X>> iR=Qinv*G;
    DiagonalMatrix<Value<X>> D(cast_exact(row_norms(iR)));

    Matrix<Bounds<X>> nG=Q*D;

    return Zonotope(c,nG);
}

template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::_cascade_over_approximation(const Zonotope<P,X,XE>& z, SizeType cs) -> Zonotope<P,X,XE>
{
    using namespace std;

    if(z.number_of_generators()<=z.dimension()*cs) { return z; }

    assert(z.number_of_generators() % z.dimension()==0);

    SizeType d=z.dimension();
    SizeType nb=z.number_of_generators()/z.dimension(); // number of generator blocks

    const Matrix<Value<X>>& G=z.generators();
    Array<PositiveUpperBound<X>> norms(nb);
    for(SizeType i=0; i!=nb; ++i) {
        Matrix<Value<X>> PG=Matrix<Value<X>>(G[range(0,d)][range(i*d,(i+1)*d)]);
        norms[i]=norm(PG);
    }

    // Compute the new number of blocks
    SizeType nnb=cs;
    PositiveUpperBound<X> sum(z.centre().zero_element().precision());
    for(SizeType i=nb-1; i!=0; --i) {
        sum+=norms[i];
        if(sum.raw()>norms[i-1].raw()) {
            nnb=i;
        }
    }
    nnb=min(nnb,cs);
    // Reduce generators
    // TODO: Is this correct?
    Matrix<Value<X>> rG(d,d*nnb);
    rG[range(0,d)][range(0,d*(nnb-1))]=G[range(0,d)][range(0,d*(nnb-1))];
    for(SizeType i=0; i!=d; ++i) {
        UpperBound<X> err=rG[i][d*(nnb-1)+i];
        for(SizeType j=d*(nnb-1); j!=G.column_size(); ++j) {
            err+=abs(G[i][j]);
        }
        rG[i][d*(nnb-1)+i]=cast_exact(err);
    }
    return Zonotope(z.centre(),rG);
}


/*
Zonotope
orthogonal_approximation(const Zonotope& z)
{

    Vector<X> c=z.centre();
    Matrix<X> J=z.generators();
    Vector<X> e=z.error();

    const SizeType m=J.row_size();
    const SizeType n=J.column_size();

    Matrix<X> G(m,m+m);

    Matrix<X> Q;
    Matrix<X> R;
    make_ltuple(cast_approximate(Q),cast_approximate(R))=orthogonal_decomposition(cast_approximate(J));

    ARIADNE_ASSERT(norm(FloatMatrix(Q*R-J))<1e-8);

    for(SizeType i=0; i!=m;++i) {
        X a=0;
        for(SizeType j=i; j!=n; ++j) {
            a+=abs(R[i][j]);
        }
        for(SizeType k=0; k!=m; ++k) {
            X b=Q[k][i]*a;
            G[k][i]=b;
        }
    }

    for(SizeType i=0; i!=m; ++i) { G[i][m+i]=e[i]; }
    return Zonotope(c,G);
}
*/

/*

Zonotope<ExactIntervalType,R>
orthogonal_over_approximation(const Zonotope<R,R>& z)
{
    // FIXME: Subdivide in zero order as well!
    static Bool warn=true;
    if(warn) {
        ARIADNE_WARN("orthogonal_over_approximation(Zonotope<I,R>) does not over-approximate roundoff errors.\n);
        warn=false;
    }
    Zonotope<R,R> oaz=over_approximation(z);

    QRMatrix< ExactIntervalType > QR(oaz.generators());
    ExactPoint< ExactIntervalType > c(oaz.centre());
    Matrix<X> G(z.size(),z.number_of_generators());

    Matrix< ExactIntervalType > q=QR.Q();
    Matrix< ExactIntervalType > r=QR.R();
    for(SizeType i=0; i!=z.size();++i) {
        ExactIntervalType a=0;
        for(SizeType j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(SizeType k=0; k!=z.size(); ++k) {
            ExactIntervalType b=q(k,i)*a;
            G(k,i)=b.midpoint();
            c[k]+=(b-b.midpoint());
        }
    }
    return Zonotope<R,R>(midpoint(c),G);
}


Zonotope<ExactIntervalType,R>
orthogonal_over_approximation(const Zonotope<ExactIntervalType,R>& z)
{
    Zonotope<R,R> oaz=over_approximation(z);

    QRMatrix< ExactIntervalType > QR(oaz.generators());
    ExactPoint< ExactIntervalType > c(oaz.centre());
    Matrix<X> G(z.size(),z.number_of_generators());

    Matrix< ExactIntervalType > q=QR.Q();
    Matrix< ExactIntervalType > r=QR.R();
    for(SizeType i=0; i!=z.size();++i) {
        ExactIntervalType a=0;
        for(SizeType j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(SizeType k=0; k!=z.size(); ++k) {
            ExactIntervalType b=q(k,i)*a;
            G(k,i)=b.midpoint();
            c[k]+=(b-b.midpoint());
        }
    }
    return Zonotope<ExactIntervalType,R>(c,G);
}


Zonotope< ExactIntervalType >
orthogonal_over_approximation(const Zonotope< ExactIntervalType >& z)
{
    Zonotope<R,R> oaz=over_approximation(z);

    QRMatrix< ExactIntervalType > QR(oaz.generators());
    ExactPoint< ExactIntervalType > c(oaz.centre());
    Matrix<X> G(z.size(),z.number_of_generators());

    Matrix< ExactIntervalType > q=QR.Q();
    Matrix< ExactIntervalType > r=QR.R();
    for(SizeType i=0; i!=z.size();++i) {
        ExactIntervalType a=0;
        for(SizeType j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(SizeType k=0; k!=z.size(); ++k) {
            ExactIntervalType b=q(k,i)*a;
            G(k,i)=b.midpoint();
            c[k]+=(b-b.midpoint());
        }
    }
    return Zonotope< ExactIntervalType >(c,G);
}
*/

template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::_apply(const ValidatedVectorMultivariateFunction& f, const ValidatedZonotope<X,XE>& z) -> ValidatedZonotope<X,XE>
{
    PR pr=z.centre().zero_element().precision();
    PRE pre=_error_precision<PRE>(pr);
    Vector<Bounds<X>> zc=z.centre();
    Matrix<Value<X>> zG=z.generators();
    Vector<Bounds<X>> ze(pm(z.error()),pr);
#warning
    Vector<Bounds<X>> zb=Vector<Bounds<X>>(cast_singleton(z.bounding_box()),pr);

    Vector<Bounds<X>> fc=f(zc);
    Matrix<Bounds<X>> fJb=jacobian(f,zb);

    Matrix<Bounds<X>> fJbzG=fJb*zG;

    Vector<X> nzc = cast_raw(midpoint(fc));
    Matrix<X> nzG = cast_raw(midpoint(fJbzG));

    Vector<Bounds<X>> zE(z.number_of_generators(),Bounds<X>(-1,+1,pr));

    Vector<Bounds<X>> nzE=(fc-Vector<Bounds<X>>(nzc)) + (fJbzG-Matrix<Bounds<X>>(nzG))*zE + fJb*ze;

    Vector<XE> nze(nzE.size()); for(SizeType i=0; i!=nze.size(); ++i) { nze[i]=XE(nzE[i].upper().raw(),up,pre); }

    return ValidatedZonotope<X,XE>(nzc,nzG,nze);
}




template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::write(OutputStream& os) const -> OutputStream&
{
    const Zonotope& z = *this;

    os << "[";
    for(SizeType i=0; i!=z.dimension(); ++i) {
            os << (i==0 ? '(' : ',') << z.centre()[i];
    }
    os << "),";
    for(SizeType j=0; j!=z.number_of_generators(); ++j) {
        for(SizeType i=0; i!=z.dimension(); ++i) {
            os << (i==0 ? '[' : ',') << z.generators()[i][j];
        }
        os << "],";
    }
    for(SizeType i=0; i!=z.dimension(); ++i) {
        os << (i==0 ? '{' : ',') << z.error()[i];
    }
    os << '}';
    os << "]";
    return os;
}






template<> inline double numeric_cast<double,FloatMP>(FloatMP const& x) { return x.get_d(); }


template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::inside(const ExactBoxType& bx) const -> ValidatedLowerKleenean
{
    return this->bounding_box().inside(bx);
}

Box<Interval<FloatDPUpperBound>> convert_box(ExactBoxType bx, DoublePrecision pr) { return bx; }

Box<Interval<FloatMPUpperBound>> convert_box(ExactBoxType bx, MultiplePrecision pr) {
    typedef MultiplePrecision PR;
    typedef RawFloatType<PR> X;
    typedef Interval<UpperBound<X>> UIvl;
//    return Box<UIvl>(bx.dimension(),[&bx,&pr](SizeType i){return Interval<FloatUpperBound<PR>>(FloatLowerBound<PR>(bx[i],pr);});
    return Box<UIvl>(bx.dimension(),
        [&bx,&pr](SizeType i){return Interval<FloatUpperBound<PR>>(X(bx[i].lower().raw(),pr),X(bx[i].upper().raw(),pr));});
}

Interval<FloatMPUpperBound> add(Interval<FloatMPUpperBound> const& ivl1, Interval<FloatMPUpperBound> const& ivl2) {
    return Interval<FloatMPUpperBound>(cast_singleton(ivl1)+cast_singleton(ivl2)); }


/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[I,z.G], b=z.c, l=[r.l,-o], u=[r.u,+o]
 */
template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::separated(const ExactBoxType& bx) const -> ValidatedLowerKleenean
{
    const Zonotope& z = *this;
    PR pr=z.centre().zero_element().precision();
    ARIADNE_ASSERT(z.dimension()==bx.dimension());
    SizeType d=z.dimension();
    SizeType ng=z.number_of_generators();
    Vector<Bounds<X>> ebx=cast_singleton(convert_box(bx,pr))+pm(z.error());
    const Vector<Value<X>>& zc=z.centre();
    const Matrix<Value<X>>& zG=z.generators();
    Matrix<Value<X>> A(d,d+ng);
    Vector<Value<X>> b(d);
    Vector<Value<X>> xl(d+ng);
    Vector<Value<X>> xu(d+ng);

    A[range(0,d)][range(0,d)]=Matrix<Value<X>>::identity(d);
    A[range(0,d)][range(d,d+ng)]=zG;
    b=zc;
    for(SizeType j=0; j!=d; ++j) {
        xl[j]=cast_exact(ebx[j].lower());
        xu[j]=cast_exact(ebx[j].upper());
    }
    for(SizeType j=0; j!=ng; ++j) {
        xl[d+j]=-1;
        xu[d+j]=+1;
    }

    return ! SimplexSolver<Value<X>>().feasible(xl,xu,A,b);
}



/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[z1.G,z2.G], b=z1.c-z2.c, l=[-o,-o], u=[+o,+o]
 * Still need to take into account errors, particularly in
 * the \a b vector.
 */
template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::_separated(const Zonotope& z1, const Zonotope& z2) -> ValidatedLowerKleenean
{
    PR pr=min(z1.centre().zero_element().precision(),z2.centre().zero_element().precision());
    ARIADNE_ASSERT(z1.dimension()==z2.dimension());
    SizeType d=z1.dimension();
    SizeType ng1=z1.number_of_generators();
    SizeType ng2=z2.number_of_generators();
    const Vector<Value<X>>& c1=z1.centre();
    const Matrix<Value<X>>& G1=z1.generators();
    const Vector<Value<X>>& c2=z2.centre();
    const Matrix<Value<X>>& G2=z2.generators();

    Matrix<Value<X>> A(d,ng1+ng2,pr);
    Vector<Value<X>> b(cast_exact(c1-c2));
    Vector<Value<X>> xl(ng1+ng2,Value<X>(-1,pr));
    Vector<Value<X>> xu(ng1+ng2,Value<X>(+1,pr));

    A[range(0,d)][range(0,ng1)]=G1;
    A[range(0,d)][range(ng1,ng1+ng2)]=G2;

    return ! SimplexSolver<Value<X>>().feasible(xl,xu,A,b);
}



void bar() { Value<FloatMP>(Dyadic(Value<FloatDP>(0)),MultiplePrecision(12)); }

/* Set up LP problem to solve \f$c+Gs=p\f$; \f$-1<=s<=1\f$.
 */
template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::contains(const ExactPoint& pt) const -> ValidatedLowerKleenean
{
    PR pr=this->centre().zero_element().precision();
    Vector<Bounds<X>> vpt(Vector<Bounds<X>>(pt.dimension(),[&](SizeType i){return Bounds<X>(Dyadic(pt[i]),pr);}));

    //std::clog << "Zonotope<ValidatedTag,X,XE>::contains(const Vector<X>& )" << std::endl;
    const Zonotope& z = *this;
    assert(z.dimension()==pt.dimension());
    DimensionType d=z.dimension();
    SizeType m=z.number_of_generators();

    Matrix<Value<X>> A(d,m+d,z.generators().zero_element());
    A[range(0,d)][range(0,m)]=z.generators();
    Vector<Bounds<X>> r=(vpt-z.centre());
    Vector<Value<X>> b(d,z.centre().zero_element());
    for(SizeType i=0; i!=d; ++i) { b[i]=r[i].value(); A[i][m+i]=cast_exact(r[i].error()+Error<X>(z.error()[i],pr)); }
    Vector<Value<X>> xl(m, -1,pr);
    Vector<Value<X>> xu(m, +1,pr);

    ValidatedKleenean result=SimplexSolver<Value<X>>().feasible(xl,xu,A,b);
    return result;
}

struct angle_less {
    Bool operator() (const Vector2d& v1, const Vector2d& v2) const {
        assert(v1.x>0 && v2.x>0);
        return (v1.y/v1.x) < (v2.y/v2.x);
    }
};

template<class X, class XE> auto
Zonotope<ValidatedTag,X,XE>::draw(CanvasInterface& c, const Projection2d& p) const -> Void
{
    const Zonotope& z=*this;
    SizeType ix=p.x_coordinate(); SizeType iy=p.y_coordinate();

    const Vector<X>& zc=cast_raw(z.centre());
    const Matrix<X>& zg=cast_raw(z.generators());
    const Vector<XE>& ze=cast_raw(z.error());

    double eps=1.0/(1ul<<31);
    Point2d pc(zc[ix],zc[iy]);
    std::vector< Vector2d > pg;
    for(SizeType j=0; j!=z.number_of_generators(); ++j) {
        Vector2d g(zg[ix][j],zg[iy][j]);
        if(g.x<0) { g=-g; }
        else if (g.x==0) { g.x=eps; }
        pg.push_back(g);
    }
    if(ze[ix]>0) { pg.push_back(Vector2d(ze[ix],0.0)); }
    if(ze[iy]>0) { pg.push_back(Vector2d(eps,ze[iy])); }

    std::sort(pg.begin(),pg.end(),angle_less());

    const SizeType npg=pg.size();
    Point2d pt=pc;
    for(SizeType i=0; i!=npg; ++i) {
        pt-=pg[i];
    }

    c.move_to(pt.x,pt.y);
    for(SizeType i=0; i!=npg; ++i) {
        pt+=2*pg[i];
        c.line_to(pt.x,pt.y);
    }
    for(SizeType i=0; i!=npg; ++i) {
        pt-=2*pg[i];
        c.line_to(pt.x,pt.y);
    }
    c.fill();
}

template class Zonotope<ValidatedTag,FloatDP,FloatDP>;
template class Zonotope<ValidatedTag,FloatMP,FloatMP>;
template class Zonotope<ValidatedTag,FloatMP,FloatDP>;

} // namespace Ariadne

