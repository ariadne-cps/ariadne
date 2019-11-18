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
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/diagonal_matrix.hpp"
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

Matrix<FloatDPApproximation>& cast_approximate(Matrix<FloatDP>& A) { return reinterpret_cast<Matrix<FloatDPApproximation>&>(A); }
Matrix<FloatDPApproximation> const& cast_approximate(Matrix<FloatDP> const& A) { return reinterpret_cast<Matrix<FloatDPApproximation>const&>(A); }

DiagonalMatrix<FloatDPValue> const& cast_exact(DiagonalMatrix<FloatDPError> const& D) { return reinterpret_cast<DiagonalMatrix<FloatDPValue>const&>(D); }

inline Vector<FloatDP> add(RoundApproximately,const Vector<FloatDP>& v1, const Vector<FloatDP>& v2) { return v1+v2; }
inline Vector<FloatDP> sub(RoundApproximately,const Vector<FloatDP>& v1, const Vector<FloatDP>& v2) { return v1-v2; }

Vector<FloatDP> add(RoundUpward upw, const Vector<FloatDP>& v1, const Vector<FloatDP>& v2) {
    Vector<FloatDP> result;
    for(SizeType i=0; i!=v1.size(); ++i) {
        result[i]=add(upw,v1[i],v2[i]);
    }
    return result;
}


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

Void
accumulate(FloatDP& value, FloatDP& error, SizeType n, const FloatDP* aptr, const FloatDP* bptr)
{
    FloatDP vl=value.raw();
    FloatDP vu=value.raw();
    for(SizeType i=0; i!=n; ++i) {
        vl=add(down,vl,mul(down,aptr[i],bptr[i]));
        vu=add(up,vu,mul(up,aptr[i],bptr[i]));
    }
    value=hlf(add(near,vl,vu));
    error=max(sub(up,vu,value),sub(up,value,vl));
}

Vector<PositiveFloatDPUpperBound>
row_norms(const Matrix<FloatDPBounds>& A)
{
    SizeType const m=A.row_size();
    SizeType const n=A.column_size();
    Vector<PositiveFloatDPUpperBound> e(m);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            e[i]+=mag(A[i][j]);
        }
    }
    return e;
}

Vector<FloatDPError>
row_norms(const DiagonalMatrix<FloatDPBounds>& A)
{
    SizeType const n=A.row_size();
    Vector<FloatDPError> e(n);
    for(SizeType i=0; i!=n; ++i) {
        e[i]=mag(A[i]);
    }
    return e;
}

Vector<FloatDPError>
row_errors(const Matrix<FloatDPBounds>& A)
{
    SizeType const m=A.row_size();
    SizeType const n=A.column_size();
    Vector<FloatDPError> e(m);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            e[i]+=A[i][j].error();
        }
    }
    return e;
}

Vector<FloatDPError>
errors(const Vector<FloatDPBounds>& pt)
{
    Vector<FloatDPError> result(pt.size());
    for(SizeType i=0; i!=pt.size(); ++i) {
        result[i]=pt[i].error();
    }
    return result;
}


Vector<FloatDPError>
row_errors(const Vector<FloatDPBounds>& b, const Matrix<FloatDPBounds>& A)
{
    assert(b.size()==A.row_size());
    Vector<FloatDPError> result(b.size());
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







Zonotope::~Zonotope()
{
}


Zonotope::Zonotope()
    : _centre(), _generators(), _error()
{
}


Zonotope::Zonotope(DimensionType d)
    : _centre(d), _generators(d,0), _error(d)
{
}


Zonotope::Zonotope(DimensionType d, SizeType m)
    : _centre(d), _generators(d,m), _error(d)
{
}

Zonotope::Zonotope(InitializerList< Tuple<FloatDP,InitializerList<FloatDP>,FloatDP> > lst)
    : _centre(lst.size()), _generators(lst.size(),lst.size()==0?0u:std::get<1>(*lst.begin()).size()), _error(lst.size())
{
    for(InitializerList< Tuple<FloatDP,InitializerList<FloatDP>,FloatDP> >::const_iterator aff_iter=lst.begin();
        aff_iter!=lst.end(); ++aff_iter)
    {
        SizeType i=static_cast<SizeType>(aff_iter-lst.begin()); // Cast needed to change signed to unsigned
        this->_centre[i]=std::get<0>(*aff_iter);
        for(SizeType j=0; j!=_generators.column_size(); ++j) {
            this->_generators[i][j]=*(std::get<1>(*aff_iter).begin()+j);
        }
        this->_error[i]=std::get<2>(*aff_iter);
        ARIADNE_ASSERT(this->_error[i]>=0);
    }
}


Zonotope::Zonotope(const Vector<FloatDP>& c, const Matrix<FloatDP>& G)
    : _centre(c), _generators(G), _error(c.size(),c.zero_element().precision())
{
    assert(c.size()==G.row_size());
}

Zonotope::Zonotope(const Vector<FloatDP>& c, const Matrix<FloatDP>& G, const Vector<FloatDP>& e)
    : _centre(c), _generators(G), _error(e)
{
    assert(c.size()==G.row_size());
    assert(c.size()==e.size());
}

Zonotope::Zonotope(const Vector<FloatDPValue>& c, const Matrix<FloatDPValue>& G, const Vector<FloatDPError>& e)
    : _centre(cast_raw(c)), _generators(cast_raw(G)), _error(cast_raw(e))
{
    assert(c.size()==G.row_size());
    assert(c.size()==e.size());
}


Zonotope::Zonotope(const Vector<FloatDPValue>& c, const Matrix<FloatDPValue>& G)
    : _centre(cast_raw(c)), _generators(cast_raw(G)), _error(c.size(),FloatDP(double_precision))
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<FloatDPBounds>& c, const Matrix<FloatDPValue>& G)
    : _centre(cast_raw(midpoint(c))), _generators(cast_raw(G)), _error(cast_raw(errors(c)))
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<FloatDPValue>& c, const Matrix<FloatDPBounds>& G)
    : _centre(cast_raw(c)), _generators(cast_raw(midpoint(G))), _error(cast_raw(row_errors(G)))
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<FloatDPBounds>& c, const Matrix<FloatDPBounds>& G)
    : _centre(cast_raw(midpoint(c))), _generators(cast_raw(midpoint(G))), _error(cast_raw(row_errors(c,G)))
{
    assert(c.size()==G.row_size());
}



Zonotope::Zonotope(const Zonotope& z)
    : _centre(z._centre), _generators(z._generators), _error(z._error)
{
}


Zonotope&
Zonotope::operator=(const Zonotope& z)
{
    if(this!=&z) {
        this->_centre=z._centre;
        this->_generators=z._generators;
        this->_error=z._error;
    }
    return *this;
}


Zonotope*
Zonotope::clone() const
{
    return new Zonotope(*this);
}

Bool
operator==(const Zonotope& z1, const Zonotope& z2)
{
    return (z1._centre==z2._centre) && (z1._generators==z2._generators)
        && (z1._error==z2._error);
}

DimensionType
Zonotope::dimension() const
{
    return this->_centre.size();
}


SizeType
Zonotope::number_of_generators() const
{
    return this->_generators.column_size();
}


const Vector<FloatDP>&
Zonotope::centre() const
{
    return this->_centre;
}


const Matrix<FloatDP>&
Zonotope::generators() const
{
    return this->_generators;
}


const Vector<FloatDP>&
Zonotope::error() const
{
    return this->_error;
}


Vector<ExactIntervalType>
Zonotope::domain() const
{
    return Vector<ExactIntervalType>(this->number_of_generators(),ExactIntervalType(-1,1));
}


UpperBoxType
Zonotope::bounding_box() const
{
    const Zonotope& z=*this;
//    std::cerr<<"zD="<<z.domain()<<"\n";
//    std::cerr<<"zc="<<cast_exact(z.centre())<<"\n";
//    std::cerr<<"zG="<<cast_exact(z.generators())<<"\n";
//    std::cerr<<"ze"<<cast_exact(z.error())*ExactIntervalType(-1,1)<<"\n";
//    std::cerr<<"zG*E"<<cast_exact(z.generators())*z.domain()<<"\n";
    UpperBoxType b=cast_exact(z.centre())+(cast_exact(z.generators())*z.domain())+cast_exact(z.error())*ExactIntervalType(-1,1);
//    std::cerr<<"bb="<<b<<"\n";
    return b;
}


PositiveFloatDPUpperBound
Zonotope::radius() const
{
    return Ariadne::radius(this->bounding_box());
}

ValidatedKleenean
Zonotope::contains(const ExactPoint& pt) const
{
    return Ariadne::contains(*this,pt);
}

ValidatedLowerKleenean
Zonotope::separated(const ExactBoxType& bx) const
{
    return Ariadne::separated(*this,ExactBoxType(bx));
}


ValidatedLowerKleenean
Zonotope::inside(const ExactBoxType& bx) const
{
    return Ariadne::inside(*this,ExactBoxType(bx));
}



OutputStream&
Zonotope::write(OutputStream& os) const
{
    return os << *this;
}





ValidatedKleenean
empty(const Zonotope& z)
{
    return false;
}


ValidatedKleenean
is_bounded(const Zonotope& z)
{
    return true;
}




PositiveFloatDPUpperBound
radius(const Zonotope& z)
{
    return Ariadne::radius(z.bounding_box());
}










UpperBoxType
bounding_box(const Zonotope& z)
{
    return z.bounding_box();
}



ListSet< Zonotope >
split(const Zonotope& z)
{
    // FIXME: Not quite guarenteed to give an over-approximation
    ListSet< Zonotope  > result;

    DimensionType d=z.dimension();
    SizeType m=z.number_of_generators();
    Vector<FloatDP> const& c=z.centre();
    Matrix<FloatDP> const& G=z.generators();
    Vector<FloatDP> const& e=z.error();

    Array<FloatDP> norms(m,0);
    for(SizeType j=0; j!=m; ++j) {
        norms[j]=norm(Vector<FloatDP>(column(G,j)));
    }

    FloatDP max_norm=0;
    SizeType longest_generator=0;
    for(SizeType j=0; j<m; ++j) {
        if(norms[j]>max_norm) {
            max_norm=norms[j];
            longest_generator=j;
        }
    }
    for(SizeType k=0; k<d; ++k) {
        if(e[k]>max_norm) {
            max_norm=e[k];
            longest_generator=m+k;
        }
    }

    if(longest_generator<m) {
        Matrix<FloatDP> new_generators=z.generators();
        SizeType j=longest_generator;
        for(SizeType i=0; i!=d; ++i) {
            new_generators[i][j]=div(up,new_generators[i][j],2);
        }

        Vector<FloatDP> v=column(new_generators,j);
        Vector<FloatDP> new_centre=sub(approx,c,v);
        result.adjoin(Zonotope(new_centre,new_generators,e));
        new_centre=add(approx,c,v);
        result.adjoin(Zonotope(new_centre,new_generators,e));
    } else {
        SizeType k=longest_generator-m;
        Vector<FloatDP> new_centre = z.centre();
        const Matrix<FloatDP>& new_generators = z.generators();
        Vector<FloatDP> new_error=e;
        new_error[k]=div(up,new_error[k],2);
        new_centre[k]=add(approx,z.centre()[k],new_error[k]);
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
        new_centre[k]=sub(approx,z.centre()[k],new_error[k]);
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
    }
    return result;
}





Zonotope::Zonotope(const ExactBoxType& r)
    : _centre(r.size()), _generators(r.size(),r.size()), _error(r.size())
{
    SizeType d=r.size();
    Vector<FloatDP>& c=this->_centre;
    Matrix<FloatDP>& G=this->_generators;
    Vector<FloatDP>& e=this->_error;
    for(SizeType i=0; i!=d; ++i) {
        c[i]=med(approx,r[i].lower().raw(),r[i].upper().raw());
        for(SizeType j=0; j!=d; ++j) {
            G[i][j]=0;
        }
        G[i][i]=rad(up,r[i].lower().raw(),r[i].upper().raw());
        e[i]=0;
    }
}







Zonotope
approximation(const Zonotope& z)
{
    return Zonotope(z.centre(),z.generators(),z.error());
}


Zonotope
over_approximation(const Zonotope& z)
{
    return z;
}


Zonotope
error_free_over_approximation(const Zonotope& z)
{
    SizeType d=z.dimension();
    SizeType m=z.number_of_generators();

    // Count number of nonzero error values
    SizeType e=0;
    for(SizeType i=0; i!=d; ++i) {
        if(z.error()[i]!=0) { ++e; }
    }
    Matrix<FloatDP> nG(d,m+e);
    project(nG,range(0,d),range(0,m))=z.generators();

    SizeType j=m;
    for(SizeType i=0; i!=d; ++i) {
        if(z.error()[i]!=0) {
            nG[i][j]=z.error()[i];
            ++j;
        }
    }
    return Zonotope(z.centre(),nG);
}





Zonotope
orthogonal_over_approximation(const Zonotope& z)
{
    //assert(iz.size()==iz.number_of_generators());
    Zonotope ez=error_free_over_approximation(z);

    const Vector<FloatDPValue>& c=cast_exact(ez.centre());
    const Matrix<FloatDPValue>& G=cast_exact(ez.generators());

    const Matrix<FloatDPApproximation> aG=cast_approximate(G);
    Matrix<FloatDPApproximation> aQ,aR;
    make_ltuple(aQ,aR)=orthogonal_decomposition(aG);
    Matrix<FloatDPValue> Q=cast_exact(Q);
    Matrix<FloatDPValue> Qinv=transpose(Q);

    Matrix<FloatDPBounds> iR=Qinv*G;
    DiagonalMatrix<FloatDPValue> D(cast_exact(row_norms(iR)));

    Matrix<FloatDPBounds> nG=Q*D;

    return Zonotope(c,nG);
}

Zonotope
cascade_over_approximation(const Zonotope& z, SizeType cs)
{
    using namespace std;

    if(z.number_of_generators()<=z.dimension()*cs) { return z; }

    assert(z.number_of_generators() % z.dimension()==0);

    SizeType d=z.dimension();
    SizeType nb=z.number_of_generators()/z.dimension(); // number of generator blocks

    const Matrix<FloatDP>& G=z.generators();
    Array<FloatDP> norms(nb);
    for(SizeType i=0; i!=nb; ++i) {
        Matrix<FloatDP> PG=Matrix<FloatDP>(project(G,range(0,d),range(i*d,(i+1)*d)));
        norms[i]=norm(PG);
    }

    // Compute the new number of blocks
    SizeType nnb=cs;
    FloatDP sum=0;
    for(SizeType i=nb-1; i!=0; --i) {
        sum=add(approx,sum,norms[i]);
        if(sum>norms[i-1]) {
            nnb=i;
        }
    }
    nnb=min(nnb,cs);
    // Reduce generators
    Matrix<FloatDP> rG(d,d*nnb);
    project(rG,range(0,d),range(0,d*(nnb-1)))=project(G,range(0,d),range(0,d*(nnb-1)));
    for(SizeType i=0; i!=d; ++i) {
        FloatDP& err=rG[i][d*(nnb-1)+i];
        for(SizeType j=d*(nnb-1); j!=G.column_size(); ++j) {
            err=add(up,err,abs(G[i][j]));
        }
    }
    return Zonotope(z.centre(),rG);
}


Zonotope
orthogonal_approximation(const Zonotope& z)
{

    Vector<FloatDP> c=z.centre();
    Matrix<FloatDP> J=z.generators();
    Vector<FloatDP> e=z.error();

    const SizeType m=J.row_size();
    const SizeType n=J.column_size();

    Matrix<FloatDP> G(m,m+m);

    Matrix< FloatDP > Q;
    Matrix< FloatDP > R;
    make_ltuple(cast_approximate(Q),cast_approximate(R))=orthogonal_decomposition(cast_approximate(J));

    ARIADNE_ASSERT(norm(FloatMatrix(Q*R-J))<1e-8);

    for(SizeType i=0; i!=m;++i) {
        FloatDP a=0;
        for(SizeType j=i; j!=n; ++j) {
            a+=abs(R[i][j]);
        }
        for(SizeType k=0; k!=m; ++k) {
            FloatDP b=Q[k][i]*a;
            G[k][i]=b;
        }
    }

    for(SizeType i=0; i!=m; ++i) { G[i][m+i]=e[i]; }
    return Zonotope(c,G);
// Choose
}


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
    Matrix<FloatDP> G(z.size(),z.number_of_generators());

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
    Matrix<FloatDP> G(z.size(),z.number_of_generators());

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
    Matrix<FloatDP> G(z.size(),z.number_of_generators());

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

Zonotope apply(const ValidatedVectorMultivariateFunction& f, const Zonotope& z) {
    DoublePrecision pr;
    Vector<FloatDPBounds> zc=cast_exact(z.centre());
    Matrix<FloatDPValue> zG=cast_exact(z.generators());
    Vector<FloatDPBounds> ze=cast_exact(z.error())*FloatDPBounds(-1,+1);
    Vector<FloatDPBounds> zb=cast_singleton(z.bounding_box());

    Vector<FloatDPBounds> fc=f(zc);
    Matrix<FloatDPBounds> fJb=jacobian(f,zb);

    Matrix<FloatDPBounds> fJbzG=fJb*zG;

    std::cerr<<"  fJb="<<fJb<<"\n";

    RawFloatVector nzc = cast_raw(midpoint(fc));
    RawFloatMatrix nzG = cast_raw(midpoint(fJbzG));

    Vector<FloatDPBounds> zE(z.number_of_generators(),FloatDPBounds(-1,+1,pr));

    Vector<FloatDPBounds> nzE=(fc-Vector<FloatDPBounds>(nzc)) + (fJbzG-Matrix<FloatDPBounds>(nzG))*zE + fJb*ze;

    Vector<FloatDP> nze(nzE.size()); for(SizeType i=0; i!=nze.size(); ++i) { nze[i]=nzE[i].upper().raw(); }
    std::cerr<<"  nzE="<<nzE<<"\n";
    std::cerr<<"  nze="<<nze<<"\n";

    return Zonotope(nzc,nzG,nze);
}




OutputStream&
operator<<(OutputStream& os, const Zonotope& z)
{
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








ValidatedLowerKleenean
inside(const Zonotope& z, const ExactBoxType& bx)
{
    // FIXME: Avoid casting 'indeterminate'
    return z.bounding_box().inside(bx) || ValidatedKleenean(indeterminate);
}

namespace {

/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[I,z.G], b=z.c, l=[r.l,-o], u=[r.u,+o]
 */
ValidatedKleenean
disjoint(const Zonotope& z, const ExactBoxType& bx)
{
    ARIADNE_ASSERT(z.dimension()==bx.dimension());
    SizeType d=z.dimension();
    SizeType ng=z.number_of_generators();
    Vector<UpperIntervalType> ebx=bx+UpperIntervalType(-1,1)*cast_exact(z.error());
    const Vector<FloatDP>& zc=z.centre();
    const Matrix<FloatDP>& zG=z.generators();
    Matrix<FloatDP> A(d,d+ng);
    Vector<FloatDP> b(d);
    Vector<FloatDP> xl(d+ng);
    Vector<FloatDP> xu(d+ng);

    project(A,range(0,d),range(0,d))=Matrix<FloatDP>::identity(d);
    project(A,range(0,d),range(d,d+ng))=zG;
    b=zc;
    for(SizeType j=0; j!=d; ++j) {
        xl[j]=ebx[j].lower().raw();
        xu[j]=ebx[j].upper().raw();
    }
    for(SizeType j=0; j!=ng; ++j) {
        xl[d+j]=-1;
        xu[d+j]=+1;
    }

    return ! SimplexSolver<FloatDP>().feasible(xl,xu,A,b);
}

}


ValidatedLowerKleenean
separated(const Zonotope& z, const ExactBoxType& bx)
{
    return disjoint(z,bx);
}

ValidatedLowerKleenean
overlaps(const Zonotope& z, const ExactBoxType& bx)
{
    return !disjoint(z,bx);
}

/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[z1.G,z2.G], b=z1.c-z2.c, l=[-o,-o], u=[+o,+o]
 * Still need to take into account errors, particularly in
 * the \a b vector.
 */
ValidatedLowerKleenean
separated(const Zonotope& z1, const Zonotope& z2)
{
    ARIADNE_ASSERT(z1.dimension()==z2.dimension());
    SizeType d=z1.dimension();
    SizeType ng1=z1.number_of_generators();
    SizeType ng2=z2.number_of_generators();
    const Vector<FloatDP>& c1=z1.centre();
    const Matrix<FloatDP>& G1=z1.generators();
    const Vector<FloatDP>& c2=z2.centre();
    const Matrix<FloatDP>& G2=z2.generators();

    Matrix<FloatDP> A(d,ng1+ng2);
    Vector<FloatDP> b(c1-c2);
    Vector<FloatDP> xl(ng1+ng2,FloatDP(-1.0,dp));
    Vector<FloatDP> xu(ng1+ng2,FloatDP(+1.0,dp));

    project(A,range(0,d),range(0,ng1))=G1;
    project(A,range(0,d),range(ng1,ng1+ng2))=G2;

    return ! SimplexSolver<FloatDP>().feasible(xl,xu,A,b);
}


/* Set up LP problem to solve \f$c+Gs=p\f$; \f$-1<=s<=1\f$.
 */
ValidatedKleenean
contains(const Zonotope& z, const ExactPoint& pt)
{
    //std::clog << "Zonotope::contains(const Vector<FloatDP>& )" << std::endl;
    assert(z.dimension()==pt.dimension());
    SizeType m=z.number_of_generators();

    const Matrix<FloatDP>& A=z.generators();
    //FIXME: Incorrect if c-p is not computed exactly, or if e!=0.
    Vector<FloatDP> b=sub(approx,cast_raw(pt),z.centre());
    Vector<FloatDP> xl(m,-1.0);
    Vector<FloatDP> xu(m,1.0);

    ValidatedKleenean result=SimplexSolver<FloatDP>().feasible(xl,xu,A,b);
    return result;
}

struct angle_less {
    Bool operator() (const Vector2d& v1, const Vector2d& v2) const {
        assert(v1.x>0 && v2.x>0);
        return (v1.y/v1.x) < (v2.y/v2.x);
    }
};

Void Zonotope::draw(CanvasInterface& c, const Projection2d& p) const {
    const Zonotope& z=*this;
    SizeType ix=p.x_coordinate(); SizeType iy=p.y_coordinate();

    const Vector<FloatDP>& zc=z.centre();
    const Matrix<FloatDP>& zg=z.generators();
    const Vector<FloatDP>& ze=z.error();

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

} // namespace Ariadne

