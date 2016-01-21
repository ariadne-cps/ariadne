/***************************************************************************
 *            zonotope.cc
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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

#include "function/functional.h"
#include "config.h"

#include "geometry/zonotope.h"

#include <iostream>
#include <vector>
#include <algorithm>

#include "utility/macros.h"
#include "utility/array.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "solvers/linear_programming.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/list_set.h"
#include "function/function.h"
#include "geometry/polytope.h"
#include "output/geometry2d.h"


namespace Ariadne {

inline Vector<Float64> add_approx(const Vector<Float64>& v1, const Vector<Float64>& v2) { return v1+v2; }
inline Vector<Float64> sub_approx(const Vector<Float64>& v1, const Vector<Float64>& v2) { return v1-v2; }

template<class X> class LinearProgram {
  public:
    LinearProgram(const Matrix<X>& A) { ARIADNE_NOT_IMPLEMENTED; }
    Bool is_feasible() { ARIADNE_NOT_IMPLEMENTED; }
};


template<class BS>
ListSet<BS>
subdivide(const BS& bs, const Float64& r)
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



Void
accumulate(Float64& value, Float64& error, Nat n, const Float64* aptr, const Float64* bptr)
{
    ExactIntervalType v=ExactIntervalType(value);
    for(Nat i=0; i!=n; ++i) {
        v+=aptr[i]*bptr[i];
    }
    value=v.midpoint().raw();
    error=add_up(error,v.radius().raw());
}

Vector<Float64>
row_norms(const Matrix<ExactIntervalType>& A)
{
    Nat const& m=A.row_size();
    Nat const& n=A.column_size();
    Vector<Float64> e(m);
    for(Nat i=0; i!=m; ++i) {
        for(Nat j=0; j!=n; ++j) {
            e[i]=add_up(e[i],mag(A[i][j]));
        }
    }
    return e;
}

Vector<Float64>
row_errors(const Matrix<ExactIntervalType>& A)
{
    Nat const& m=A.row_size();
    Nat const& n=A.column_size();
    Vector<Float64> e(m);
    for(Nat i=0; i!=m; ++i) {
        for(Nat j=0; j!=n; ++j) {
            e[i]=add_up(e[i],A[i][j].radius().raw());
        }
    }
    return e;
}

Vector<Float64>
errors(const Vector<ExactIntervalType>& pt)
{
    Vector<Float64> result(pt.size());
    for(Nat i=0; i!=pt.size(); ++i) {
        result[i]=pt[i].radius().raw();
    }
    return result;
}


Vector<Float64>
row_errors(const Vector<ExactIntervalType>& pt, const Matrix<ExactIntervalType>& A)
{
    assert(pt.size()==A.row_size());
    Vector<Float64> result(pt.size());
    for(Nat i=0; i!=A.row_size(); ++i) {
        result[i]=pt[i].radius().raw();
        for(Nat j=0; j!=A.column_size(); ++j) {
            result[i]=add_up(result[i],A[i][j].radius().raw());
        }
    }
    return result;
}

Vector<Float64>
add_up(const Vector<Float64>& v1, const Vector<Float64>& v2)
{
    Vector<Float64> result;
    for(Nat i=0; i!=v1.size(); ++i) {
        result[i]=add_up(v1[i],v2[i]);
    }
    return result;
}

ValidatedKleenean
norm_grtr(const Vector<Float64>& v1, const Vector<Float64>& v2)
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


Zonotope::Zonotope(Nat d)
    : _centre(d), _generators(d,0), _error(d)
{
}


Zonotope::Zonotope(Nat d, Nat m)
    : _centre(d), _generators(d,m), _error(d)
{
}

Zonotope::Zonotope(InitializerList< std::tuple<Float64,InitializerList<Float64>,Float64> > lst)
    : _centre(lst.size()), _generators(lst.size(),lst.size()==0?0u:std::get<1>(*lst.begin()).size()), _error(lst.size())
{
    for(InitializerList< std::tuple<Float64,InitializerList<Float64>,Float64> >::ConstIterator aff_iter=lst.begin();
        aff_iter!=lst.end(); ++aff_iter)
    {
        Nat i=aff_iter-lst.begin();
        this->_centre[i]=std::get<0>(*aff_iter);
        for(Nat j=0; j!=_generators.column_size(); ++j) {
            this->_generators[i][j]=*(std::get<1>(*aff_iter).begin()+j);
        }
        this->_error[i]=std::get<2>(*aff_iter);
        ARIADNE_ASSERT(this->_error[i]>=0);
    }
}


Zonotope::Zonotope(const Vector<Float64>& c, const Matrix<Float64>& G, const Vector<Float64>& e)
    : _centre(c), _generators(G), _error(e)
{
    assert(c.size()==G.row_size());
    assert(c.size()==e.size());
}


Zonotope::Zonotope(const Vector<Float64>& c, const Matrix<Float64>& G)
    : _centre(c), _generators(G), _error(c.size())
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<ExactIntervalType>& c, const Matrix<Float64>& G)
    : _centre(midpoint(c)), _generators(G), _error(errors(c))
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<Float64>& c, const Matrix<ExactIntervalType>& G)
    : _centre(c), _generators(midpoint(G)), _error(row_errors(G))
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<ExactIntervalType>& c, const Matrix<ExactIntervalType>& G)
    : _centre(midpoint(c)), _generators(midpoint(G)), _error(row_errors(c,G))
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

Nat
Zonotope::dimension() const
{
    return this->_centre.size();
}


Nat
Zonotope::number_of_generators() const
{
    return this->_generators.column_size();
}


const Vector<Float64>&
Zonotope::centre() const
{
    return this->_centre;
}


const Matrix<Float64>&
Zonotope::generators() const
{
    return this->_generators;
}


const Vector<Float64>&
Zonotope::error() const
{
    return this->_error;
}


Vector<ExactIntervalType>
Zonotope::domain() const
{
    return Vector<ExactIntervalType>(this->number_of_generators(),ExactIntervalType(-1,1));
}


ExactBoxType
Zonotope::bounding_box() const
{
    const Zonotope& z=*this;
//    std::cerr<<"zD="<<z.domain()<<"\n";
//    std::cerr<<"zc="<<cast_exact(z.centre())<<"\n";
//    std::cerr<<"zG="<<cast_exact(z.generators())<<"\n";
//    std::cerr<<"ze"<<cast_exact(z.error())*ExactIntervalType(-1,1)<<"\n";
//    std::cerr<<"zG*E"<<cast_exact(z.generators())*z.domain()<<"\n";
    ExactBoxType b=cast_exact(z.centre())+(cast_exact(z.generators())*z.domain())+cast_exact(z.error())*ExactIntervalType(-1,1);
//    std::cerr<<"bb="<<b<<"\n";
    return b;
}


Float64
Zonotope::radius() const
{
    return Ariadne::radius(this->bounding_box());
}


ValidatedKleenean
Zonotope::contains(const ExactPoint& pt) const
{
    return Ariadne::contains(*this,pt);
}


ValidatedKleenean
Zonotope::separated(const ExactBoxType& bx) const
{
    return Ariadne::separated(*this,ExactBoxType(bx));
}


ValidatedKleenean
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




Float64
radius(const Zonotope& z)
{
    return Ariadne::radius(z.bounding_box());
}










ExactBoxType
bounding_box(const Zonotope& z)
{
    return z.bounding_box();
}



ListSet< Zonotope >
split(const Zonotope& z)
{
    // FIXME: Not quite guarenteed to give an over-approximation
    typedef ExactIntervalType I;


    ListSet< Zonotope  > result;

    Nat d=z.dimension();
    Nat m=z.number_of_generators();
    Vector<Float64> const& c=z.centre();
    Matrix<Float64> const& G=z.generators();
    Vector<Float64> const& e=z.error();

    Array<Float64> norms(m,0);
    for(Nat j=0; j!=m; ++j) {
        norms[j]=norm(Vector<Float64>(column(G,j)));
    }

    Float64 max_norm=0;
    Nat longest_generator=0;
    for(Nat j=0; j<m; ++j) {
        if(norms[j]>max_norm) {
            max_norm=norms[j];
            longest_generator=j;
        }
    }
    for(Nat k=0; k<d; ++k) {
        if(e[k]>max_norm) {
            max_norm=e[k];
            longest_generator=m+k;
        }
    }

    if(longest_generator<m) {
        Matrix<Float64> new_generators=z.generators();
        Nat j=longest_generator;
        for(Nat i=0; i!=d; ++i) {
            new_generators[i][j]=div_up(new_generators[i][j],2);
        }

        Vector<Float64> v=column(new_generators,j);
        Vector<Float64> new_centre=sub_approx(c,v);
        result.adjoin(Zonotope(new_centre,new_generators,e));
        new_centre=add_approx(c,v);
        result.adjoin(Zonotope(new_centre,new_generators,e));
    } else {
        Nat k=longest_generator-m;
        Vector<Float64> new_centre = z.centre();
        const Matrix<Float64>& new_generators = z.generators();
        Vector<Float64> new_error=e;
        new_error[k]=div_up(new_error[k],2);
        new_centre[k]=add_approx(z.centre()[k],new_error[k]);
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
        new_centre[k]=sub_approx(z.centre()[k],new_error[k]);
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
    }
    return result;
}





Zonotope::Zonotope(const ExactBoxType& r)
    : _centre(r.size()), _generators(r.size(),r.size()), _error(r.size())
{
    Nat d=r.size();
    Vector<Float64>& c=this->_centre;
    Matrix<Float64>& G=this->_generators;
    Vector<Float64>& e=this->_error;
    for(Nat i=0; i!=d; ++i) {
        c[i]=med_approx(r[i].lower().raw(),r[i].upper().raw());
        for(Nat j=0; j!=d; ++j) {
            G[i][j]=0;
        }
        G[i][i]=rad_up(r[i].lower().raw(),r[i].upper().raw());
        e[i]=0;
    }
}







Zonotope
approximation(const Zonotope& z)
{
    return Zonotope(z.centre(),z.generators());
}


Zonotope
over_approximation(const Zonotope& z)
{
    return z;
}


Zonotope
error_free_over_approximation(const Zonotope& z)
{
    Nat d=z.dimension();
    Nat m=z.number_of_generators();

    // Count number of nonzero error values
    Nat e=0;
    for(Nat i=0; i!=d; ++i) {
        if(z.error()[i]!=0) { ++e; }
    }
    Matrix<Float64> nG(d,m+e);
    project(nG,range(0,d),range(0,m))=z.generators();

    Nat j=m;
    for(Nat i=0; i!=d; ++i) {
        if(z.error()[i]!=0) {
            nG[i][j]=z.error()[i];
            ++j;
        }
    }
    return Zonotope(z.centre(),nG);
}





/*
Zonotope
orthogonal_over_approximation(const Zonotope& z)
{
    //assert(iz.size()==iz.number_of_generators());
    typedef ExactIntervalType I;
    Zonotope ez=error_free_over_approximation(z);

    const Vector<Float64>& c=ez.centre();
    const Matrix<Float64>& G=ez.generators();

    Matrix<Float64> aQ,aR;
    make_lpair(aQ,aR)=qr_approx(G);

    Matrix<ExactIntervalType> aQinv=inverse(aQ);
    Matrix<ExactIntervalType> iR=aQinv*G;
    DiagonalMatrix<Float64> aD(::row_norms(iR));

    Matrix<ExactIntervalType> niG=aQ*aD;

    return Zonotope(c,niG);
}
*/


Zonotope
cascade_over_approximation(const Zonotope& z, Nat cs)
{
    using namespace std;

    if(z.number_of_generators()<=z.dimension()*cs) { return z; }

    assert(z.number_of_generators()%z.dimension()==0);

    Nat d=z.dimension();
    Nat nb=z.number_of_generators()/z.dimension(); // number of generator blocks


    const Matrix<Float64>& G=z.generators();
    Array<Float64> norms(nb);
    for(Nat i=0; i!=nb; ++i) {
        norms[i]=Ariadne::norm(Matrix<Float64>(project(G,range(0,d),range(i*d,(i+1)*d))));
    }

    // Compute the new number of blocks
    Nat nnb=cs;
    Float64 sum=0;
    for(Nat i=nb-1; i!=0; --i) {
        sum=add_approx(sum,norms[i]);
        if(sum>norms[i-1]) {
            nnb=i;
        }
    }
    nnb=min(nnb,cs);
    // Reduce generators
    Matrix<Float64> rG(d,d*nnb);
    project(rG,range(0,d),range(0,d*(nnb-1)))=project(G,range(0,d),range(0,d*(nnb-1)));
    for(Nat i=0; i!=d; ++i) {
        Float64& err=rG[i][d*(nnb-1)+i];
        for(Nat j=d*(nnb-1); j!=G.column_size(); ++j) {
            err=add_up(err,abs(G[i][j]));
        }
    }
    return Zonotope(z.centre(),rG);
}



Zonotope
orthogonal_over_approximation(const Zonotope& z)
{

    Zonotope r=error_free_over_approximation(z);
    Matrix<Float64> J=r.generators();

    const Nat m=J.row_size();
    const Nat n=J.column_size();

    Matrix<Float64> R(m,n);

    Array<Float64> column_norm_squares(n);

    for(Nat j=0; j!=n; ++j) { }

    return z;
// Choose
}

Tuple< Matrix<Float64>, Matrix<Float64>, PivotMatrix > orthogonal_decomposition(const Matrix<Float64>& A, Bool allow_pivoting=true) {
    Matrix<Float64Approximation> approximate_matrix=reinterpret_cast<Matrix<Float64Approximation>const&>(A);
    Tuple< Matrix<Float64Approximation>, Matrix<Float64Approximation>, PivotMatrix >
        approximate_decomposition=orthogonal_decomposition(approximate_matrix);
    return reinterpret_cast<Tuple<Matrix<Float64>, Matrix<Float64>, PivotMatrix >const&>(approximate_decomposition);
}


Zonotope
orthogonal_approximation(const Zonotope& z)
{

    Vector<Float64> c=z.centre();
    Matrix<Float64> J=z.generators();
    Vector<Float64> e=z.error();

    const Nat m=J.row_size();
    const Nat n=J.column_size();

    Matrix<Float64> G(m,m+m);

    Matrix< Float64 > Q;
    Matrix< Float64 > R;
    PivotMatrix P;
    make_ltuple(Q,R,P)=orthogonal_decomposition(J,false);

    ARIADNE_ASSERT(norm(FloatMatrix(Q*R-J))<1e-8);

    for(Nat i=0; i!=m;++i) {
        Float64 a=0;
        for(Nat j=i; j!=n; ++j) {
            a+=abs(R[i][j]);
        }
        for(Nat k=0; k!=m; ++k) {
            Float64 b=Q[k][i]*a;
            G[k][i]=b;
        }
    }

    for(Nat i=0; i!=m; ++i) { G[i][m+i]=e[i]; }
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
    Matrix<Float64> G(z.size(),z.number_of_generators());

    Matrix< ExactIntervalType > q=QR.Q();
    Matrix< ExactIntervalType > r=QR.R();
    for(Nat i=0; i!=z.size();++i) {
        ExactIntervalType a=0;
        for(Nat j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(Nat k=0; k!=z.size(); ++k) {
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
    Matrix<Float64> G(z.size(),z.number_of_generators());

    Matrix< ExactIntervalType > q=QR.Q();
    Matrix< ExactIntervalType > r=QR.R();
    for(Nat i=0; i!=z.size();++i) {
        ExactIntervalType a=0;
        for(Nat j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(Nat k=0; k!=z.size(); ++k) {
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
    Matrix<Float64> G(z.size(),z.number_of_generators());

    Matrix< ExactIntervalType > q=QR.Q();
    Matrix< ExactIntervalType > r=QR.R();
    for(Nat i=0; i!=z.size();++i) {
        ExactIntervalType a=0;
        for(Nat j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(Nat k=0; k!=z.size(); ++k) {
            ExactIntervalType b=q(k,i)*a;
            G(k,i)=b.midpoint();
            c[k]+=(b-b.midpoint());
        }
    }
    return Zonotope< ExactIntervalType >(c,G);
}
*/

Zonotope apply(const ValidatedVectorFunction& f, const Zonotope& z) {
    std::cerr<<"Zonotope apply(ValidatedVectorFunction,Zonotope)\n";
    ExactIntervalVectorType zc=z.centre();
    ExactIntervalMatrixType zG=z.generators();
    ExactIntervalVectorType ze=z.error()*ExactIntervalType(-1,+1);
    ExactIntervalVectorType zb=z.bounding_box();

    ExactIntervalVectorType fc=apply(f,zc);
    ExactIntervalMatrixType fJb=jacobian(f,zb);
    ExactIntervalMatrixType fJbzG=fJb*zG;

    std::cerr<<"  fJb="<<fJb<<"\n";

    RawFloatVector nzc = midpoint(fc);
    FloatMatrix nzG = midpoint(fJbzG);

    ExactIntervalVectorType zE(z.number_of_generators(),ExactIntervalType(-1,+1));

    ExactIntervalVectorType nzE=(fc-ExactIntervalVectorType(nzc)) + (fJbzG-ExactIntervalMatrixType(nzG))*zE + fJb*ze;

    RawFloatVector nze(nzE.size()); for(Nat i=0; i!=nze.size(); ++i) { nze[i]=nzE[i].upper(); }
    std::cerr<<"  nzE="<<nzE<<"\n";
    std::cerr<<"  nze="<<nze<<"\n";

    return Zonotope(nzc,nzG,nze);
}





OutputStream&
operator<<(OutputStream& os, const Zonotope& z)
{
    os << "[";
    for(Nat i=0; i!=z.dimension(); ++i) {
            os << (i==0 ? '(' : ',') << z.centre()[i];
    }
    os << "),";
    for(Nat j=0; j!=z.number_of_generators(); ++j) {
        for(Nat i=0; i!=z.dimension(); ++i) {
            os << (i==0 ? '[' : ',') << z.generators()[i][j];
        }
        os << "],";
    }
    for(Nat i=0; i!=z.dimension(); ++i) {
        os << (i==0 ? '{' : ',') << z.error()[i];
    }
    os << '}';
    os << "]";
    return os;
}




InputStream&
operator>>(InputStream& is, Zonotope& z)
{
    Vector<Float64> centre;
    Matrix<Float64> generators;
    char c0,c1,c2;
    is >> c0 >> centre >> c1 >> generators >> c2;
    z = Zonotope(centre,generators);
    return is;
}





ValidatedKleenean
inside(const Zonotope& z, const ExactBoxType& bx)
{
    return z.bounding_box().inside(bx) || indeterminate;
}


ValidatedKleenean
overlaps(const Zonotope& z, const ExactBoxType& bx)
{
    return !separated(z,bx);
}


/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[I,z.G], b=z.c, l=[r.l,-o], u=[r.u,+o]
 */
ValidatedKleenean
separated(const Zonotope& z, const ExactBoxType& bx)
{
    ARIADNE_ASSERT(z.dimension()==bx.dimension());
    SizeType d=z.dimension();
    SizeType ng=z.number_of_generators();
    Vector<ExactIntervalType> ebx=bx+ExactIntervalType(-1,1)*cast_exact(z.error());
    const Vector<Float64>& zc=z.centre();
    const Matrix<Float64>& zG=z.generators();
    Matrix<Float64> A(d,d+ng);
    Vector<Float64> b(d);
    Vector<Float64> xl(d+ng);
    Vector<Float64> xu(d+ng);

    project(A,range(0,d),range(0,d))=Matrix<Float64>::identity(d);
    project(A,range(0,d),range(d,d+ng))=zG;
    b=zc;
    for(SizeType j=0; j!=d; ++j) {
        xl[j]=ebx[j].lower();
        xu[j]=ebx[j].upper();
    }
    for(SizeType j=0; j!=ng; ++j) {
        xl[d+j]=-1;
        xu[d+j]=+1;
    }

    return ! SimplexSolver<Float64>().feasible(xl,xu,A,b);
}


/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[z1.G,z2.G], b=z1.c-z2.c, l=[-o,-o], u=[+o,+o]
 * Still need to take into account errors, particularly in
 * the \a b vector.
 */
ValidatedKleenean
separated(const Zonotope& z1, const Zonotope& z2)
{
    ARIADNE_ASSERT(z1.dimension()==z2.dimension());
    SizeType d=z1.dimension();
    SizeType ng1=z1.number_of_generators();
    SizeType ng2=z2.number_of_generators();
    const Vector<Float64>& c1=z1.centre();
    const Matrix<Float64>& G1=z1.generators();
    const Vector<Float64>& c2=z2.centre();
    const Matrix<Float64>& G2=z2.generators();

    Matrix<Float64> A(d,ng1+ng2);
    Vector<Float64> b(c1-c2);
    Vector<Float64> xl(ng1+ng2,-1.0);
    Vector<Float64> xu(ng1+ng2,+1.0);

    project(A,range(0,d),range(0,ng1))=G1;
    project(A,range(0,d),range(ng1,ng1+ng2))=G2;

    return ! SimplexSolver<Float64>().feasible(xl,xu,A,b);
}


/* Set up LP problem to solve \f$c+Ge=p\f$; \f$-1<=e<=1\f$.
 */
ValidatedKleenean
contains(const Zonotope& z, const ExactPoint& pt)
{
    //std::clog << "Zonotope::contains(const Vector<Float64>& )" << std::endl;
    assert(z.dimension()==pt.dimension());
    Nat m=z.number_of_generators();

    const Matrix<Float64>& A=z.generators();
    Vector<Float64> b=pt-z.centre();
    Vector<Float64> xl(m,-1.0);
    Vector<Float64> xu(m,1.0);

    ValidatedKleenean result=SimplexSolver<Float64>().feasible(xl,xu,A,b);
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
    Nat ix=p.x_coordinate(); Nat iy=p.y_coordinate();

    const Vector<Float64>& zc=z.centre();
    const Matrix<Float64>& zg=z.generators();
    const Vector<Float64>& ze=z.error();

    double eps=1.0/(1ul<<31);
    Point2d pc(zc[ix],zc[iy]);
    std::vector< Vector2d > pg;
    for(Nat j=0; j!=z.number_of_generators(); ++j) {
        Vector2d g(zg[ix][j],zg[iy][j]);
        if(g.x<0) { g=-g; }
        else if (g.x==0) { g.x=eps; }
        pg.push_back(g);
    }
    if(ze[ix]>0) { pg.push_back(Vector2d(ze[ix],0.0)); }
    if(ze[iy]>0) { pg.push_back(Vector2d(eps,ze[iy])); }

    std::sort(pg.begin(),pg.end(),angle_less());

    const Nat npg=pg.size();
    Point2d pt=pc;
    for(Nat i=0; i!=npg; ++i) {
        pt-=pg[i];
    }

    c.move_to(pt.x,pt.y);
    for(Nat i=0; i!=npg; ++i) {
        pt+=2*pg[i];
        c.line_to(pt.x,pt.y);
    }
    for(Nat i=0; i!=npg; ++i) {
        pt-=2*pg[i];
        c.line_to(pt.x,pt.y);
    }
    c.fill();
}



} // namespace Ariadne

