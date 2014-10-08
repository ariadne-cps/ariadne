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

inline Vector<Float> add_approx(const Vector<Float>& v1, const Vector<Float>& v2) { return v1+v2; }
inline Vector<Float> sub_approx(const Vector<Float>& v1, const Vector<Float>& v2) { return v1-v2; }

template<class X> class LinearProgram {
  public:
    LinearProgram(const Matrix<X>& A) { ARIADNE_NOT_IMPLEMENTED; }
    bool is_feasible() { ARIADNE_NOT_IMPLEMENTED; }
};


template<class BS>
ListSet<BS>
subdivide(const BS& bs, const Float& r)
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



void
accumulate(Float& value, Float& error, uint n, const Float* aptr, const Float* bptr)
{
    ExactInterval v=ExactInterval(value);
    for(uint i=0; i!=n; ++i) {
        v+=aptr[i]*bptr[i];
    }
    value=v.midpoint().raw();
    error=add_up(error,v.radius().raw());
}

Vector<Float>
row_norms(const Matrix<ExactInterval>& A)
{
    uint const& m=A.row_size();
    uint const& n=A.column_size();
    Vector<Float> e(m);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            e[i]=add_up(e[i],mag(A[i][j]));
        }
    }
    return e;
}

Vector<Float>
row_errors(const Matrix<ExactInterval>& A)
{
    uint const& m=A.row_size();
    uint const& n=A.column_size();
    Vector<Float> e(m);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            e[i]=add_up(e[i],A[i][j].radius().raw());
        }
    }
    return e;
}

Vector<Float>
errors(const Vector<ExactInterval>& pt)
{
    Vector<Float> result(pt.size());
    for(uint i=0; i!=pt.size(); ++i) {
        result[i]=pt[i].radius().raw();
    }
    return result;
}


Vector<Float>
row_errors(const Vector<ExactInterval>& pt, const Matrix<ExactInterval>& A)
{
    assert(pt.size()==A.row_size());
    Vector<Float> result(pt.size());
    for(uint i=0; i!=A.row_size(); ++i) {
        result[i]=pt[i].radius().raw();
        for(uint j=0; j!=A.column_size(); ++j) {
            result[i]=add_up(result[i],A[i][j].radius().raw());
        }
    }
    return result;
}

Vector<Float>
add_up(const Vector<Float>& v1, const Vector<Float>& v2)
{
    Vector<Float> result;
    for(uint i=0; i!=v1.size(); ++i) {
        result[i]=add_up(v1[i],v2[i]);
    }
    return result;
}

tribool
norm_grtr(const Vector<Float>& v1, const Vector<Float>& v2)
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


Zonotope::Zonotope(uint d)
    : _centre(d), _generators(d,0), _error(d)
{
}


Zonotope::Zonotope(uint d, uint m)
    : _centre(d), _generators(d,m), _error(d)
{
}

Zonotope::Zonotope(std::initializer_list< std::tuple<Float,std::initializer_list<Float>,Float> > lst)
    : _centre(lst.size()), _generators(lst.size(),lst.size()==0?0u:std::get<1>(*lst.begin()).size()), _error(lst.size())
{
    for(std::initializer_list< std::tuple<Float,std::initializer_list<Float>,Float> >::const_iterator aff_iter=lst.begin();
        aff_iter!=lst.end(); ++aff_iter)
    {
        uint i=aff_iter-lst.begin();
        this->_centre[i]=std::get<0>(*aff_iter);
        for(uint j=0; j!=_generators.column_size(); ++j) {
            this->_generators[i][j]=*(std::get<1>(*aff_iter).begin()+j);
        }
        this->_error[i]=std::get<2>(*aff_iter);
        ARIADNE_ASSERT(this->_error[i]>=0);
    }
}


Zonotope::Zonotope(const Vector<Float>& c, const Matrix<Float>& G, const Vector<Float>& e)
    : _centre(c), _generators(G), _error(e)
{
    assert(c.size()==G.row_size());
    assert(c.size()==e.size());
}


Zonotope::Zonotope(const Vector<Float>& c, const Matrix<Float>& G)
    : _centre(c), _generators(G), _error(c.size())
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<ExactInterval>& c, const Matrix<Float>& G)
    : _centre(midpoint(c)), _generators(G), _error(errors(c))
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<Float>& c, const Matrix<ExactInterval>& G)
    : _centre(c), _generators(midpoint(G)), _error(row_errors(G))
{
    assert(c.size()==G.row_size());
}


Zonotope::Zonotope(const Vector<ExactInterval>& c, const Matrix<ExactInterval>& G)
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

bool
operator==(const Zonotope& z1, const Zonotope& z2)
{
    return (z1._centre==z2._centre) && (z1._generators==z2._generators)
        && (z1._error==z2._error);
}

uint
Zonotope::dimension() const
{
    return this->_centre.size();
}


uint
Zonotope::number_of_generators() const
{
    return this->_generators.column_size();
}


const Vector<Float>&
Zonotope::centre() const
{
    return this->_centre;
}


const Matrix<Float>&
Zonotope::generators() const
{
    return this->_generators;
}


const Vector<Float>&
Zonotope::error() const
{
    return this->_error;
}


Vector<ExactInterval>
Zonotope::domain() const
{
    return Vector<ExactInterval>(this->number_of_generators(),ExactInterval(-1,1));
}


ExactBox
Zonotope::bounding_box() const
{
    const Zonotope& z=*this;
//    std::cerr<<"zD="<<z.domain()<<"\n";
//    std::cerr<<"zc="<<make_exact(z.centre())<<"\n";
//    std::cerr<<"zG="<<make_exact(z.generators())<<"\n";
//    std::cerr<<"ze"<<make_exact(z.error())*ExactInterval(-1,1)<<"\n";
//    std::cerr<<"zG*E"<<make_exact(z.generators())*z.domain()<<"\n";
    ExactBox b=make_exact(z.centre())+(make_exact(z.generators())*z.domain())+make_exact(z.error())*ExactInterval(-1,1);
//    std::cerr<<"bb="<<b<<"\n";
    return b;
}


Float
Zonotope::radius() const
{
    return Ariadne::radius(this->bounding_box());
}


tribool
Zonotope::contains(const ExactPoint& pt) const
{
    return Ariadne::contains(*this,pt);
}


tribool
Zonotope::separated(const ExactBox& bx) const
{
    return Ariadne::separated(*this,ExactBox(bx));
}


tribool
Zonotope::inside(const ExactBox& bx) const
{
    return Ariadne::inside(*this,ExactBox(bx));
}



std::ostream&
Zonotope::write(std::ostream& os) const
{
    return os << *this;
}





tribool
empty(const Zonotope& z)
{
    return false;
}


tribool
bounded(const Zonotope& z)
{
    return true;
}




Float
radius(const Zonotope& z)
{
    return Ariadne::radius(z.bounding_box());
}










ExactBox
bounding_box(const Zonotope& z)
{
    return z.bounding_box();
}



ListSet< Zonotope >
split(const Zonotope& z)
{
    // FIXME: Not quite guarenteed to give an over-approximation
    typedef ExactInterval I;


    ListSet< Zonotope  > result;

    uint d=z.dimension();
    uint m=z.number_of_generators();
    Vector<Float> const& c=z.centre();
    Matrix<Float> const& G=z.generators();
    Vector<Float> const& e=z.error();

    Array<Float> norms(m,0);
    for(uint j=0; j!=m; ++j) {
        norms[j]=norm(Vector<Float>(column(G,j)));
    }

    Float max_norm=0;
    uint longest_generator=0;
    for(uint j=0; j<m; ++j) {
        if(norms[j]>max_norm) {
            max_norm=norms[j];
            longest_generator=j;
        }
    }
    for(uint k=0; k<d; ++k) {
        if(e[k]>max_norm) {
            max_norm=e[k];
            longest_generator=m+k;
        }
    }

    if(longest_generator<m) {
        Matrix<Float> new_generators=z.generators();
        uint j=longest_generator;
        for(uint i=0; i!=d; ++i) {
            new_generators[i][j]=div_up(new_generators[i][j],2);
        }

        Vector<Float> v=column(new_generators,j);
        Vector<Float> new_centre=sub_approx(c,v);
        result.adjoin(Zonotope(new_centre,new_generators,e));
        new_centre=add_approx(c,v);
        result.adjoin(Zonotope(new_centre,new_generators,e));
    } else {
        uint k=longest_generator-m;
        Vector<Float> new_centre = z.centre();
        const Matrix<Float>& new_generators = z.generators();
        Vector<Float> new_error=e;
        new_error[k]=div_up(new_error[k],2);
        new_centre[k]=add_approx(z.centre()[k],new_error[k]);
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
        new_centre[k]=sub_approx(z.centre()[k],new_error[k]);
        result.adjoin(Zonotope(new_centre,new_generators,new_error));
    }
    return result;
}





Zonotope::Zonotope(const ExactBox& r)
    : _centre(r.size()), _generators(r.size(),r.size()), _error(r.size())
{
    uint d=r.size();
    Vector<Float>& c=this->_centre;
    Matrix<Float>& G=this->_generators;
    Vector<Float>& e=this->_error;
    for(uint i=0; i!=d; ++i) {
        c[i]=med_approx(r[i].lower().raw(),r[i].upper().raw());
        for(uint j=0; j!=d; ++j) {
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
    uint d=z.dimension();
    uint m=z.number_of_generators();

    // Count number of nonzero error values
    uint e=0;
    for(uint i=0; i!=d; ++i) {
        if(z.error()[i]!=0) { ++e; }
    }
    Matrix<Float> nG(d,m+e);
    project(nG,range(0,d),range(0,m))=z.generators();

    uint j=m;
    for(uint i=0; i!=d; ++i) {
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
    typedef ExactInterval I;
    Zonotope ez=error_free_over_approximation(z);

    const Vector<Float>& c=ez.centre();
    const Matrix<Float>& G=ez.generators();

    Matrix<Float> aQ,aR;
    make_lpair(aQ,aR)=qr_approx(G);

    Matrix<ExactInterval> aQinv=inverse(aQ);
    Matrix<ExactInterval> iR=aQinv*G;
    DiagonalMatrix<Float> aD(::row_norms(iR));

    Matrix<ExactInterval> niG=aQ*aD;

    return Zonotope(c,niG);
}
*/


Zonotope
cascade_over_approximation(const Zonotope& z, uint cs)
{
    using namespace std;

    if(z.number_of_generators()<=z.dimension()*cs) { return z; }

    assert(z.number_of_generators()%z.dimension()==0);

    uint d=z.dimension();
    uint nb=z.number_of_generators()/z.dimension(); // number of generator blocks


    const Matrix<Float>& G=z.generators();
    Array<Float> norms(nb);
    for(uint i=0; i!=nb; ++i) {
        norms[i]=Ariadne::norm(Matrix<Float>(project(G,range(0,d),range(i*d,(i+1)*d))));
    }

    // Compute the new number of blocks
    uint nnb=cs;
    Float sum=0;
    for(uint i=nb-1; i!=0; --i) {
        sum=add_approx(sum,norms[i]);
        if(sum>norms[i-1]) {
            nnb=i;
        }
    }
    nnb=min(nnb,cs);
    // Reduce generators
    Matrix<Float> rG(d,d*nnb);
    project(rG,range(0,d),range(0,d*(nnb-1)))=project(G,range(0,d),range(0,d*(nnb-1)));
    for(uint i=0; i!=d; ++i) {
        Float& err=rG[i][d*(nnb-1)+i];
        for(uint j=d*(nnb-1); j!=G.column_size(); ++j) {
            err=add_up(err,abs(G[i][j]));
        }
    }
    return Zonotope(z.centre(),rG);
}



Zonotope
orthogonal_over_approximation(const Zonotope& z)
{

    Zonotope r=error_free_over_approximation(z);
    Matrix<Float> J=r.generators();

    const uint m=J.row_size();
    const uint n=J.column_size();

    Matrix<Float> R(m,n);

    Array<Float> column_norm_squares(n);

    for(uint j=0; j!=n; ++j) { }

    return z;
// Choose
}

Tuple< Matrix<Float>, Matrix<Float>, PivotMatrix > orthogonal_decomposition(const Matrix<Float>& A, bool allow_pivoting=true) {
    Matrix<ApproximateFloat> approximate_matrix=reinterpret_cast<Matrix<ApproximateFloat>const&>(A);
    Tuple< Matrix<ApproximateFloat>, Matrix<ApproximateFloat>, PivotMatrix >
        approximate_decomposition=orthogonal_decomposition(approximate_matrix);
    return reinterpret_cast<Tuple<Matrix<Float>, Matrix<Float>, PivotMatrix >const&>(approximate_decomposition);
}


Zonotope
orthogonal_approximation(const Zonotope& z)
{

    Vector<Float> c=z.centre();
    Matrix<Float> J=z.generators();
    Vector<Float> e=z.error();

    const uint m=J.row_size();
    const uint n=J.column_size();

    Matrix<Float> G(m,m+m);

    Matrix< Float > Q;
    Matrix< Float > R;
    PivotMatrix P;
    make_ltuple(Q,R,P)=orthogonal_decomposition(J,false);

    ARIADNE_ASSERT(norm(FloatMatrix(Q*R-J))<1e-8);

    for(uint i=0; i!=m;++i) {
        Float a=0;
        for(uint j=i; j!=n; ++j) {
            a+=abs(R[i][j]);
        }
        for(uint k=0; k!=m; ++k) {
            Float b=Q[k][i]*a;
            G[k][i]=b;
        }
    }

    for(uint i=0; i!=m; ++i) { G[i][m+i]=e[i]; }
    return Zonotope(c,G);
// Choose
}


/*

Zonotope<ExactInterval,R>
orthogonal_over_approximation(const Zonotope<R,R>& z)
{
    // FIXME: Subdivide in zero order as well!
    static bool warn=true;
    if(warn) {
        ARIADNE_WARN("orthogonal_over_approximation(Zonotope<I,R>) does not over-approximate roundoff errors.\n);
        warn=false;
    }
    Zonotope<R,R> oaz=over_approximation(z);

    QRMatrix< ExactInterval > QR(oaz.generators());
    ExactPoint< ExactInterval > c(oaz.centre());
    Matrix<Float> G(z.size(),z.number_of_generators());

    Matrix< ExactInterval > q=QR.Q();
    Matrix< ExactInterval > r=QR.R();
    for(uint i=0; i!=z.size();++i) {
        ExactInterval a=0;
        for(uint j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(uint k=0; k!=z.size(); ++k) {
            ExactInterval b=q(k,i)*a;
            G(k,i)=b.midpoint();
            c[k]+=(b-b.midpoint());
        }
    }
    return Zonotope<R,R>(midpoint(c),G);
}


Zonotope<ExactInterval,R>
orthogonal_over_approximation(const Zonotope<ExactInterval,R>& z)
{
    Zonotope<R,R> oaz=over_approximation(z);

    QRMatrix< ExactInterval > QR(oaz.generators());
    ExactPoint< ExactInterval > c(oaz.centre());
    Matrix<Float> G(z.size(),z.number_of_generators());

    Matrix< ExactInterval > q=QR.Q();
    Matrix< ExactInterval > r=QR.R();
    for(uint i=0; i!=z.size();++i) {
        ExactInterval a=0;
        for(uint j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(uint k=0; k!=z.size(); ++k) {
            ExactInterval b=q(k,i)*a;
            G(k,i)=b.midpoint();
            c[k]+=(b-b.midpoint());
        }
    }
    return Zonotope<ExactInterval,R>(c,G);
}


Zonotope< ExactInterval >
orthogonal_over_approximation(const Zonotope< ExactInterval >& z)
{
    Zonotope<R,R> oaz=over_approximation(z);

    QRMatrix< ExactInterval > QR(oaz.generators());
    ExactPoint< ExactInterval > c(oaz.centre());
    Matrix<Float> G(z.size(),z.number_of_generators());

    Matrix< ExactInterval > q=QR.Q();
    Matrix< ExactInterval > r=QR.R();
    for(uint i=0; i!=z.size();++i) {
        ExactInterval a=0;
        for(uint j=i; j!=z.number_of_generators(); ++j) {
            a+=r[i][j];
        }
        for(uint k=0; k!=z.size(); ++k) {
            ExactInterval b=q(k,i)*a;
            G(k,i)=b.midpoint();
            c[k]+=(b-b.midpoint());
        }
    }
    return Zonotope< ExactInterval >(c,G);
}
*/

Zonotope apply(const ValidatedVectorFunction& f, const Zonotope& z) {
    std::cerr<<"Zonotope apply(ValidatedVectorFunction,Zonotope)\n";
    ExactIntervalVector zc=z.centre();
    ExactIntervalMatrix zG=z.generators();
    ExactIntervalVector ze=z.error()*ExactInterval(-1,+1);
    ExactIntervalVector zb=z.bounding_box();

    ExactIntervalVector fc=apply(f,zc);
    ExactIntervalMatrix fJb=jacobian(f,zb);
    ExactIntervalMatrix fJbzG=fJb*zG;

    std::cerr<<"  fJb="<<fJb<<"\n";

    RawFloatVector nzc = midpoint(fc);
    FloatMatrix nzG = midpoint(fJbzG);

    ExactIntervalVector zE(z.number_of_generators(),ExactInterval(-1,+1));

    ExactIntervalVector nzE=(fc-ExactIntervalVector(nzc)) + (fJbzG-ExactIntervalMatrix(nzG))*zE + fJb*ze;

    RawFloatVector nze(nzE.size()); for(uint i=0; i!=nze.size(); ++i) { nze[i]=nzE[i].upper(); }
    std::cerr<<"  nzE="<<nzE<<"\n";
    std::cerr<<"  nze="<<nze<<"\n";

    return Zonotope(nzc,nzG,nze);
}





std::ostream&
operator<<(std::ostream& os, const Zonotope& z)
{
    os << "[";
    for(uint i=0; i!=z.dimension(); ++i) {
            os << (i==0 ? '(' : ',') << z.centre()[i];
    }
    os << "),";
    for(uint j=0; j!=z.number_of_generators(); ++j) {
        for(uint i=0; i!=z.dimension(); ++i) {
            os << (i==0 ? '[' : ',') << z.generators()[i][j];
        }
        os << "],";
    }
    for(uint i=0; i!=z.dimension(); ++i) {
        os << (i==0 ? '{' : ',') << z.error()[i];
    }
    os << '}';
    os << "]";
    return os;
}




std::istream&
operator>>(std::istream& is, Zonotope& z)
{
    Vector<Float> centre;
    Matrix<Float> generators;
    char c0,c1,c2;
    is >> c0 >> centre >> c1 >> generators >> c2;
    z = Zonotope(centre,generators);
    return is;
}





tribool
inside(const Zonotope& z, const ExactBox& bx)
{
    return z.bounding_box().inside(bx) || indeterminate;
}


tribool
overlaps(const Zonotope& z, const ExactBox& bx)
{
    return !separated(z,bx);
}


/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[I,z.G], b=z.c, l=[r.l,-o], u=[r.u,+o]
 */
tribool
separated(const Zonotope& z, const ExactBox& bx)
{
    ARIADNE_ASSERT(z.dimension()==bx.dimension());
    size_t d=z.dimension();
    size_t ng=z.number_of_generators();
    Vector<ExactInterval> ebx=bx+ExactInterval(-1,1)*make_exact(z.error());
    const Vector<Float>& zc=z.centre();
    const Matrix<Float>& zG=z.generators();
    Matrix<Float> A(d,d+ng);
    Vector<Float> b(d);
    Vector<Float> xl(d+ng);
    Vector<Float> xu(d+ng);

    project(A,range(0,d),range(0,d))=Matrix<Float>::identity(d);
    project(A,range(0,d),range(d,d+ng))=zG;
    b=zc;
    for(size_t j=0; j!=d; ++j) {
        xl[j]=ebx[j].lower();
        xu[j]=ebx[j].upper();
    }
    for(size_t j=0; j!=ng; ++j) {
        xl[d+j]=-1;
        xu[d+j]=+1;
    }

    return ! SimplexSolver<Float>().feasible(xl,xu,A,b);
}


/* Set up constrained linear program Ax=b, l\leq x\leq u.
 * Here, A=[z1.G,z2.G], b=z1.c-z2.c, l=[-o,-o], u=[+o,+o]
 * Still need to take into account errors, particularly in
 * the \a b vector.
 */
tribool
separated(const Zonotope& z1, const Zonotope& z2)
{
    ARIADNE_ASSERT(z1.dimension()==z2.dimension());
    size_t d=z1.dimension();
    size_t ng1=z1.number_of_generators();
    size_t ng2=z2.number_of_generators();
    const Vector<Float>& c1=z1.centre();
    const Matrix<Float>& G1=z1.generators();
    const Vector<Float>& c2=z2.centre();
    const Matrix<Float>& G2=z2.generators();

    Matrix<Float> A(d,ng1+ng2);
    Vector<Float> b(c1-c2);
    Vector<Float> xl(ng1+ng2,-1.0);
    Vector<Float> xu(ng1+ng2,+1.0);

    project(A,range(0,d),range(0,ng1))=G1;
    project(A,range(0,d),range(ng1,ng1+ng2))=G2;

    return ! SimplexSolver<Float>().feasible(xl,xu,A,b);
}


/* Set up LP problem to solve \f$c+Ge=p\f$; \f$-1<=e<=1\f$.
 */
tribool
contains(const Zonotope& z, const ExactPoint& pt)
{
    //std::clog << "Zonotope::contains(const Vector<Float>& )" << std::endl;
    assert(z.dimension()==pt.dimension());
    uint m=z.number_of_generators();

    const Matrix<Float>& A=z.generators();
    Vector<Float> b=pt-z.centre();
    Vector<Float> xl(m,-1.0);
    Vector<Float> xu(m,1.0);

    tribool result=SimplexSolver<Float>().feasible(xl,xu,A,b);
    return result;
}


struct angle_less {
    bool operator() (const Vector2d& v1, const Vector2d& v2) const {
        assert(v1.x>0 && v2.x>0);
        return (v1.y/v1.x) < (v2.y/v2.x);
    }
};

void Zonotope::draw(CanvasInterface& c, const Projection2d& p) const {
    const Zonotope& z=*this;
    uint ix=p.x_coordinate(); uint iy=p.y_coordinate();

    const Vector<Float>& zc=z.centre();
    const Matrix<Float>& zg=z.generators();
    const Vector<Float>& ze=z.error();

    double eps=1.0/(1ul<<31);
    Point2d pc(zc[ix],zc[iy]);
    std::vector< Vector2d > pg;
    for(uint j=0; j!=z.number_of_generators(); ++j) {
        Vector2d g(zg[ix][j],zg[iy][j]);
        if(g.x<0) { g=-g; }
        else if (g.x==0) { g.x=eps; }
        pg.push_back(g);
    }
    if(ze[ix]>0) { pg.push_back(Vector2d(ze[ix],0.0)); }
    if(ze[iy]>0) { pg.push_back(Vector2d(eps,ze[iy])); }

    std::sort(pg.begin(),pg.end(),angle_less());

    const uint npg=pg.size();
    Point2d pt=pc;
    for(uint i=0; i!=npg; ++i) {
        pt-=pg[i];
    }

    c.move_to(pt.x,pt.y);
    for(uint i=0; i!=npg; ++i) {
        pt+=2*pg[i];
        c.line_to(pt.x,pt.y);
    }
    for(uint i=0; i!=npg; ++i) {
        pt-=2*pg[i];
        c.line_to(pt.x,pt.y);
    }
    c.fill();
}



} // namespace Ariadne

