/***************************************************************************
 *            linear_programming.cc
 *
 *  Copyright 2008-10  Pieter Collins
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

#include "config.h"

#include "tuple.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "affine.h"
#include "linear_programming.h"

#include "macros.h"
#include "logging.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>

static const int verbosity=0;

namespace Ariadne {


template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& x) {
    return prod(A,x); }
template<class X> inline Vector<X> operator*(const Vector<X>& y, const Matrix<X>& A) {
    return prod(y,A); }
template<class X> inline Matrix<X> operator*(const Matrix<X>& A, const Matrix<X>& B) {
    return prod(A,B); }


// Compute A A^T
template<class X> inline Matrix<X> sprod(const Matrix<X>& A) {
    Matrix<X> S(A.row_size(),A.row_size());
    // Compute lower triangle
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=i+1u; ++j) {
            S[i][j]=0.0;
            for(uint k=0; k!=A.column_size(); ++k) {
                S[i][j]+=A[i][k]*A[j][k];
            }
        }
    }
    for(uint j=1; j!=A.row_size(); ++j) {
        for(uint i=0; i!=j; ++i) {
            S[i][j]=S[j][i];
        }
    }
    return S;
}




// Compute the product A:=XA, where X is diagonal
void ldmul(Matrix<Float>& A, const Vector<Float>& x) {
    ARIADNE_ASSERT(A.row_size()==x.size());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.row_size(); ++j) {
            A[i][j]*=x[i];
        }
    }
}

// Compute the product A:=AX, where X is diagonal
void rdmul(Matrix<Float>& A, const Vector<Float>& x) {
    ARIADNE_ASSERT(A.column_size()==x.size());
    for(uint j=0; j!=A.row_size(); ++j) {
        for(uint i=0; i!=A.row_size(); ++i) {
            A[i][j]*=x[j];
        }
    }
}

// Set r[i]=x[i]-c for i=1,...,n
Vector<Float> esub(const Vector<Float>& x, const Float& c) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]-c;
    }
    return r;
}

// Set r[i]=x[i]*y[i] for i=1,...,n
void _emul(Vector<Float>& r, const Vector<Float>& x, const Vector<Float>& y) {
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i];
    }
}

// Set r[i]=x[i]*y[i] for i=1,...,n
Vector<Float> emul(const Vector<Float>& x, const Vector<Float>& y) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i];
    }
    return r;
}

// Set r[i]=x[i]*y[i]+z for i=1,...,n
Vector<Float> efma(const Vector<Float>& x, const Vector<Float>& y, const Float& z) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i]+z;
    }
    return r;
}

// Return r[i]=x[i]/z[i] for i=1,...,n
void _ediv(Vector<Float>& r, const Vector<Float>& x, const Vector<Float>& z) {
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]/z[i];
    }
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<Float> ediv(const Vector<Float>& x, const Vector<Float>& z) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]/z[i];
    }
    return r;
}

// Return r[i]=1/y[i] for i=1,...,n
Vector<Float> erec(const Vector<Float>& v) {
    Vector<Float> r(v.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=1.0/v[i];
    }
    return r;
}

// Compute R=ADA^T for diagonal D
Matrix<Float> mdtmul(const Matrix<Float>& A, const Vector<Float>& D)
{
    Matrix<Float> R(A.row_size(),A.row_size());
    ARIADNE_ASSERT(D.size()==A.column_size());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint k=0; k!=A.column_size(); ++k) {
            Float ADik=A[i][k]*D[k];
            for(uint j=0; j!=A.row_size(); ++j) {
                R[i][j]+=ADik*A[j][k];
            }
        }
    }
    return R;
}

bool all_greater(const Vector<Float>& x, const Float& e) {
    for(uint i=0; i!=x.size(); ++i) {
        if(x[i]<=e) { return false; }
    }
    return true;
}

bool all_less(const Vector<Float>& x, const Float& e) {
    for(uint i=0; i!=x.size(); ++i) {
        if(x[i]>=e) { return false; }
    }
    return true;
}

bool all_greater(const Vector<Float>& x1, const Vector<Float>& x2) {
    for(uint i=0; i!=x1.size(); ++i) {
        if(x1[i]<=x2[i]) { return false; }
    }
    return true;
}



Interval dot(const Vector<Float>& c, const Vector<Interval>& x) {
    Interval r=0.0;
    for(uint i=0; i!=x.size(); ++i) {
        r+=c[i]*x[i];
    }
    return r;
}

void _set_es_matrix(Matrix<Float>& S, const Matrix<Float>& A, const Vector<Float>& D) {
    // D=diag(Db,Dc,Dl,Du)
    // S=[ A*(Dc+Db)*AT + (Du+Dl),     A*(Dc-Db)*e + (Du-Dl)*e ]
    //   [ eT*(Dc-Db)*AT + eT*(Du-Dl), eT*(Dc+Db+Du+Dl)*e      ]
    static const uint m=A.row_size();
    static const uint n=A.column_size();

    ARIADNE_ASSERT(S.row_size()==m+1);

    const Float* Db=&D[0];
    const Float* Dc=Db+n;
    const Float* Dl=Dc+n;
    const Float* Du=Dl+m;
    Float* Sptr=&S[0][0];

    //Clear matrix
    for(uint i=0; i!=(m+1u)*(m+1u); ++i) {
        Sptr[i]=0.0;
    }

    // Set upper-triangular values
    for(uint i=0; i!=m; ++i) {
        for(uint k=0; k!=n; ++k) {
            Float ADik=A[i][k]*(Db[k]+Dc[k]);
            for(uint j=i; j!=m; ++j) {
                S[i][j]+=ADik*A[j][k];
            }
            S[i][m]+=A[i][k]*(Dc[k]-Db[k]);
        }
        S[i][i]+=(Du[i]+Dl[i]);
        S[i][m]+=(Du[i]-Dl[i]);
    }
    for(uint h=0; h!=(2*(m+n)); ++h) {
        S[m][m]+=D[h];
    }

    // Set lower part
    for(uint i=1; i!=m+1; ++i) {
        for(uint j=0; j!=i; ++j) {
            S[i][j]=S[j][i];
        }
    }
}


void _set_s_matrix(Matrix<Float>& S, const Matrix<Float>& A, const Vector<Float>& Dc, const Vector<Float>& Dd) {
    // S=A*Dc*AT+Dd

    static const uint m=A.row_size();
    static const uint n=A.column_size();

    ARIADNE_ASSERT(S.row_size()==m);
    ARIADNE_ASSERT(S.column_size()==m);
    ARIADNE_ASSERT(Dc.size()==n);
    ARIADNE_ASSERT(Dd.size()==m);

    //Clear matrix
    Float* Sptr=&S[0][0];
    for(uint i=0; i!=m*m; ++i) {
        Sptr[i]=0.0;
    }

    // Set upper-triangular values
    for(uint i=0; i!=m; ++i) {
        S[i][i]=Dd[i];
        for(uint k=0; k!=n; ++k) {
            Float ADik=A[i][k]*Dc[k];
            for(uint j=i; j!=m; ++j) {
                S[i][j]+=ADik*A[j][k];
            }
        }
    }
    // Set lower part
    for(uint i=1; i!=m; ++i) {
        for(uint j=0; j!=i; ++j) {
            S[i][j]=S[j][i];
        }
    }
}



Float compute_mu(const Vector<Float>& xl, const Vector<Float>& xu,
                 const Vector<Float>& x, const Vector<Float>& zl, const Vector<Float>& zu)
{
    const uint n=x.size();
    Float mu = 0.0;
    for(uint i=0; i!=n; ++i) {
        if(xl[i]!=-infty) { mu += ((x[i]-xl[i])*zl[i]); }
        if(xu[i]!=+infty) { mu += ((xu[i]-x[i])*zu[i]); }
    }
    mu /= (2.0*n);
    return mu;
}

Float compute_mu(const Vector<Float>& x, const Vector<Float>& z)
{
    const uint n=x.size();
    Float mu = 0.0;
    for(uint i=0; i!=n; ++i) {
        mu += x[i]*z[i];
    }
    mu /= n;
    return mu;
}


Interval
InteriorPointSolver::validate(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c, const Vector<Float>& x, const Vector<Float>& y) const
{
    // Find a primal-dual feasible point (x,y) suitable as an input
    // to an optimization solver. The slack variables are A
    ARIADNE_ASSERT(A.row_size()==b.size());
    ARIADNE_ASSERT(A.column_size()==c.size());

    Matrix<Interval> IA=A;
    Matrix<Interval> IS=prod(transpose(IA),IA);
    Matrix<Interval> ISinv=inverse(IS);
    Vector<Interval> dx=IA*Vector<Interval>(x); dx=IS*dx; dx=dx*IA;
    Vector<Interval> Ix=x-dx;
    Vector<Interval> Iy(y);
    Vector<Interval> Iz=prod(y,IA)-c;
    for(uint i=0; i!=Ix.size(); ++i) {
        ARIADNE_ASSERT(Ix[i].lower()>=0.0);
    }
    for(uint i=0; i!=Iz.size(); ++i) {
        ARIADNE_ASSERT(Iz[i].lower()>=0.0);
    }
    return hull(dot(c,Ix),dot(b,Iy));
}


tribool
InteriorPointSolver::validate_feasibility(const Matrix<Float>& A, const Vector<Float>& b,
                                          const Vector<Float>& x, const Vector<Float>& y) const
{
    const uint m=A.row_size();
    const uint n=A.column_size();

    ARIADNE_LOG(2,"InteriorPointSolver::validate_feasibility(A,b,x,y)\n");
    ARIADNE_LOG(3,"A="<<A<<" b="<<b<<"\n");
    ARIADNE_LOG(3,"x="<<x<<" y="<<y<<"\n");
    ARIADNE_LOG(4,"Ax-b="<<Vector<Float>(A*x-b)<<" yA="<<y*A<<" yb="<<dot(y,b)<<"\n");

    // If yA < 0 and yb > 0, then problem is infeasible
    tribool result = false;

    set_rounding_upward();
    for(uint j=0; j!=x.size(); ++j) {
        Float yAj=0.0;
        for(uint i=0; i!=m; ++i) {
            yAj += y[i]*A[i][j];
        }
        if(yAj>=0.0) { result=indeterminate; break; }
    }
    set_rounding_downward();
    if(result==false) {
        Float yb=0.0;
        for(uint i=0; i!=y.size(); ++i) {
            yb += b[i]*y[i];
        }
        if(yb <= 0.0) { result=indeterminate; }
    }
    set_rounding_to_nearest();
    if(definitely(!result)) { return result; }

    // x should be an approximate solution to Ax=b
    // Use the fact that for any x, x'=(x + A^T (AA^T)^{-1}(b-Ax)) satisfies Ax'=0
    Vector<Interval> ivlx = x;
    Vector<Interval> ivle = b-A*ivlx;

    ARIADNE_LOG(4,"[e] = "<<ivle << "\n");
    Matrix<Interval> ivlS(m,m);
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=i1; i2!=m; ++i2) {
            for(uint j=0; j!=n; ++j) {
                ivlS[i1][i2] += mul_ivl(A[i1][j],A[i2][j]);
            }
        }
    }
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=0; i2!=i1; ++i2) {
            ivlS[i1][i2]=ivlS[i2][i1];
        }
    }

    Vector<Interval> ivld = solve(ivlS,ivle) * A;
    ARIADNE_LOG(4,"[d] = "<<ivld << "\n");

    ivlx = x + ivld;

    ARIADNE_LOG(3,"[x] = "<<ivlx << " A[x]-b="<<Vector<Interval>(A*ivlx-b)<<"\n");
    result=true;
    for(uint i=0; i!=n; ++i) {
        if(ivlx[i].lower()<=0.0) {
            result = indeterminate; break;
        }
    }

    return result;

}


tribool
InteriorPointSolver::validate_constrained_feasibility(const Matrix<Float>& A, const Vector<Float>& b,
                                                      const Vector<Float>& xl, const Vector<Float>& xu,
                                                      const Vector<Float>& x, const Vector<Float>& y) const
{
    const uint m=A.row_size();
    const uint n=A.column_size();

    // x should be an approximate solution to Ax=b
    // Use the fact that for any x, x'=(x + A^T (AA^T)^{-1}(b-Ax)) satisfies Ax'=0
    Vector<Interval> ivlx = x;
    Vector<Interval> ivle = b-A*ivlx;

    Matrix<Interval> ivlS(m,m);
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=i1; i2!=m; ++i2) {
            for(uint j=0; j!=n; ++j) {
                ivlS[i1][i2] += mul_ivl(A[i1][j],A[i2][j]);
            }
        }
    }
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=0; i2!=i1; ++i2) {
            ivlS[i1][i2]=ivlS[i2][i1];
        }
    }

    Vector<Interval> ivld = solve(ivlS,ivle) * A;

    ivlx += ivld;

    ARIADNE_LOG(2,"[x] = "<<ivlx);
    tribool result=true;
    for(uint i=0; i!=n; ++i) {
        if(ivlx[i].lower()<=xl[i] || ivlx[i].upper()>=xu[i]) {
            result = indeterminate; break;
        }
    }

    if(result==true) { return result; }

    // If yb - max(yA,0) xu + min(yA,0) xl > 0, then problem is infeasible
    // Evaluate lower bound for yb - max(z,0) xu + min(z,0) xl
    Vector<Interval> z=Vector<Interval>(y)*A;
    Float mx = 0.0;
    set_rounding_downward();
    for(uint i=0; i!=y.size(); ++i) {
        mx += b[i]*y[i];
    }
    for(uint i=0; i!=x.size(); ++i) {
        Float neg_xui = -xu[i];
        if(z[i].upper()>0.0) { mx += z[i].upper() * neg_xui; }
        if(z[i].lower()<0.0) { mx += z[i].lower() * xl[i]; }
    }
    set_rounding_to_nearest();
    if(mx>0.0) { return false; }

    return indeterminate;
}


tuple< Vector<Float>, Vector<Float>, Vector<Float> >
InteriorPointSolver::optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c) const
{
    Vector<Float> x,y,z;
    return this->hotstarted_optimize(A,b,c,x,y,z);
}


tuple< Vector<Float>, Vector<Float>, Vector<Float> >
InteriorPointSolver::hotstarted_optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                         Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    ARIADNE_ASSERT(A.row_size()==b.size());
    ARIADNE_ASSERT(A.column_size()==c.size());
    ARIADNE_ASSERT(A.column_size()==x.size());
    ARIADNE_ASSERT(A.row_size()==y.size());
    ARIADNE_ASSERT(A.column_size()==z.size());

    const double maxerror=1e-3;
    const uint maxsteps=10;

    Float cx,yb;
    cx=dot(c,x);
    yb=dot(y,b);
    ARIADNE_ASSERT(yb<=cx);

    ARIADNE_LOG(2,"A="<<A<<" b="<<b<<" c="<<c<<"\n");
    ARIADNE_LOG(2,"x="<<x<<" y="<<y<<" z="<<z<<"\n");

    uint steps=0;
    while(steps<maxsteps && (cx-yb)>maxerror) {
        this->_optimization_step(A,b,c,x,y,z);
        ++steps;
    }

    return make_tuple(x,y,z);

}


tuple<Float, Vector<Float>, Vector<Float> >
InteriorPointSolver::constrained_optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                          const Vector<Float>& xl, const Vector<Float>& xu) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::constrained_optimize(A,b,c,xl,xu)\n");
    ARIADNE_LOG(2,"A="<<A<<", b="<<b<<", c="<<c<<"\n");
    ARIADNE_LOG(2,"xl="<<xl<<", xu="<<xu<<"\n");

    const uint m = b.size();
    const uint n = c.size();
    Vector<Float> y(m, 0.0);
    Vector<Float> x(n);
    Vector<Float> zl(n);
    Vector<Float> zu(n);
    for(uint i=0; i!=n; ++i) {
        if(xl[i]==-infty) {
            if(xu[i]==+infty) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+infty) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-infty) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+infty) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }

    bool done;
    do {
        done = this->_constrained_optimization_step(A,b,c,xl,xu, x,y,zl,zu);
    } while(!done);

    do {
        this->_constrained_optimization_step(A,b,c,xl,xu, x,y,zl,zu);
    } while(dot(c,x)-dot(y,b)>1e-4);

    // Todo: check for optimality
    return make_tuple(dot(c,x),x,y);
}




// Consider feasibility of  bl<=yB<=bu; yC<=c; dl<=y<=du
// Augment system by considering b<=yA-ve; yA+ve<=c; max(t)
tribool
InteriorPointSolver::primal_feasible(const Matrix<Float>& A, const Vector<Float>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::primal_feasible(A,b)\n");
    ARIADNE_LOG(3,"A="<<A<<" b="<<b<<"\n");
    Vector<Float> x(A.column_size(),1.0);
    Vector<Float> y(A.row_size(),0.0);
    Vector<Float> z(A.column_size(),1.0);

    const double THRESHOLD = 1e-6;

    while(true) {
        tribool result=this->_primal_feasibility_step(A,b,x,y,z);
        if(!indeterminate(result)) { return result; }
        if(compute_mu(x,z)<THRESHOLD) { return indeterminate; }
    }

}



tribool
InteriorPointSolver::constrained_feasible(const Matrix<Float>& A, const Vector<Float>& b,
                                          const Vector<Float>& xl, const Vector<Float>& xu) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::constrained_feasible(A,b,xl,xu)\n");
    ARIADNE_LOG(3,"A="<<A<<", b="<<b<<"\n");
    ARIADNE_LOG(3,"xl="<<xl<<", xu="<<xu<<"\n");

    const uint m = A.row_size();
    const uint n = A.column_size();
    Vector<Float> c(n, 0.0);
    Vector<Float> y(m, 0.0);
    Vector<Float> x(n);
    Vector<Float> zl(n);
    Vector<Float> zu(n);
    for(uint i=0; i!=n; ++i) {
        if(xl[i]==-infty) {
            if(xu[i]==+infty) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+infty) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-infty) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+infty) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }
    Vector<Interval> X=hull(xl,xu);

    const double THRESHOLD = 1e-8;
    uint step=0;
    while(step++<24) {
        LinearProgramStatus result=this->_constrained_feasibility_step(A,b,xl,xu, x,y,zl,zu);
        if(result==PRIMAL_DUAL_FEASIBLE || result==PRIMAL_FEASIBLE) {
            tribool validated_feasible=this->validate_constrained_feasibility(A,b,xl,xu, x,y);
            if(definitely(validated_feasible)) { return true; }
        }
        Interval yb=dot(IntervalVector(y),IntervalVector(b));
        Interval yAX = dot( (IntervalVector(y)*A), X );
        if(disjoint(yb,yAX)) { return false; }
        if(result==DEGENERATE_FEASIBILITY) { ARIADNE_LOG(2,"  degenerate\n"); return indeterminate; }
        if(compute_mu(xl,xu, x,zl,zu)<THRESHOLD ) { ARIADNE_LOG(2,"  threshold\n"); return indeterminate; }
    }
    return indeterminate;
}

tribool
InteriorPointSolver::constrained_dual_feasible(const Matrix<Float>& A, const Vector<Float>& cl, const Vector<Float>& cu,
                                               const Vector<Float>& l, const Vector<Float>& u) const
{
    ARIADNE_NOT_IMPLEMENTED;
}









LinearProgramStatus
InteriorPointSolver::_optimization_step(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                        Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    ARIADNE_LOG(4,"InteriorPointSolver::_optimization_step(A,b,c, x,y,z)\n");
    ARIADNE_LOG(5,"x="<<x<<" y="<<y<<" z="<<z<<"\n");
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> dx(n),dy(n),dz(n);
    Vector<Float> nx(n),ny(n),nz(n),nxz(n);
    Vector<Float> rx(m),ry(n),rz(n),r(n);
    Matrix<Float> S(m,m),Sinv(m,m);

    Float mu=sigma*dot(x,z)/n; // duality constant
    ARIADNE_LOG(5,"mu="<<mu<<"\n");

    // rx = Ax-b; ry=yA+z-c; rz=x.z-mu
    rx=A*x-b;
    ry=y*A+z-c;
    rz=esub(emul(x,z),mu);
    ARIADNE_LOG(5,"rx="<<rx<<" ry="<<ry<<" rz="<<rz<<"\n");

    r=rx-A*(ediv(rz-emul(x,ry),z));

    S=mdtmul(A,ediv(x,z));
    Sinv=inverse(S);
    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // Adx = rx; ATdy + dz = ry; Zdx + X dz = rz
    // dx = (rz - X dz)/Z
    // dz = ry - AT dy
    // A(X/Z)AT dy = rx - A Zinv (rz - X ry)
    dy=Sinv*r;
    dz=ry-dy*A;
    dx=ediv(rz-emul(x,dz),z);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<" dz="<<dz<<"\n");

    // Try to enforce feasibility or dual feasibility
    Float alphax=1.0;
    nx=x-dx;
    while ( !all_greater(emul(nx,z),gamma*mu) ) {
        alphax=alphax*scale;
        nx=(x-alphax*dx);
    }

    Float alphaz=1.0;
    nz=z-dz;
    while ( !all_greater(emul(x,nz),gamma*mu) ) {
        alphaz=alphaz*scale;
        nz=(z-alphaz*dz);
    }
    ny=(y-alphaz*dy);
    ARIADNE_LOG(6,"alphax="<<alphax<<" nx="<<nx<<" alphaz="<<alphaz<<" ny="<<ny<<" nz="<<nz<<"\n");

    x=nx; y=ny; z=nz;
    ARIADNE_LOG(5,"cx="<<dot(c,x)<<" yb="<<dot(y,b)<<" Ax-b="<<Vector<Float>(prod(A,x)-b)<<" yA+z-c="<<Vector<Float>(y*A+z-c)<<"\n");

    if(alphax==1.0) { return alphaz==1.0 ? PRIMAL_DUAL_FEASIBLE : PRIMAL_FEASIBLE; }
    else { return alphaz==1.0 ? DUAL_FEASIBLE : INDETERMINATE_FEASIBILITY; }
}


LinearProgramStatus
InteriorPointSolver::_constrained_optimization_step(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                                    const Vector<Float>& xl, const Vector<Float>& xu,
                                                    Vector<Float>& x, Vector<Float>& y, Vector<Float>& zl, Vector<Float>& zu) const
{
    ARIADNE_LOG(4,"InteriorPointSolver::_constrained_optimization_step(A,b,c,xl,xu, x,y,zl,zu)\n");
    ARIADNE_LOG(4,"x="<<x<<", y="<<y<<", zl="<<zl<<", zu="<<zu<<"\n");

    static const double gamma=1.0/256;
    static const double sigma=1.0/8;
    static const double scale=0.75;


    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> dx(n),dy(m),dzl(n),dzu(n);
    Vector<Float> nx,ny,nzl,nzu;
    Vector<Float> rx(m),ry(n),rzl(n),rzu(n);
    Matrix<Float> S(m,m),Sinv(m,m);
    DiagonalMatrix<Float> Xl(xl), Xu(xu), X(x), Zl(zl), Zu(zu);

    Float mu = compute_mu(xl,xu, x,zl,zu) * sigma;
    mu=1.0;
    ARIADNE_LOG(4,"mu="<<mu<<"\n");

    // rx = Ax-b; ry=yA+zl-zu-c; rzl=(x-xl).zl-mu; rzu=(xu-x).zu-mu.
    rx=A*x-b;
    ry=y*A+(zl-zu)-c;
    for(uint i=0; i!=n; ++i) {
        if(xl[i]!=-infty) { rzl[i] = (x[i]-xl[i])*zl[i] - mu; } else { rzl[i]=0.0; }
        if(xu[i]!=+infty) { rzu[i] = (xu[i]-x[i])*zu[i] - mu; } else { rzu[i]=0.0; }
    }
    ARIADNE_LOG(5,"rx="<<rx<<", ry="<<ry<<", rzl="<<rzl<<", rzu="<<rzu<<"\n");

    // A dx = rx;  AT dy + dzl - dzu = ry;  Zl dx + (X-Xl) dzl = rzl; -Zu dx + (Xu-X) dzu = rzu;

    // Derivation
    //   AT dy + (rzl-Zl*dx)/(X-Xl) - (rzu+Zu*dx)/(Xu-X) = ry
    //   AT dy - ( Zl/(X-Xl) + Zu/(Xu-X) ) dx = ry - rzl/(X-Xl) + rzu/(Xu-X)
    //   AT dy - Dinv dx = ryz
    //   D AT dy - dx = D ryz
    //      dx = D (AT dy - ryz)
    //   A D AT dy - A dx = A D ryz
    //   S dy - rx = A D ryz
    //    S dy = A D ryz + rx

    //   ryz = ( ry - rzl/(X-Xl) + rzu/(Xu-X) )
    //   D = 1/( Zu/(Xu-X) + Zl/(X-Xl) )
    DiagonalMatrix<Float> D(erec(ediv(zu,xu-x)+ediv(zl,x-xl)));
    Vector<Float> ryz = ry - ediv(rzl,Vector<Float>(x-xl)) + ediv(rzu,Vector<Float>(xu-x));
    S=mdtmul(A,D.diagonal());;
    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // dzl = (rzl - Zl dx) / (X-Xl)
    // dzu = (rzu + Zu dx) / (Xu-X)
    // dx = D (AT dy - ryz)
    // (A D AT) dy = rx + A D ryz

    dy = solve(S, Vector<Float>( rx + A * (D * ryz) ) );
    dx = D * Vector<Float>(dy*A - ryz);
    dzl = Vector<Float>(rzl-Zl*dx)/(X-Xl);
    dzu = Vector<Float>(rzu+Zu*dx)/(Xu-X);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<" dzl="<<dzl<<" dzu="<<dzu<<"\n");

    ARIADNE_LOG(7,"A*dx="<<(A*dx)<<" AT*dy+dzl-dzu="<<(dy*A+dzl-dzu)<<" Zl*dx+(X-Xl)*dzl="<<(Zl*dx+(X-Xl)*dzl)<<" -Zu*dx+(Xu-X)*dzu="<<(-(Zu*dx)+(Xu-X)*dzu)<<"\n");
    ARIADNE_LOG(7,"A*dx-rx="<<(A*dx-rx)<<" AT*dy+dzl-dzu-ry="<<(dy*A+dzl-dzu-ry)<<" Zl*dx+(X-Xl)*dzl-rzl="<<(Zl*dx+(X-Xl)*dzl-rzl)<<" -Zu*dx+(Xu-X)*dzu-rzu="<<(-(Zu*dx)+(Xu-X)*dzu-rzu)<<"\n");
    // Try to enforce feasibility or dual feasibility
    Float alphax=1.0;
    nx=x-dx;
    while ( !all_greater(emul(nx-xl,zl),gamma*mu) || !all_greater(emul(xu-nx,zu),gamma*mu) ) {
        alphax=alphax*scale;
        nx=(x-alphax*dx);
        if(alphax<gamma*mu/4096) { return DEGENERATE_FEASIBILITY; }
    }

    Float alphaz=1.0;
    nzl=zl-dzl; nzu=zu-dzu;
    while ( !all_greater(emul(nx-xl,nzl),gamma*mu) || !all_greater(emul(xu-nx,nzu),gamma*mu) ) {
        alphaz=alphaz*scale;
        nzl=(zl-alphaz*dzl);
        nzu=(zu-alphaz*dzu);
        if(alphaz<gamma*mu/4096) { return DEGENERATE_FEASIBILITY; }
        //ARIADNE_LOG(9,"alphaz="<<alphaz<<" nzl="<<nzl<<" nzu="<<nzu<<"\n");
    }
    ny=(y-alphaz*dy);
    ARIADNE_LOG(5,"alphax="<<alphax<<" nx="<<nx<<" alphaz="<<alphaz<<" ny="<<ny<<" nzl="<<nzl<<" nzu="<<nzu<<"\n");

    x=nx; y=ny; zl=nzl; zu=nzu;
    ARIADNE_LOG(5,"cx="<<dot(c,x)<<" yb="<<dot(y,b)<<" Ax-b="<<Vector<Float>(prod(A,x)-b)<<" yA+(zl-zu)-c="<<Vector<Float>(y*A+(zl-zu)-c)<<"\n");

    if(alphax==1.0 && alphaz==1.0) { return PRIMAL_DUAL_FEASIBLE; }
    if(alphax==1.0) { return PRIMAL_FEASIBLE; }
    if(alphaz==1.0) { return DUAL_FEASIBLE; }
    return INDETERMINATE_FEASIBILITY;
}



LinearProgramStatus
InteriorPointSolver::_primal_feasibility_step(const Matrix<Float>& A, const Vector<Float>& b,
                                              Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    Vector<Float> c(x.size(),0.0);
    return this->_optimization_step(A,b,c, x,y,z);

    ARIADNE_LOG(4,"primal_feasibility_step(A,b,x,y)\n");
    ARIADNE_LOG(5,"x="<<x<<" y="<<y<<" z="<<z<<"\n");
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> dx(n),dy(n),dz(n);
    Vector<Float> nx(n),ny(n),nz(n);
    Vector<Float> rx(m),ry(n),rz(n);
    Matrix<Float> S(m,m),Sinv(m,m);

    DiagonalMatrix<Float> Z(z);
    DiagonalMatrix<Float> X(x);

    Float mu = dot(x,z)/n * sigma;
    ARIADNE_LOG(5,"mu="<<mu<<"\n");

    rx = y*A+z;
    ry = A*x-b;
    rz = efma(x,z,-mu);
    ARIADNE_LOG(5,"rx=yA+z="<<rx<<" ry=Ax-b="<<ry<<" rz=x.z-mu="<<rz<<"\n");

    S=mdtmul(A,ediv(x,z));
    Sinv=inverse(S);

    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // S dy = ry - A Zinv (rz - X rx)
    // dz = rx - AT dy
    // dx = Zinv * (rz - X dz)
    dy = Sinv * Vector<Float>( ry - A * Vector<Float>(rz-X*rx)/Z );
    dz = rx - dy * A;
    dx = Vector<Float>(rz-X*dz)/Z;
    ARIADNE_LOG(4,"  dx="<<dx<<" dy="<<dy<<" dz="<<dz<<"\n");

    // Check for feasibility
    nx=x-dx;
    if(all_greater(nx,0.0)) {
        ARIADNE_LOG(4,"nx="<<nx<<" Ax-b="<<Vector<Float>(A*nx-b)<<"\n");
        x=nx; return PRIMAL_FEASIBLE;
    }

    // Check for infeasibility
    nz=z-dz; ny=y-dy;
    if(all_greater(nz,0.0) && dot(y,b)>0.0) {
        ARIADNE_LOG(4,"ny="<<ny<<" nz="<<nz<<" yA="<<Vector<Float>(ny*A)<<" yb="<<dot(ny,b)<<" yA+z="<<Vector<Float>(ny*A+nz)<<"\n");
        x=nx; y=ny; z=nz; return DUAL_FEASIBLE;
    }

    Float alpha=1.0;
    while ( !all_greater(nx,0.0) || !all_greater(emul(nx,nz),gamma*mu) ) {
        alpha *= scale;
        nx=(x-alpha*dx);
        nz=(z-alpha*dz);
        ARIADNE_LOG(6,"    alpha="<<alpha<<" nx="<<nx<<" nz="<<nz<<"\n");
    }
    ny=(y-alpha*dy);

    x=nx; y=ny; z=nz;

    ARIADNE_LOG(4,"nx="<<x<<" ny="<<y<<" nz="<<z<<"\n");
    ARIADNE_LOG(5,"yb="<<dot(y,b)<<" yA="<<y*A<<" Ax-b="<<Vector<Float>(A*x-b)<<" yA+z="<<Vector<Float>(y*A+z)<<" x.z="<<emul(x,z)<<"\n");

    return INDETERMINATE_FEASIBILITY;
}



LinearProgramStatus
InteriorPointSolver::_constrained_feasibility_step(const Matrix<Float>& A, const Vector<Float>& b,
                                                   const Vector<Float>& xl, const Vector<Float>& xu,
                                                   Vector<Float>& x, Vector<Float>& y, Vector<Float>& zl, Vector<Float>& zu) const
{
    Vector<Float> c(x.size(),0.0);
    return this->_constrained_optimization_step(A,b,c,xl,xu,x,y,zl,zu);


    ARIADNE_LOG(4,"InteriorPointSolver::constrained_feasibility_step(A,b,xl,xu, x,y,zl,zu)\n");
    ARIADNE_LOG(5,"x="<<x<<" y="<<y<<" zl="<<zl<<" zu="<<zu<<"\n");
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> wl(n),wu(n);
    Vector<Float> dx(n),dy(n),dz(n),dzl(n),dzu(n);
    Vector<Float> nx(n),ny(n),nzl(n),nzu(n);
    Vector<Float> rx(m),ry(n),rzl(n),rzu(n);
    Matrix<Float> S(m,m),Sinv(m,m);

    DiagonalMatrix<Float> Xl(xl), Xu(xu), X(x), Zl(zl), Zu(zu);

    // Compute central path parameter mu
    Float mu = 0.0;
    for(uint i=0; i!=n; ++i) {
        if(xl[i]!=-infty) { mu += (x[i]-xl[i])*zl[i]; }
        if(xu[i]!=+infty) { mu += (xu[i]-x[i])*zu[i]; }
    }
    mu *= ( sigma / (2*n) );
    ARIADNE_LOG(5,"mu="<<mu<<"\n");

    rx = y*A+(zl-zu);
    ry = A*x-b;
    rzl = efma(x-xl,zl,-mu);
    rzu = efma(xu-x,zu,-mu);
    //Vector<Float> ryz = ry - ediv(rzl,x-xl) + ediv(rzu,xu-x);
    ARIADNE_LOG(5,"rx=yA+z="<<rx<<" ry=Ax-b="<<ry<<" rzl=(x-xl).zl-mu="<<rzl<<" rzu=(xu-x).zu-mu="<<rzu<<"\n");

    DiagonalMatrix<Float> D(erec(ediv(wl,zl)+ediv(wu,zu)));
    S=mdtmul(A,D.diagonal());
    Sinv=inverse(S);

    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // Ax-b=rx=0; yA+zl-zu=ry=0; (x-xl).zl-mu=rzl=0; (xu-x).zu-mu=rzu=0
    //
    // dzl = (rzl - Zl dx) / (X-Xl)
    // dzu = (rzu + Zu dx) / (Xu-X)
    // (Zu/(Xu-X) + Zl/(X-Xl)) dx = AT dy - ( ry - rzl/(X-Xl) + rzu/(Xu-X) )
    // (A D AT) dy = rx - A D ( ry - rzl/(X-Xl) + rzu/(Xu-X) )

    //dy = Sinv * Vector<Float>( rx - A * (D * ryz) );
    //dx = D * (dy*A - ryz);
    //dzl = Vector<Float>(rzl-Zl*dx)/(X-Xl);
    //dzu = Vector<Float>(rzu-Zu*dx)/(Xu-X);
    ARIADNE_LOG(4,"  dx="<<dx<<" dy="<<dy<<" dzl="<<dzl<<" dzu="<<dzu<<"\n");

    // Check for feasibility
    nx=x-dx;
    if(all_greater(nx,0.0)) {
        ARIADNE_LOG(4,"nx="<<nx<<" Ax-b="<<Vector<Float>(A*nx-b)<<"\n");
        x=nx;
        return PRIMAL_FEASIBLE;
    }

    // Check for infeasibility
    nzl=zl-dzl; nzu=zu-dzu; ny=y-dy;
    if(all_greater(nzl,0.0) && all_greater(nzu,0.0) && dot(y,b)>0.0) {
        ARIADNE_LOG(4,"ny="<<ny<<" nzl="<<nzl<<" nzu="<<nzu<<" yA="<<Vector<Float>(ny*A)<<" yb="<<dot(ny,b)<<" yA+zl-zu="<<Vector<Float>(ny*A+nzl-nzu)<<"\n");
        x=nx; y=ny; zl=nzl; zu=nzu;
        return DUAL_FEASIBLE;
    }

    Float alpha=1.0;
    while ( !all_greater(nx,0.0) || !all_greater(emul(nx-xl,nzl),gamma*mu) || !all_greater(emul(xu-nx,nzu),gamma*mu) ) {
        alpha *= scale;
        nx=(x-alpha*dx);
        nzl=(zl-alpha*dzl);
        nzu=(zu-alpha*dzu);
        ARIADNE_LOG(6,"    alpha="<<alpha<<" nx="<<nx<<" nzl="<<nzl<<" nzu="<<nzu<<"\n");
    }
    ny=(y-alpha*dy);

    x=nx; y=ny; zl=nzl; zu=nzu;

    ARIADNE_LOG(4,"nx="<<x<<" ny="<<y<<" nzl="<<zl<<" nzu="<<zu<<"\n");
    ARIADNE_LOG(5,"yb="<<dot(y,b)<<" yA="<<y*A<<" Ax-b="<<Vector<Float>(A*x-b)<<" yA+zl-zu="<<Vector<Float>(y*A+zl-zu)<<
                  " (x-xl).zl="<<emul(x-xl,zl)<<" (xu-x).zu"<<emul(xu-x,zu)<<"\n");

    return INDETERMINATE_FEASIBILITY;
}



} // namespace Ariadne

