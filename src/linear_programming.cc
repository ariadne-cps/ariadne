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

#include "functional.h"
#include "config.h"

#include "tuple.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "affine.h"
#include "linear_programming.h"

#include "macros.h"
#include "logging.h"

static const int verbosity=0;

namespace Ariadne {



// Return r[i]=x[i]-c for i=1,...,n
Vector<Float> esub(const Vector<Float>& x, const Float& c) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]-c;
    }
    return r;
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<Float> emul(const Vector<Float>& x, const Vector<Float>& y) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i];
    }
    return r;
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<Float> ediv(const Vector<Float>& x, const Vector<Float>& z) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]/z[i];
    }
    return r;
}

// Return r[i]=x[i]*y[i]+z for i=1,...,n
Vector<Float> efma(const Vector<Float>& x, const Vector<Float>& y, const Float& z) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i]+z;
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
Matrix<Float> adat(const Matrix<Float>& A, const Vector<Float>& D)
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


Float compute_mu(const Vector<Float>& xl, const Vector<Float>& xu,
                 const Vector<Float>& x, const Vector<Float>& zl, const Vector<Float>& zu)
{
    const uint n=x.size();
    Float mu = 0.0;
    for(uint i=0; i!=n; ++i) {
        if(xl[i]!=-inf) { mu += ((x[i]-xl[i])*zl[i]); }
        if(xu[i]!=+inf) { mu += ((xu[i]-x[i])*zu[i]); }
    }
    mu /= (2.0*n);
    return mu;
}





tribool InteriorPointSolver::
validate_feasibility(const Vector<Float>& xl, const Vector<Float>& xu,
                     const Matrix<Float>& A, const Vector<Float>& b,
                     const Vector<Float>& x, const Vector<Float>& y) const
{
    const uint m=A.row_size();
    const uint n=A.column_size();

    // x should be an approximate solution to Ax=b
    // Use the fact that for any x, x'=(x + A^T (AA^T)^{-1}(b-Ax)) satisfies Ax'=0
    Vector<Interval> ivlx = x;
    Vector<Interval> ivle = make_exact(b)-make_exact(A)*ivlx;

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

    Vector<Interval> ivld = solve(ivlS,ivle) * make_exact(A);

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
    Vector<Interval> z=make_exact(y)*make_exact(A);
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


Tuple< Float, Vector<Float>, Vector<Float> >
InteriorPointSolver::minimise(const Vector<Float>& c,
                              const Vector<Float>& xl, const Vector<Float>& xu,
                              const Matrix<Float>& A, const Vector<Float>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::minimise(c,xl,xu,A,b)\n");
    ARIADNE_LOG(2,"A="<<A<<", b="<<b<<", c="<<c<<"\n");
    ARIADNE_LOG(2,"xl="<<xl<<", xu="<<xu<<"\n");

    const uint m = b.size();
    const uint n = c.size();
    Vector<Float> y(m, 0.0);
    Vector<Float> x(n);
    Vector<Float> zl(n);
    Vector<Float> zu(n);
    for(uint i=0; i!=n; ++i) {
        if(xl[i]==-inf) {
            if(xu[i]==+inf) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+inf) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-inf) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+inf) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }

    bool done;
    do {
        done = this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
    } while(!done);

    do {
        this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
    } while(dot(c,x)-dot(y,b)>1e-4);

    // Todo: check for optimality
    return make_tuple(dot(c,x),x,y);
}



Tuple< Float, Vector<Float>, Vector<Float> >
InteriorPointSolver::
hotstarted_minimise(const Vector<Float>& c,
                    const Vector<Float>& xl, const Vector<Float>& xu,
                    const Matrix<Float>& A, const Vector<Float>& b,
                    Vector<Float>& x, Vector<Float>& y, Vector<Float>& zl, Vector<Float>& zu) const
{
    ARIADNE_ASSERT(A.column_size()==c.size());
    ARIADNE_ASSERT(A.column_size()==xl.size());
    ARIADNE_ASSERT(A.column_size()==xu.size());
    ARIADNE_ASSERT(A.row_size()==b.size());
    ARIADNE_ASSERT(A.column_size()==x.size());
    ARIADNE_ASSERT(A.row_size()==y.size());
    ARIADNE_ASSERT(A.column_size()==zl.size());
    ARIADNE_ASSERT(A.column_size()==zu.size());

    const double maxerror=1e-3;
    const uint maxsteps=10;

    Float cx,yb;
    cx=dot(c,x);
    yb=dot(y,b);
    ARIADNE_ASSERT(yb<=cx);

    ARIADNE_LOG(2,"xl="<<xl<<" xu="<<xu<<" A="<<A<<" b="<<b<<" c="<<c<<"\n");
    ARIADNE_LOG(2,"x="<<x<<" y="<<y<<" z="<<zl<<" zu="<<zu<<"\n");

    uint steps=0;
    while(steps<maxsteps && (cx-yb)>maxerror) {
        this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
        ++steps;
    }

    return make_tuple(dot(c,x),x,y);

}


tribool
InteriorPointSolver::
feasible(const Vector<Float>& xl, const Vector<Float>& xu,
         const Matrix<Float>& A, const Vector<Float>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::feasible(xl,xu,A,b)\n");
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
        if(xl[i]==-inf) {
            if(xu[i]==+inf) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+inf) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<=xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-inf) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+inf) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }
    Vector<Interval> X=hull(xl,xu);

    const double THRESHOLD = 1e-8;
    uint step=0;
    while(step++<24) {
        LinearProgramStatus result=this->_feasibility_step(xl,xu,A,b, x,y,zl,zu);
        if(result==PRIMAL_DUAL_FEASIBLE || result==PRIMAL_FEASIBLE) {
            tribool validated_feasible=this->validate_feasibility(xl,xu,A,b, x,y);
            if(definitely(validated_feasible)) { return true; }
        }
        Interval yb=dot(IntervalVector(y),IntervalVector(b));
        Interval yAX = dot( make_exact(y)*make_exact(A), X );
        if(disjoint(yb,yAX)) { return false; }
        if(result==DEGENERATE_FEASIBILITY) { ARIADNE_LOG(2,"  degenerate\n"); return indeterminate; }
        if(compute_mu(xl,xu, x,zl,zu)<THRESHOLD ) { ARIADNE_LOG(2,"  threshold\n"); return indeterminate; }
    }
    return indeterminate;
}












LinearProgramStatus InteriorPointSolver::
_minimisation_step(const Vector<Float>& c, const Vector<Float>& xl, const Vector<Float>& xu, const Matrix<Float>& A, const Vector<Float>& b,
                   Vector<Float>& x, Vector<Float>& y, Vector<Float>& zl, Vector<Float>& zu) const
{
    ARIADNE_LOG(4,"InteriorPointSolver::_minimisation_step(c,xl,xu,A,b, x,y,zl,zu)\n");
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
        if(xl[i]!=-inf) { rzl[i] = (x[i]-xl[i])*zl[i] - mu; } else { rzl[i]=0.0; }
        if(xu[i]!=+inf) { rzu[i] = (xu[i]-x[i])*zu[i] - mu; } else { rzu[i]=0.0; }
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
    S=adat(A,D.diagonal());;
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
    ARIADNE_LOG(5,"cx="<<dot(c,x)<<" yb="<<dot(y,b)<<" Ax-b="<<(A*x-b)<<" yA+(zl-zu)-c="<<(y*A+(zl-zu)-c)<<"\n");

    if(alphax==1.0 && alphaz==1.0) { return PRIMAL_DUAL_FEASIBLE; }
    if(alphax==1.0) { return PRIMAL_FEASIBLE; }
    if(alphaz==1.0) { return DUAL_FEASIBLE; }
    return INDETERMINATE_FEASIBILITY;
}



LinearProgramStatus InteriorPointSolver::
_feasibility_step(const Vector<Float>& xl, const Vector<Float>& xu, const Matrix<Float>& A, const Vector<Float>& b,
                  Vector<Float>& x, Vector<Float>& y, Vector<Float>& zl, Vector<Float>& zu) const
{
    Vector<Float> c(x.size(),0.0);
    return this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
}



} // namespace Ariadne

