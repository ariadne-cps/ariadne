/***************************************************************************
 *            linear_programming.cc
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "utility/tuple.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/diagonal_matrix.h"
#include "function/affine.h"
#include "solvers/linear_programming.h"

#include "utility/macros.h"
#include "utility/logging.h"

static const int verbosity=0;

namespace Ariadne {

typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Matrix<UpperIntervalType> UpperIntervalMatrixType;


// Return r[i]=x[i]-c for i=1,...,n
Vector<Float64> esub(const Vector<Float64>& x, const Float64& c) {
    Vector<Float64> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]-c;
    }
    return r;
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<Float64> emul(const Vector<Float64>& x, const Vector<Float64>& y) {
    Vector<Float64> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i];
    }
    return r;
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<Float64> ediv(const Vector<Float64>& x, const Vector<Float64>& z) {
    Vector<Float64> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]/z[i];
    }
    return r;
}

// Return r[i]=x[i]*y[i]+z for i=1,...,n
Vector<Float64> efma(const Vector<Float64>& x, const Vector<Float64>& y, const Float64& z) {
    Vector<Float64> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i]+z;
    }
    return r;
}

// Return r[i]=1/y[i] for i=1,...,n
Vector<Float64> erec(const Vector<Float64>& v) {
    Vector<Float64> r(v.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=1.0/v[i];
    }
    return r;
}

Bool is_nan(Vector<Float64> const& v) {
    for(Nat i=0; i!=v.size(); ++i) {
        if(is_nan(v[i])) { return true; }
    }
    return false;
}


// Compute R=ADA^T for diagonal D
Matrix<Float64> adat(const Matrix<Float64>& A, const Vector<Float64>& D)
{
    Matrix<Float64> R(A.row_size(),A.row_size());
    ARIADNE_ASSERT(D.size()==A.column_size());
    for(Nat i=0; i!=A.row_size(); ++i) {
        for(Nat k=0; k!=A.column_size(); ++k) {
            Float64 ADik=A[i][k]*D[k];
            for(Nat j=0; j!=A.row_size(); ++j) {
                R[i][j]+=ADik*A[j][k];
            }
        }
    }
    return R;
}

Bool all_greater(const Vector<Float64>& x, const Float64& e) {
    for(Nat i=0; i!=x.size(); ++i) {
        if(x[i]<=e) { return false; }
    }
    return true;
}

Bool all_less(const Vector<Float64>& x, const Float64& e) {
    for(Nat i=0; i!=x.size(); ++i) {
        if(x[i]>=e) { return false; }
    }
    return true;
}

Bool all_greater(const Vector<Float64>& x1, const Vector<Float64>& x2) {
    for(Nat i=0; i!=x1.size(); ++i) {
        if(x1[i]<=x2[i]) { return false; }
    }
    return true;
}


Float64 compute_mu(const Vector<Float64>& xl, const Vector<Float64>& xu,
                 const Vector<Float64>& x, const Vector<Float64>& zl, const Vector<Float64>& zu)
{
    const Nat n=x.size();
    Float64 mu = 0.0;
    for(Nat i=0; i!=n; ++i) {
        if(xl[i]!=-inf) { mu += ((x[i]-xl[i])*zl[i]); }
        if(xu[i]!=+inf) { mu += ((xu[i]-x[i])*zu[i]); }
    }
    mu /= (2.0*n);
    return mu;
}



inline Float64Bounds mul_val(Float64 x1, Float64 x2) { return Float64Bounds(mul_down(x1,x2),mul_up(x1,x2)); }
inline Vector<Float64Value> const& cast_exact(Vector<Float64> const& v) { return reinterpret_cast<Vector<Float64Value>const&>(v); }
inline Matrix<Float64Value> const& cast_exact(Matrix<Float64> const& A) { return reinterpret_cast<Matrix<Float64Value>const&>(A); }


ValidatedKleenean InteriorPointSolver::
validate_feasibility(const Vector<Float64>& xl, const Vector<Float64>& xu,
                     const Matrix<Float64>& A, const Vector<Float64>& b,
                     const Vector<Float64>& x, const Vector<Float64>& y) const
{
    const Nat m=A.row_size();
    const Nat n=A.column_size();

    // x should be an approximate solution to Ax=b
    // Use the fact that for any x, x'=(x + A^T (AA^T)^{-1}(b-Ax)) satisfies Ax'=0
    Vector<ValidatedNumericType> ivlx = cast_exact(x);
    Vector<ValidatedNumericType> ivle = cast_exact(b)-cast_exact(A)*ivlx;

    Matrix<ValidatedNumericType> ivlS(m,m);
    for(Nat i1=0; i1!=m; ++i1) {
        for(Nat i2=i1; i2!=m; ++i2) {
            for(Nat j=0; j!=n; ++j) {
                ivlS[i1][i2] += mul_val(A[i1][j],A[i2][j]);
            }
        }
    }
    for(Nat i1=0; i1!=m; ++i1) {
        for(Nat i2=0; i2!=i1; ++i2) {
            ivlS[i1][i2]=ivlS[i2][i1];
        }
    }

    Vector<ValidatedNumericType> ivld =  transpose(cast_exact(A)) * solve(ivlS,ivle);

    ivlx += ivld;

    ARIADNE_LOG(2,"[x] = "<<ivlx);
    ValidatedKleenean result=true;
    for(Nat i=0; i!=n; ++i) {
        if(ivlx[i].lower().raw()<=xl[i] || ivlx[i].upper().raw()>=xu[i]) {
            result = indeterminate; break;
        }
    }

    if(definitely(result)) { return result; }

    // If yb - max(yA,0) xu + min(yA,0) xl > 0, then problem is infeasible
    // Evaluate lower bound for yb - max(z,0) xu + min(z,0) xl
    Vector<ValidatedNumericType> z=transpose(cast_exact(A)) * cast_exact(y);
    Float64 mx = 0.0;
    Float64::set_rounding_downward();
    for(Nat i=0; i!=y.size(); ++i) {
        mx += (b[i]*y[i]).raw();
    }
    for(Nat i=0; i!=x.size(); ++i) {
        Float64 neg_xui = -xu[i];
        if(z[i].upper().raw()>0.0) { mx += z[i].upper().raw() * neg_xui; }
        if(z[i].lower().raw()<0.0) { mx += z[i].lower().raw() * xl[i]; }
    }
    Float64::set_rounding_to_nearest();
    if(mx>0.0) { return false; }

    return indeterminate;
}


Tuple< Float64, Vector<Float64>, Vector<Float64> >
InteriorPointSolver::minimise(const Vector<Float64>& c,
                              const Vector<Float64>& xl, const Vector<Float64>& xu,
                              const Matrix<Float64>& A, const Vector<Float64>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::minimise(c,xl,xu,A,b)\n");
    ARIADNE_LOG(2,"A="<<A<<", b="<<b<<", c="<<c<<"\n");
    ARIADNE_LOG(2,"xl="<<xl<<", xu="<<xu<<"\n");

    const Nat m = b.size();
    const Nat n = c.size();
    Vector<Float64> y(m, 0.0);
    Vector<Float64> x(n);
    Vector<Float64> zl(n);
    Vector<Float64> zu(n);
    for(Nat i=0; i!=n; ++i) {
        if(xl[i]==-inf) {
            if(xu[i]==+inf) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+inf) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-inf) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+inf) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }

    Bool done;
    do {
        done = this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
    } while(!done);

    do {
        this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
    } while(dot(c,x)-dot(y,b)>1e-4);

    // Todo: check for optimality
    return make_tuple(dot(c,x),x,y);
}



Tuple< Float64, Vector<Float64>, Vector<Float64> >
InteriorPointSolver::
hotstarted_minimise(const Vector<Float64>& c,
                    const Vector<Float64>& xl, const Vector<Float64>& xu,
                    const Matrix<Float64>& A, const Vector<Float64>& b,
                    Vector<Float64>& x, Vector<Float64>& y, Vector<Float64>& zl, Vector<Float64>& zu) const
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
    const Nat maxsteps=10;

    Float64 cx,yb;
    cx=dot(c,x);
    yb=dot(y,b);
    ARIADNE_ASSERT(yb<=cx);

    ARIADNE_LOG(2,"xl="<<xl<<" xu="<<xu<<" A="<<A<<" b="<<b<<" c="<<c<<"\n");
    ARIADNE_LOG(2,"x="<<x<<" y="<<y<<" z="<<zl<<" zu="<<zu<<"\n");

    Nat steps=0;
    while(steps<maxsteps && (cx-yb)>maxerror) {
        this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
        ++steps;
    }

    return make_tuple(dot(c,x),x,y);

}


ValidatedKleenean
InteriorPointSolver::
feasible(const Vector<Float64>& xl, const Vector<Float64>& xu,
         const Matrix<Float64>& A, const Vector<Float64>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::feasible(xl,xu,A,b)\n");
    ARIADNE_LOG(3,"A="<<A<<", b="<<b<<"\n");
    ARIADNE_LOG(3,"xl="<<xl<<", xu="<<xu<<"\n");

    const Nat m = A.row_size();
    const Nat n = A.column_size();
    Vector<Float64> c(n, 0.0);
    Vector<Float64> y(m, 0.0);
    Vector<Float64> x(n);
    Vector<Float64> zl(n);
    Vector<Float64> zu(n);
    for(Nat i=0; i!=n; ++i) {
        if(xl[i]==-inf) {
            if(xu[i]==+inf) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+inf) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<=xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-inf) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+inf) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }
    Vector<Float64Bounds> X(n);
    for(SizeType i=0; i!=n; ++i) {
        X[i]=Float64Bounds(xl[i],xu[i]);
    }

    const double THRESHOLD = 1e-8;
    Nat step=0;
    while(step++<24) {
        LinearProgramStatus result=this->_feasibility_step(xl,xu,A,b, x,y,zl,zu);
        if(result==PRIMAL_DUAL_FEASIBLE || result==PRIMAL_FEASIBLE) {
            ValidatedKleenean validated_feasible=this->validate_feasibility(xl,xu,A,b, x,y);
            if(definitely(validated_feasible)) { return true; }
        }
        Float64Bounds yb=dot(Vector<Float64Bounds>(y),Vector<Float64Bounds>(b));
        // NOTE: Must compute y*A first, as A*X may give NaN.
        Float64Bounds yAX = dot( transpose(Matrix<Float64Bounds>(A)) * Vector<Float64Bounds>(y), X );
        if(inconsistent(yb,yAX)) { return false; }
        if(result==DEGENERATE_FEASIBILITY) { ARIADNE_LOG(2,"  degenerate\n"); return indeterminate; }
        if(compute_mu(xl,xu, x,zl,zu)<THRESHOLD ) { ARIADNE_LOG(2,"  threshold\n"); return indeterminate; }
    }
    return indeterminate;
}












LinearProgramStatus InteriorPointSolver::
_minimisation_step(const Vector<Float64>& c, const Vector<Float64>& xl, const Vector<Float64>& xu, const Matrix<Float64>& A, const Vector<Float64>& b,
                   Vector<Float64>& x, Vector<Float64>& y, Vector<Float64>& zl, Vector<Float64>& zu) const
{
    ARIADNE_LOG(4,"InteriorPointSolver::_minimisation_step(c,xl,xu,A,b, x,y,zl,zu)\n");
    ARIADNE_LOG(4,"x="<<x<<", y="<<y<<", zl="<<zl<<", zu="<<zu<<"\n");

    static const double gamma=1.0/256;
    static const double sigma=1.0/8;
    static const double scale=0.75;


    const Nat m=A.row_size();
    const Nat n=A.column_size();

    Vector<Float64> dx(n),dy(m),dzl(n),dzu(n);
    Vector<Float64> nx,ny,nzl,nzu;
    Vector<Float64> rx(m),ry(n),rzl(n),rzu(n);
    Matrix<Float64> S(m,m);
    DiagonalMatrix<Float64> Xl(xl), Xu(xu), X(x), Zl(zl), Zu(zu);

    Float64 mu = compute_mu(xl,xu, x,zl,zu) * sigma;
    mu=1.0;
    ARIADNE_LOG(4,"mu="<<mu<<"\n");

    // rx = Ax-b; ry=yA+zl-zu-c; rzl=(x-xl).zl-mu; rzu=(xu-x).zu-mu.
    rx=A*x-b;
    ry=transpose(A)*y+(zl-zu)-c;
    for(Nat i=0; i!=n; ++i) {
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
    DiagonalMatrix<Float64> D(erec(ediv(zu,xu-x)+ediv(zl,x-xl)));
    Vector<Float64> ryz = ry - ediv(rzl,Vector<Float64>(x-xl)) + ediv(rzu,Vector<Float64>(xu-x));
    S=adat(A,D.diagonal());;
    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<inverse(S)<<"\n");

    // dzl = (rzl - Zl dx) / (X-Xl)
    // dzu = (rzu + Zu dx) / (Xu-X)
    // dx = D (AT dy - ryz)
    // (A D AT) dy = rx + A D ryz

    dy = solve(S, Vector<Float64>( rx + A * (D * ryz) ) );
    if(is_nan(dy)) { return DEGENERATE_FEASIBILITY; }
    dx = D * Vector<Float64>(transpose(A)*dy - ryz);
    dzl = Vector<Float64>(rzl-Zl*dx)/(X-Xl);
    dzu = Vector<Float64>(rzu+Zu*dx)/(Xu-X);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<" dzl="<<dzl<<" dzu="<<dzu<<"\n");

    ARIADNE_LOG(7,"A*dx="<<(A*dx)<<" AT*dy+dzl-dzu="<<(transpose(A)*dy+dzl-dzu)<<" Zl*dx+(X-Xl)*dzl="<<(Zl*dx+(X-Xl)*dzl)<<" -Zu*dx+(Xu-X)*dzu="<<(-(Zu*dx)+(Xu-X)*dzu)<<"\n");
    ARIADNE_LOG(7,"A*dx-rx="<<(A*dx-rx)<<" AT*dy+dzl-dzu-ry="<<(transpose(A)*dy+dzl-dzu-ry)<<" Zl*dx+(X-Xl)*dzl-rzl="<<(Zl*dx+(X-Xl)*dzl-rzl)<<" -Zu*dx+(Xu-X)*dzu-rzu="<<(-(Zu*dx)+(Xu-X)*dzu-rzu)<<"\n");
    // Try to enforce feasibility or dual feasibility
    Float64 alphax=1.0;
    nx=x-dx;
    while ( !all_greater(emul(nx-xl,zl),gamma*mu) || !all_greater(emul(xu-nx,zu),gamma*mu) ) {
        alphax=alphax*scale;
        nx=(x-alphax*dx);
        if(alphax<gamma*mu/4096) { return DEGENERATE_FEASIBILITY; }
    }

    Float64 alphaz=1.0;
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
    ARIADNE_LOG(5,"cx="<<dot(c,x)<<" yb="<<dot(y,b)<<" Ax-b="<<(A*x-b)<<" yA+(zl-zu)-c="<<(transpose(A)*y+(zl-zu)-c)<<"\n");

    if(alphax==1.0 && alphaz==1.0) { return PRIMAL_DUAL_FEASIBLE; }
    if(alphax==1.0) { return PRIMAL_FEASIBLE; }
    if(alphaz==1.0) { return DUAL_FEASIBLE; }
    return INDETERMINATE_FEASIBILITY;
}



LinearProgramStatus InteriorPointSolver::
_feasibility_step(const Vector<Float64>& xl, const Vector<Float64>& xu, const Matrix<Float64>& A, const Vector<Float64>& b,
                  Vector<Float64>& x, Vector<Float64>& y, Vector<Float64>& zl, Vector<Float64>& zu) const
{
    Vector<Float64> c(x.size(),0.0);
    return this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
}



} // namespace Ariadne

