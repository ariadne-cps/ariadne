/***************************************************************************
 *            solvers/linear_programming.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../utility/tuple.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/diagonal_matrix.hpp"
#include "../function/affine.hpp"
#include "../solvers/linear_programming.hpp"

#include "../utility/macros.hpp"
#include "../output/logging.hpp"

namespace Ariadne {

typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Matrix<UpperIntervalType> UpperIntervalMatrixType;

// Return r[i]=x[i]-c for i=1,...,n
Vector<FloatDP> esub(const Vector<FloatDP>& x, const FloatDP& c);
// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<FloatDP> emul(const Vector<FloatDP>& x, const Vector<FloatDP>& y);
// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<FloatDP> ediv(const Vector<FloatDP>& x, const Vector<FloatDP>& z);
// Return r[i]=x[i]*y[i]+z for i=1,...,n
Vector<FloatDP> efma(const Vector<FloatDP>& x, const Vector<FloatDP>& y, const FloatDP& z);
// Return r[i]=1/y[i] for i=1,...,n
Vector<FloatDP> erec(const Vector<FloatDP>& v);

Bool is_nan(Vector<FloatDP> const& v);

// Compute R=ADA^T for diagonal D
Matrix<FloatDP> adat(const Matrix<FloatDP>& A, const Vector<FloatDP>& D);

Bool all_greater(const Vector<FloatDP>& x, const FloatDP& e);
Bool all_less(const Vector<FloatDP>& x, const FloatDP& e);
Bool all_greater(const Vector<FloatDP>& x1, const Vector<FloatDP>& x2);

FloatDP compute_mu(const Vector<FloatDP>& xl, const Vector<FloatDP>& xu,
                 const Vector<FloatDP>& x, const Vector<FloatDP>& zl, const Vector<FloatDP>& zu);


// Return r[i]=x[i]-c for i=1,...,n
Vector<FloatDP> esub(const Vector<FloatDP>& x, const FloatDP& c) {
    Vector<FloatDP> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]-c;
    }
    return r;
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<FloatDP> emul(const Vector<FloatDP>& x, const Vector<FloatDP>& y) {
    Vector<FloatDP> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i];
    }
    return r;
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<FloatDP> ediv(const Vector<FloatDP>& x, const Vector<FloatDP>& z) {
    Vector<FloatDP> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]/z[i];
    }
    return r;
}

// Return r[i]=x[i]*y[i]+z for i=1,...,n
Vector<FloatDP> efma(const Vector<FloatDP>& x, const Vector<FloatDP>& y, const FloatDP& z) {
    Vector<FloatDP> r(x.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i]+z;
    }
    return r;
}

// Return r[i]=1/y[i] for i=1,...,n
Vector<FloatDP> erec(const Vector<FloatDP>& v) {
    Vector<FloatDP> r(v.size());
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=1.0/v[i];
    }
    return r;
}

Bool is_nan(Vector<FloatDP> const& v) {
    for(Nat i=0; i!=v.size(); ++i) {
        if(is_nan(v[i])) { return true; }
    }
    return false;
}


// Compute R=ADA^T for diagonal D
Matrix<FloatDP> adat(const Matrix<FloatDP>& A, const Vector<FloatDP>& D)
{
    Matrix<FloatDP> R(A.row_size(),A.row_size());
    ARIADNE_ASSERT(D.size()==A.column_size());
    for(Nat i=0; i!=A.row_size(); ++i) {
        for(Nat k=0; k!=A.column_size(); ++k) {
            FloatDP ADik=A[i][k]*D[k];
            for(Nat j=0; j!=A.row_size(); ++j) {
                R[i][j]+=ADik*A[j][k];
            }
        }
    }
    return R;
}

Bool all_greater(const Vector<FloatDP>& x, const FloatDP& e) {
    for(Nat i=0; i!=x.size(); ++i) {
        if(x[i]<=e) { return false; }
    }
    return true;
}

Bool all_less(const Vector<FloatDP>& x, const FloatDP& e) {
    for(Nat i=0; i!=x.size(); ++i) {
        if(x[i]>=e) { return false; }
    }
    return true;
}

Bool all_greater(const Vector<FloatDP>& x1, const Vector<FloatDP>& x2) {
    for(Nat i=0; i!=x1.size(); ++i) {
        if(x1[i]<=x2[i]) { return false; }
    }
    return true;
}


FloatDP compute_mu(const Vector<FloatDP>& xl, const Vector<FloatDP>& xu,
                 const Vector<FloatDP>& x, const Vector<FloatDP>& zl, const Vector<FloatDP>& zu)
{
    const Nat n=x.size();
    FloatDP mu = 0.0;
    for(Nat i=0; i!=n; ++i) {
        if(xl[i]!=-inf) { mu += ((x[i]-xl[i])*zl[i]); }
        if(xu[i]!=+inf) { mu += ((xu[i]-x[i])*zu[i]); }
    }
    mu /= (2.0*n);
    return mu;
}



inline FloatDPBounds mul_val(FloatDP x1, FloatDP x2) { return FloatDPBounds(mul(down,x1,x2),mul(up,x1,x2)); }
inline Vector<FloatDPValue> const& cast_exact(Vector<FloatDP> const& v) { return reinterpret_cast<Vector<FloatDPValue>const&>(v); }
inline Matrix<FloatDPValue> const& cast_exact(Matrix<FloatDP> const& A) { return reinterpret_cast<Matrix<FloatDPValue>const&>(A); }


ValidatedKleenean InteriorPointSolver::
validate_feasibility(const Vector<FloatDP>& axl, const Vector<FloatDP>& axu,
                     const Matrix<FloatDP>& aA, const Vector<FloatDP>& ab,
                     const Vector<FloatDP>& ax, const Vector<FloatDP>& ay) const
{
    Vector<FloatDPValue> const& xl = cast_exact(axl);
    Vector<FloatDPValue> const& xu = cast_exact(axu);
    Matrix<FloatDPValue> const& A = cast_exact(aA);
    Vector<FloatDPValue> const& b = cast_exact(ab);

    Vector<FloatDPBounds> x = cast_exact(ax);
    Vector<FloatDPBounds> y = cast_exact(ay);

    FloatDPValue zero;

    const Nat n=A.column_size();

    // x should be an approximate solution to Ax=b
    // Use the fact that for any x, x'=(x + A^T (AA^T)^{-1}(b-Ax)) satisfies Ax'=0
    Vector<FloatDPBounds> e = b-A*x;

    Matrix<FloatDPBounds> S=A*transpose(A);
    Vector<FloatDPBounds> d =  transpose(A) * solve(S,e);

    x += d;

    ARIADNE_LOG(2,"[x] = "<<x);
    ValidatedKleenean result=true;
    for(Nat i=0; i!=n; ++i) {
        if(x[i].lower().raw()<=xl[i].raw() || x[i].upper().raw()>=xu[i].raw()) {
            result = indeterminate; break;
        }
    }

    if(definitely(result)) { return result; }

    // If yb - max(yA,0) xu + min(yA,0) xl > 0, then problem is infeasible
    // Evaluate lower bound for yb - max(z,0) xu + min(z,0) xl
    Vector<FloatDPBounds> z=transpose(A) * y;
    FloatDPLowerBound mx = zero;
    for(Nat i=0; i!=y.size(); ++i) {
        mx += (b[i]*y[i]);
    }
    for(Nat i=0; i!=x.size(); ++i) {
        FloatDPValue zil=cast_exact(z[i].lower());
        FloatDPValue ziu=cast_exact(z[i].upper());
        if(ziu>0) { mx -= ziu * xu[i]; }
        if(zil<0) { mx += zil * xl[i]; }
    }
    if(definitely(mx>0)) { return false; }

    return indeterminate;
}


Tuple< FloatDP, Vector<FloatDP>, Vector<FloatDP> >
InteriorPointSolver::minimise(const Vector<FloatDP>& c,
                              const Vector<FloatDP>& xl, const Vector<FloatDP>& xu,
                              const Matrix<FloatDP>& A, const Vector<FloatDP>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::minimise(c,xl,xu,A,b)\n");
    ARIADNE_LOG(2,"A="<<A<<", b="<<b<<", c="<<c<<"\n");
    ARIADNE_LOG(2,"xl="<<xl<<", xu="<<xu<<"\n");

    const Nat m = b.size();
    const Nat n = c.size();
    Vector<FloatDP> y(m, 0.0);
    Vector<FloatDP> x(n);
    Vector<FloatDP> zl(n);
    Vector<FloatDP> zu(n);
    for(Nat i=0; i!=n; ++i) {
        if(xl[i]==-inf) {
            if(xu[i]==+inf) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+inf) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-inf) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+inf) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }

    LinearProgramStatus status = LinearProgramStatus::INDETERMINATE_FEASIBILITY;
    while(status == LinearProgramStatus::INDETERMINATE_FEASIBILITY) {
        status = this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
    }

    do {
        this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
    } while(dot(c,x)-dot(y,b)>1e-4);

    // Todo: check for optimality
    return make_tuple(dot(c,x),x,y);
}



Tuple< FloatDP, Vector<FloatDP>, Vector<FloatDP> >
InteriorPointSolver::
hotstarted_minimise(const Vector<FloatDP>& c,
                    const Vector<FloatDP>& xl, const Vector<FloatDP>& xu,
                    const Matrix<FloatDP>& A, const Vector<FloatDP>& b,
                    Vector<FloatDP>& x, Vector<FloatDP>& y, Vector<FloatDP>& zl, Vector<FloatDP>& zu) const
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

    FloatDP cx,yb;
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
feasible(const Vector<FloatDP>& xl, const Vector<FloatDP>& xu,
         const Matrix<FloatDP>& A, const Vector<FloatDP>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::feasible(xl,xu,A,b)\n");
    ARIADNE_LOG(3,"A="<<A<<", b="<<b<<"\n");
    ARIADNE_LOG(3,"xl="<<xl<<", xu="<<xu<<"\n");

    const Nat m = A.row_size();
    const Nat n = A.column_size();
    Vector<FloatDP> c(n, 0.0);
    Vector<FloatDP> y(m, 0.0);
    Vector<FloatDP> x(n);
    Vector<FloatDP> zl(n);
    Vector<FloatDP> zu(n);
    for(Nat i=0; i!=n; ++i) {
        if(xl[i]==-inf) {
            if(xu[i]==+inf) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+inf) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<=xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-inf) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+inf) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }
    Vector<FloatDPBounds> X(n);
    for(SizeType i=0; i!=n; ++i) {
        X[i]=FloatDPBounds(xl[i],xu[i]);
    }

    const double THRESHOLD = 1e-8;
    Nat step=0;
    while(step++<24) {
        LinearProgramStatus result=this->_feasibility_step(xl,xu,A,b, x,y,zl,zu);
        if(result==LinearProgramStatus::PRIMAL_DUAL_FEASIBLE || result==LinearProgramStatus::PRIMAL_FEASIBLE) {
            ValidatedKleenean validated_feasible=this->validate_feasibility(xl,xu,A,b, x,y);
            if(definitely(validated_feasible)) { return true; }
        }
        FloatDPBounds yb=dot(Vector<FloatDPBounds>(y),Vector<FloatDPBounds>(b));
        // NOTE: Must compute y*A first, as A*X may give NaN.
        FloatDPBounds yAX = dot( transpose(Matrix<FloatDPBounds>(A)) * Vector<FloatDPBounds>(y), X );
        if(inconsistent(yb,yAX)) { return false; }
        if(result==LinearProgramStatus::DEGENERATE_FEASIBILITY) { ARIADNE_LOG(2,"  degenerate\n"); return indeterminate; }
        if(compute_mu(xl,xu, x,zl,zu)<THRESHOLD ) { ARIADNE_LOG(2,"  threshold\n"); return indeterminate; }
    }
    return indeterminate;
}












LinearProgramStatus InteriorPointSolver::
_minimisation_step(const Vector<FloatDP>& c, const Vector<FloatDP>& xl, const Vector<FloatDP>& xu, const Matrix<FloatDP>& A, const Vector<FloatDP>& b,
                   Vector<FloatDP>& x, Vector<FloatDP>& y, Vector<FloatDP>& zl, Vector<FloatDP>& zu) const
{
    ARIADNE_LOG(4,"InteriorPointSolver::_minimisation_step(c,xl,xu,A,b, x,y,zl,zu)\n");
    ARIADNE_LOG(4,"x="<<x<<", y="<<y<<", zl="<<zl<<", zu="<<zu<<"\n");

    static const double gamma=1.0/256;
    static const double sigma=1.0/8;
    static const double scale=0.75;


    const Nat m=A.row_size();
    const Nat n=A.column_size();

    Vector<FloatDP> dx(n),dy(m),dzl(n),dzu(n);
    Vector<FloatDP> nx,ny,nzl,nzu;
    Vector<FloatDP> rx(m),ry(n),rzl(n),rzu(n);
    Matrix<FloatDP> S(m,m);
    DiagonalMatrix<FloatDP> Xl(xl), Xu(xu), X(x), Zl(zl), Zu(zu);

    FloatDP mu = compute_mu(xl,xu, x,zl,zu) * sigma;
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
    DiagonalMatrix<FloatDP> D(erec(ediv(zu,xu-x)+ediv(zl,x-xl)));
    Vector<FloatDP> ryz = ry - ediv(rzl,Vector<FloatDP>(x-xl)) + ediv(rzu,Vector<FloatDP>(xu-x));
    S=adat(A,D.diagonal());
    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<inverse(S)<<"\n");

    // dzl = (rzl - Zl dx) / (X-Xl)
    // dzu = (rzu + Zu dx) / (Xu-X)
    // dx = D (AT dy - ryz)
    // (A D AT) dy = rx + A D ryz

    dy = solve(S, Vector<FloatDP>( rx + A * (D * ryz) ) );
    if(is_nan(dy)) { return LinearProgramStatus::DEGENERATE_FEASIBILITY; }
    dx = D * Vector<FloatDP>(transpose(A)*dy - ryz);
    dzl = Vector<FloatDP>(rzl-Zl*dx)/(X-Xl);
    dzu = Vector<FloatDP>(rzu+Zu*dx)/(Xu-X);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<" dzl="<<dzl<<" dzu="<<dzu<<"\n");

    ARIADNE_LOG(7,"A*dx="<<(A*dx)<<" AT*dy+dzl-dzu="<<(transpose(A)*dy+dzl-dzu)<<" Zl*dx+(X-Xl)*dzl="<<(Zl*dx+(X-Xl)*dzl)<<" -Zu*dx+(Xu-X)*dzu="<<(-(Zu*dx)+(Xu-X)*dzu)<<"\n");
    ARIADNE_LOG(7,"A*dx-rx="<<(A*dx-rx)<<" AT*dy+dzl-dzu-ry="<<(transpose(A)*dy+dzl-dzu-ry)<<" Zl*dx+(X-Xl)*dzl-rzl="<<(Zl*dx+(X-Xl)*dzl-rzl)<<" -Zu*dx+(Xu-X)*dzu-rzu="<<(-(Zu*dx)+(Xu-X)*dzu-rzu)<<"\n");
    // Try to enforce feasibility or dual feasibility
    FloatDP alphax=1.0;
    nx=x-dx;
    while ( !all_greater(emul(nx-xl,zl),gamma*mu) || !all_greater(emul(xu-nx,zu),gamma*mu) ) {
        alphax=alphax*scale;
        nx=(x-alphax*dx);
        if(alphax<gamma*mu/4096) { return LinearProgramStatus::DEGENERATE_FEASIBILITY; }
    }

    FloatDP alphaz=1.0;
    nzl=zl-dzl; nzu=zu-dzu;
    while ( !all_greater(emul(nx-xl,nzl),gamma*mu) || !all_greater(emul(xu-nx,nzu),gamma*mu) ) {
        alphaz=alphaz*scale;
        nzl=(zl-alphaz*dzl);
        nzu=(zu-alphaz*dzu);
        if(alphaz<gamma*mu/4096) { return LinearProgramStatus::DEGENERATE_FEASIBILITY; }
        //ARIADNE_LOG(9,"alphaz="<<alphaz<<" nzl="<<nzl<<" nzu="<<nzu<<"\n");
    }
    ny=(y-alphaz*dy);
    ARIADNE_LOG(5,"alphax="<<alphax<<" nx="<<nx<<" alphaz="<<alphaz<<" ny="<<ny<<" nzl="<<nzl<<" nzu="<<nzu<<"\n");

    x=nx; y=ny; zl=nzl; zu=nzu;
    ARIADNE_LOG(5,"cx="<<dot(c,x)<<" yb="<<dot(y,b)<<" Ax-b="<<(A*x-b)<<" yA+(zl-zu)-c="<<(transpose(A)*y+(zl-zu)-c)<<"\n");

    if(alphax==1.0 && alphaz==1.0) { return LinearProgramStatus::PRIMAL_DUAL_FEASIBLE; }
    if(alphax==1.0) { return LinearProgramStatus::PRIMAL_FEASIBLE; }
    if(alphaz==1.0) { return LinearProgramStatus::DUAL_FEASIBLE; }
    return LinearProgramStatus::INDETERMINATE_FEASIBILITY;
}



LinearProgramStatus InteriorPointSolver::
_feasibility_step(const Vector<FloatDP>& xl, const Vector<FloatDP>& xu, const Matrix<FloatDP>& A, const Vector<FloatDP>& b,
                  Vector<FloatDP>& x, Vector<FloatDP>& y, Vector<FloatDP>& zl, Vector<FloatDP>& zu) const
{
    Vector<FloatDP> c(x.size(),0.0);
    return this->_minimisation_step(c,xl,xu,A,b, x,y,zl,zu);
}



} // namespace Ariadne

