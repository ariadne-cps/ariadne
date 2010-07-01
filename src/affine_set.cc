/***************************************************************************
 *            affine_set.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "linear_programming.h"
#include "function.h"
#include "affine.h"

#include "box.h"
#include "grid_set.h"
#include "affine_set.h"

#include "graphics_interface.h"


namespace Ariadne {

class DegenerateFeasibilityProblemException : public std::exception { };


typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;


template<class X>
struct LinearProgram {
    Matrix<X> A;
    Vector<X> b;
    Vector<X> c;
    Vector<X> l;
    Vector<X> u;
    array<Slackness> vt;
    array<size_t> p;
    Matrix<X> B;
    Vector<X> x;
    Vector<X> y;
};




AffineSet::AffineSet(const Vector<Interval>& D, const Matrix<Float>& G, const Vector<Float>& h)
{
    this->construct(D,G,h);
}

AffineSet::AffineSet(const Matrix<Float>& G, const Vector<Float>& h)
{
    this->construct(Vector<Interval>(G.column_size(),Interval(-1,+1)),G,h);
}

void AffineSet::construct(const Vector<Interval>& D, const Matrix<Float>& G, const Vector<Float>& h)
{
    ARIADNE_ASSERT_MSG(G.row_size()==h.size() && G.row_size()>0,"G="<<G<<", h="<<h);
    this->_domain=D;
    for(uint i=0; i!=G.row_size(); ++i) {
        Affine<Float> x(G.column_size());
        x=h[i];
        for(uint j=0; j!=G.column_size(); ++j) {
            x[j]=G[i][j];
        }
        this->_function.append(x);
    }
}

void
AffineSet::new_inequality_constraint(const Vector<Float>& a, const Float& b)
{
    // Re-write the constraint ax<=b as b-ax>=0
    ARIADNE_ASSERT_MSG(a.size()==this->_function[0].argument_size(),"a="<<a<<" f="<<this->_function);
    Affine<Float> c(a.size());
    c=b;
    for(uint j=0; j!=a.size(); ++j) {
        c[j]=-a[j];
    }
    _constraints.append(c);
}

void
AffineSet::new_equality_constraint(const Vector<Float>& a, const Float& b)
{
    // Re-write the constraint ax<=b as b-ax>=0
    ARIADNE_ASSERT_MSG(a.size()==this->_function[0].argument_size(),"a="<<a<<" f="<<this->_function);
    Affine<Float> c(a.size());
    c=b;
    for(uint j=0; j!=a.size(); ++j) {
        c[j]=-a[j];
    }
    _equations.append(c);
}


AffineSet*
AffineSet::clone() const
{
    return new AffineSet(*this);
}


uint
AffineSet::dimension() const
{
    return this->_function.size();
}

uint
AffineSet::number_of_parameters() const
{
    return this->_function[0].argument_size();
}

uint
AffineSet::number_of_constraints() const
{
    return this->_constraints.size();
}

Vector<Interval>
AffineSet::domain() const
{
    return this->_domain;
}

tribool AffineSet::bounded() const {
    return Box(this->domain()).bounded() || indeterminate;
}

Box AffineSet::bounding_box() const {
    Box result(this->dimension());
    Vector<Interval> domain=this->domain();
    for(uint i=0; i!=this->dimension(); ++i) {
        //result[i]=evaluate(this->_function[i],domain);
        result[i]=this->_function[i].evaluate(domain);
    }
    return result;
}


tribool AffineSet::disjoint(const Box& bx) const {
    LinearProgram<Float> lp;
    this->construct_linear_program(lp);
    tribool feasible;
    try {
        feasible=SimplexSolver<Float>().constrained_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);
    }
    catch(DegenerateFeasibilityProblemException e) {
        feasible=indeterminate;
    }
    return !feasible;
}

tribool AffineSet::empty() const {
    return this->disjoint(this->bounding_box());
}

template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& x) {
    return prod(A,x);
}


GridTreeSet
AffineSet::outer_approximation(const Grid& g, int d) const {
    GridTreeSet r(g);
    this->adjoin_outer_approximation_to(r,d);
    return r;
}


void _adjoin_outer_approximation_to(GridTreeSet& paving, LinearProgram<Float>& lp, GridCell& cell, int depth)
{

    // No need to check if cell is already part of the set
    if(subset(cell,paving)) {
        return;
    }

    // Find concrete cell box
    const Box& bx=cell.box();

    // Make part of linear program dependent on cell
    for(uint i=0; i!=cell.dimension(); ++i) {
        lp.l[i]=bx[i].lower();
        lp.u[i]=bx[i].upper();
    }

    int cell_tree_depth=(cell.depth()-cell.height());
    int maximum_tree_depth=depth*cell.dimension();

    // Check for disjointness using linear program
    tribool feasible=SimplexSolver<Float>().constrained_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);
    //feasible=verify_constrained_feasibility(lp.A,lp.b,lp.l,lp.u,lp.vt);
    if(!feasible) { return; }

    // If not disjoint, and we are at maximum depth, adjoin the set.
    // Otherwise, split and continue.

    // FIXME: The cell depth should be an absolute depth, so this conversion should not be needed
    if(cell_tree_depth>=maximum_tree_depth) {
        paving.adjoin(cell);
    } else {
        GridCell subcell1 = cell.split(0);
        GridCell subcell2 = cell.split(1);
        _adjoin_outer_approximation_to(paving,lp,subcell1,depth);
        _adjoin_outer_approximation_to(paving,lp,subcell2,depth);
    }

}



void
AffineSet::construct_linear_program(LinearProgram<Float>& lp) const
{
    // Set up linear programming problem.
    // We have parameter e and point x, which need to satisfy
    //  x=Ge+h;  Ce+d>=0; lb<=x<=ub, -1<=e<=+1
    // Add slack variables for the inequality constraints Ce+s=d; s>=0
    // The standard form is then
    //  x-Ge=h,  -Ce+s=d; -1<=e<=+1, lb<=x<=ub, 0<=s
    // The only dependence on the cell is in the inequality constraints for x

    // Set up linear program Ax=b; l<=x<=u.
    // Order variables as x,e,s

    // Spacial dimension nx; parameter dimension ne; number of constraints nc
    const uint nx=this->dimension();
    const uint ne=this->number_of_parameters();
    const uint nc=this->number_of_constraints();

    lp.A.resize(nx+nc,nx+ne+nc);
    lp.b.resize(nx+nc);
    lp.c.resize(nx+ne+nc);
    lp.l.resize(nx+ne+nc);
    lp.u.resize(nx+ne+nc);
    lp.vt.resize(ne+nx+nc);
    lp.p.resize(ne+nx+nc);

    // Make part of linear program only dependent on set
    // Need to set all values since matrix is uninitialised
    for(uint i=0; i!=nx; ++i) {
        for(uint j=0; j!=nx; ++j) {
            lp.A[i][j]=0;
        }
        lp.A[i][i]=+1;
        for(uint j=0; j!=ne; ++j) {
            lp.A[i][nx+j]=-this->_function[i].gradient(j);
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[i][nx+ne+j]=0;
        }
        lp.b[i]=this->_function[i].value();
    }
    for(uint i=0; i!=nc; ++i) {
        for(uint j=0; j!=nx; ++j) {
            lp.A[nx+i][j]=0;
        }
        for(uint j=0; j!=ne; ++j) {
            lp.A[nx+i][nx+j]=-this->_constraints[i].gradient(j);
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[nx+i][nx+ne+j]=0;
        }
        lp.A[nx+i][ne+nx+i]=+1;
        lp.b[nx+i]=this->_constraints[i].value();
    }
    for(uint i=0; i!=ne; ++i) {
        lp.l[nx+i]=this->_domain[i].lower();
        lp.u[nx+i]=this->_domain[i].upper();
    }
    for(uint i=0; i!=nc; ++i) {
        lp.l[nx+ne+i]=0;
        lp.u[nx+ne+i]=inf<Float>();
    }

    // Take x and s variables to be basic, so the initial basis matrix is the
    // identity
    for(uint i=0; i!=nx; ++i) {
        lp.vt[i]=BASIS;
        lp.p[i]=i;
    }
    for(uint i=0; i!=ne; ++i) {
        lp.vt[nx+i]=LOWER;
        lp.p[nx+nc+i]=nx+i;
    }
    for(uint i=0; i!=nc; ++i) {
        lp.vt[nx+ne+i]=BASIS;
        lp.p[nx+i]=nx+ne+i;
    }
    lp.B=Matrix<Float>::identity(nx+nc);

    // Make part of linear program dependent on cell be +/-infinity
    for(uint i=0; i!=nx; ++i) {
        lp.l[i]=-inf<Float>();
        lp.u[i]=+inf<Float>();
    }

}


void
AffineSet::adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    ARIADNE_ASSERT(this->dimension()==paving.dimension());
    ARIADNE_ASSERT(this->_equations.empty());

    GridCell bounding_cell=GridCell::smallest_enclosing_primary_cell(this->bounding_box(),paving.grid());

    // Create linear program
    LinearProgram<Float> lp;
    this->construct_linear_program(lp);

    _adjoin_outer_approximation_to(paving,lp,bounding_cell,depth);
}




void _robust_adjoin_outer_approximation_to(GridTreeSet& paving, LinearProgram<Float>& lp, GridCell& cell, int depth)
{
    SimplexSolver<Float> lpsolver;

    const uint nx=cell.dimension();
    const uint nc=lp.A.row_size()-nx;
    const uint ne=lp.A.column_size()-lp.A.row_size()-1u;


    // No need to check if cell is already part of the set
    if(subset(cell,paving)) {
        return;
    }

    // Find concrete cell box
    const Box& bx=cell.box();

    // Make part of linear program dependent on cell
    for(uint i=0; i!=nx; ++i) {
        lp.l[ne+i]=bx[i].lower();
        lp.u[ne+i]=bx[i].upper();
    }

    int cell_tree_depth=(cell.depth()-cell.height());
    int maximum_tree_depth=depth*cell.dimension();

    // Check for disjointness using linear program
    tribool feasible=lpsolver.constrained_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);

    bool done=false;
    while(!done && lp.x[ne+nx+nc]<0.0) {
        done=lpsolver.lpstep(lp.A,lp.b,lp.c,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x);
    }
    Matrix<Interval> B=lpsolver.compute_B<Interval>(lp.A,lp.p);
    Vector<Interval> x=lpsolver.compute_x(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,B);
    if(x[ne+nx+nc].upper()<0.0) { return; } // No feasible solution

    //feasible=verify_constrained_feasibility(lp.A,lp.b,lp.l,lp.u,lp.vt);
    //if(!feasible) { return; }

    // If not disjoint, and we are at maximum depth, adjoin the set.
    // Otherwise, split and continue.

    // FIXME: The cell depth should be an absolute depth, so this conversion should not be needed
    if(cell_tree_depth>=maximum_tree_depth) {
        paving.adjoin(cell);
    } else {
        GridCell subcell1 = cell.split(0);
        GridCell subcell2 = cell.split(1);
        _adjoin_outer_approximation_to(paving,lp,subcell1,depth);
        _adjoin_outer_approximation_to(paving,lp,subcell2,depth);
    }

}



void
AffineSet::robust_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const {
    ARIADNE_ASSERT(this->dimension()==paving.dimension());
    ARIADNE_ASSERT(this->_equations.empty());

    SimplexSolver<Float> lpsolver;

    GridCell bounding_cell=GridCell::smallest_enclosing_primary_cell(this->bounding_box(),paving.grid());

    // Set up linear programming problem.
    // We have parameter e and point x, which need to satisfy
    //  x=Ge+h;  Ce+d>=0; lb<=x<=ub, -1<=e<=+1
    // Add slack variables for the inequality constraints Ce+s=d; s>=0
    // The standard form is then
    //  -Ge+x=h,  -Ce+s=d; -1<=e<=+1, lb<=x<=ub, 0<=s
    // The only dependence on the cell is in the inequality constraints for x

    // Set up linear program Ax=b; l<=x<=u.
    // Order variables as e,x,s

    // Spacial dimension nx; parameter dimension ne; number of constraints nc
    const uint nx=this->dimension();
    const uint ne=this->number_of_parameters();
    const uint nc=this->number_of_constraints();

    // Create linear program
    LinearProgram<Float> lp;
    lp.A.resize(nx+nc,nx+ne+nc+1u);
    lp.b.resize(nx+nc);
    lp.c.resize(nx+ne+nc+1u);
    lp.l.resize(nx+ne+nc+1u);
    lp.u.resize(nx+ne+nc+1u);
    lp.vt.resize(ne+nx+nc+1u);
    lp.p.resize(ne+nx+nc+1u);

    // Make part of linear program only dependent on set
    // Need to set all values since matrix is uninitialised
    for(uint i=0; i!=nx; ++i) {
        for(uint j=0; j!=ne; ++j) {
            lp.A[i][j]=-this->_function[i].gradient(j);
        }
        for(uint j=0; j!=nx; ++j) {
            lp.A[i][ne+j]=0;
        }
        lp.A[i][ne+i]=+1;
        for(uint j=0; j!=nc; ++j) {
            lp.A[i][ne+nx+j]=0;
        }
        lp.b[i]=this->_function[i].value();
    }
    for(uint i=0; i!=nc; ++i) {
        for(uint j=0; j!=ne; ++j) {
            lp.A[nx+i][j]=-this->_constraints[i].gradient(j);
        }
        for(uint j=0; j!=nx; ++j) {
            lp.A[nx+i][ne+j]=0;
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[nx+i][ne+nx+j]=0;
        }
        lp.A[nx+i][nx+ne+i]=+1;
        lp.b[nx+i]=this->_constraints[i].value();
    }
    for(uint i=0; i!=ne; ++i) {
        lp.l[i]=-1;
        lp.u[i]=+1;
    }
    for(uint i=0; i!=nc; ++i) {
        lp.l[ne+nx+i]=0;
        lp.u[ne+nx+i]=inf<Float>();
    }


    // Make part of linear program used for robustness
    for(uint i=0; i!=ne+nx+nc; ++i) {
        lp.c[i]=0;
    }
    for(uint i=0; i!=nx+nc; ++i) {
        lp.A[i][ne+nx+nc]=1;
    }
    lp.c[ne+nx+nc]=-1;
    //lp.l[ne+nx+nc]=0;
    lp.l[ne+nx+nc]=-1;
    lp.u[ne+nx+nc]=inf<Float>();
    lp.vt[ne+nx+nc]=LOWER;
    lp.p[ne+nx+nc]=ne+nx+nc;


    // Make part of linear program dependent on cell
    const Box bx=bounding_cell.box();
    for(uint i=0; i!=nx; ++i) {
        lp.l[ne+i]=bx[i].lower();
        lp.u[ne+i]=bx[i].upper();
    }

    // Take x and s variables to be basic, so the initial basis matrix is the
    // identity
    for(uint i=0; i!=ne; ++i) {
        lp.vt[i]=LOWER;
        lp.p[nx+nc+i]=i;
    }
    for(uint i=0; i!=nx+nc; ++i) {
        lp.vt[ne+i]=BASIS;
        lp.p[i]=ne+i;
    }
    lp.B=Matrix<Float>::identity(nx+nc);

    ARIADNE_LOG(9,"A="<<lp.A<<"\nb="<<lp.b<<"\nl="<<lp.l<<"\nu="<<lp.u<<"\n");
    tribool feasible=lpsolver.constrained_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);
    ARIADNE_LOG(9,"  vt="<<lp.vt<<"\nx="<<lp.x<<"\n");
    if(!feasible) { return; } // no intersection

    _adjoin_outer_approximation_to(paving,lp,bounding_cell,depth);
}



List<Point2d>
AffineSet::boundary(uint xc, uint yc) const
{
    ARIADNE_ASSERT_MSG(this->_equations.empty(),*this);
    if(this->empty()) { return List<Point2d>(); }

    SimplexSolver<Float> lpsolver;

    static const double EQUATION_WIDENING = 1e-8;
    List< Affine<Float> > constraints=this->_constraints;
    for(uint i=0; i!=this->_equations.size(); ++i) {
        constraints.append(this->_equations[i]-EQUATION_WIDENING);
        constraints.append(-this->_equations[i]-EQUATION_WIDENING);
    }
    Box const& domain=this->domain();
    const size_t m=constraints.size();
    const size_t n=domain.size();

    Affine<Float> xa=this->_function[xc];
    Affine<Float> ya=this->_function[yc];

    // For now, assume bottom-left corner of domain is feasible and gives left-most point of
    // image

    // We always want to turn left, so we need to choose a direction d such that
    // (v0*d0+v1*d1) is positive and the ratio (v0*d1-v1*d0)/(v0*d0+v1*d1) is maximised.
    // Note that in in all cases, v0*d1-v1*d0 must be positive.

    List<Point2d> vertices;
//    vertices.append(vertex);

    // Set up linear programming problem over domain variables e
    // We have constraints l<=e<=u; Ae<=b; Introduce slack variables which start in basis
    Matrix<Float> G(2,n+m);
    for(uint j=0; j!=n; ++j) { G[0][j]=xa.gradient(j); G[1][j]=ya.gradient(j); }
    Vector<Float> h(2);
    h[0]=xa.value(); h[1]=ya.value();
    Vector<Float> pt(2);

    static const double ANTIDEGENERACY_PERTURBATION = 1e-8;
    Matrix<Float> A(m,n+m);
    Vector<Float> b(m);
    Vector<Float> l(n+m);
    Vector<Float> u(n+m);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            A[i][j]=-constraints[i].gradient(j) + ANTIDEGENERACY_PERTURBATION;
        }
        A[i][i+n]=1.0;
        b[i]=constraints[i].value()+1e-5;
    }
    for(uint j=0; j!=n; ++j) {
        l[j]=domain[j].lower() - ANTIDEGENERACY_PERTURBATION; // Add a small constant to help Simplex solver
        u[j]=domain[j].upper() - ANTIDEGENERACY_PERTURBATION;
    }
    for(uint i=0; i!=m; ++i) {
        l[n+i]=0.0 - ANTIDEGENERACY_PERTURBATION;
        u[n+i]=inf<Float>();
        //u[n+i]=+std::numeric_limits<double>::max();
    }

    Matrix<Float> B=Matrix<Float>::identity(m);

    array<Slackness> vt(n+m);
    for(uint j=0; j!=n; ++j) { vt[j]=LOWER; }
    for(uint i=0; i!=m; ++i) { vt[n+i]=BASIS; }

    array<size_t> p(m);
    for(uint i=0; i!=m; ++i) { p[i]=n+i; }

    Vector<Float> c(n+m); for(uint j=0; j!=n; ++j) { c[j]=xa[j]; }

    ARIADNE_LOG(3," A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<"\n");
    Vector<Float> e=lpsolver.optimize(A,b,c,l,u,vt,p,B);
    ARIADNE_LOG(3,"\nA="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" e="<<e<<"\n\n");

    int MAX_STEPS=10;
    int STEPS=0;
    array<Slackness> initial_variable_type=vt;

    Vector<Float> x=prod(G,e)+h;
    Vector<Float> v(2,0.0,-1.0);
    Vector<Float> v_max(2);
    Vector<Float> w(2);
    Vector<Float> Aj(m);
    Vector<Float> BAj(m);

    uint lb=m+n; // Last variable to exit the basis

    do {
        ++STEPS;
        Float cos_theta_max=-inf<Float>();
        uint s=m+n;

        // Compute direction the point Gx moves in when variable i leaves the basis
        // This direction is given by Gv where v_B=-A_B^{-1} A_N e_i
        for(uint k=m; k!=m+n; ++k) {
            // Test variable to enter basis; there are m to test, one for each dimension of the domain
            uint j=p[k];
            if(j!=lb) {
                uint j=p[k];
                Aj=column(A,j);
                BAj=prod(B,Aj);

                Vector<Float> d(m+n);
                d[j] = (vt[j]==LOWER ? +1 : -1);
                for(uint i=0; i!=m; ++i) {
                    d[p[i]]=-BAj[i]*d[j];
                }
                w=prod(G,d);

                Float dot=v[0]*w[0]+v[1]*w[1];
                Float cross=v[0]*w[1]-v[1]*w[0];
                Float cos_theta=dot/cross;
                ARIADNE_LOG(5,"  p="<<p<<" k="<<k<<" p[k]="<<p[k]<<" d="<<d<<" w="<<w<<" v="<<v<<" dot="<<dot<<" cross="<<cross<<" cos="<<cos_theta<<"\n");
                if(cross<0) { ARIADNE_WARN("AffineSet::boundary(...): cross product is="<<cross<<"; should be positive."); if(cross>-1e-8) { cross=+0.0; } }
                if(cross==0.0) { cos_theta=(dot>0.0) ? +inf<Float>() : -inf<Float>(); }
                //ARIADNE_ASSERT(cross>0);
                if(cos_theta>cos_theta_max) {
                    cos_theta_max=cos_theta;
                    v_max=w;
                    s=k;
                }
            } // k!=r
        }

        ARIADNE_ASSERT_MSG(s<m+n,"Could not find direction to move along boundary of AffineSet.");
        ARIADNE_ASSERT(vt[p[s]]!=BASIS);
        lpsolver.lpstep(A,b,l,u,vt,p,B,e,s);
        lb=p[s];
        x=prod(G,e)+h;
        v=v_max;
        ARIADNE_LOG(3,"A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" e="<<e<<" v="<<v<<" x="<<x<<"\n");

        vertices.push_back(Point2d(x[0],x[1]));

    } while(STEPS<MAX_STEPS && vt!=initial_variable_type);

    ARIADNE_LOG(3,"vertices="<<vertices<<"\n");
    return vertices;

}

void AffineSet::draw(CanvasInterface& canvas) const {
    //std::cerr<<"AffineSet::draw: "<<*this<<"\n";
    List<Point2d> boundary=this->boundary(canvas.x_coordinate(),canvas.y_coordinate());
    if(boundary.empty()) { return; }

    // Trace boundary
    Point2d prv=boundary[boundary.size()-1];
    canvas.move_to(prv.x,prv.y);
    for(uint i=0; i!=boundary.size(); ++i) {
        prv=boundary[i];
        canvas.line_to(prv.x,prv.y);
    }
    canvas.fill();
}


std::ostream& AffineSet::write(std::ostream& os) const {
    return os << "AffineSet( domain=" << this->_domain << ", function=" << this->_function << ", constraints=" << this->_constraints << ", equations=" << this->_equations <<" )";
}


} // namespace Ariadne
