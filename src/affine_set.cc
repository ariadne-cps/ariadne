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
#include "geometry2d.h"


namespace Ariadne {

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;




template<class X>
struct LinearProgram {
    Matrix<X> A;
    Vector<X> b;
    Vector<X> c;
    Vector<X> l;
    Vector<X> u;
    Array<Slackness> vt;
    Array<size_t> p;
    Matrix<X> B;
    Vector<X> x;
    Vector<X> y;
    Vector<X> z;
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
    ARIADNE_ASSERT(this->_function.size()>0);
    ARIADNE_ASSERT_MSG(a.size()==this->_function[0].argument_size(),"a="<<a<<" f="<<this->_function);
    Affine<Float> c(a.size());
    c=-b;
    for(uint j=0; j!=a.size(); ++j) {
        c[j]=a[j];
    }
    _constraints.append(c);
}

void
AffineSet::new_equality_constraint(const Vector<Float>& a, const Float& b)
{
    // Re-write the constraint ax<=b as b-ax>=0
    ARIADNE_ASSERT(this->_function.size()>0);
    ARIADNE_ASSERT_MSG(a.size()==this->_function[0].argument_size(),"a="<<a<<" f="<<this->_function);
    Affine<Float> c(a.size());
    c=-b;
    for(uint j=0; j!=a.size(); ++j) {
        c[j]=a[j];
    }
    _equations.append(c);
}


void
AffineSet::new_inequality_constraint(const Affine<Float>& a)
{
    ARIADNE_ASSERT_MSG(a.argument_size()==this->number_of_parameters(),"a="<<a<<" f="<<this->_function);
    this->_constraints.append(a);
}

void
AffineSet::new_equality_constraint(const Affine<Float>& a)
{
    ARIADNE_ASSERT_MSG(a.argument_size()==this->number_of_parameters(),"a="<<a<<" f="<<this->_function);
    this->_equations.append(a);
}


bool AffineSet::operator==(const AffineSet& other) const {
    return this->_domain==other._domain
        && this->_function==other._function
        && this->_constraints==other._constraints
        && this->_equations==other._equations;
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
    ARIADNE_ASSERT(this->_function.size()>0);
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
    ARIADNE_PRECONDITION_MSG(this->dimension()==bx.dimension(),"set="<<*this<<", box="<<bx);
    LinearProgram<Float> lp;
    this->construct_linear_program(lp);
    for(uint i=0; i!=bx.size(); ++i) {
        lp.l[i]=bx[i].lower();
        lp.u[i]=bx[i].upper();
    }
    tribool feasible=indeterminate;
    try {
        InteriorPointSolver optimiser;
        feasible=optimiser.feasible(lp.l,lp.u,lp.A,lp.b);
        //feasible=SimplexSolver<Float>().hotstarted_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);
    }
    catch(DegenerateFeasibilityProblemException e) {
        feasible=indeterminate;
    }
    return !feasible;
}

tribool AffineSet::empty() const {
    return this->disjoint(this->bounding_box());
}


GridTreeSet
AffineSet::outer_approximation(const Grid& g, int d) const {
    GridTreeSet r(g);
    this->adjoin_outer_approximation_to(r,d);
    return r;
}


void AffineSet::_adjoin_outer_approximation_to(GridTreeSet& paving, LinearProgram<Float>& lp, GridCell& cell, int depth)
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
    //tribool feasible=SimplexSolver<Float>().hotstarted_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);
    InteriorPointSolver optimiser;
    tribool feasible=optimiser.feasible(lp.l,lp.u,lp.A,lp.b);
    //feasible=verify_feasibility(lp.A,lp.b,lp.l,lp.u,lp.vt);
    if(definitely(!feasible)) { return; }

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
    //  x=Gp+h;  Cp+d<=0; Ep+f = 0; lb<=x<=ub, -1<=p<=+1
    // Add slack variables for the inequality constraints Cp+s=d; s>=0
    // The standard form is then
    //  -x+Gp=-h,  Cp+s=-d; Ep = -f; lb<=x<=ub, -1<=p<=+1, 0<=s
    // The only dependence on the cell is in the inequality constraints for x

    // Set up linear program of form Ax=b; l<=x<=u.
    //  Order variables as x,e,s
    //  A=( -I G 0 )  b=( -h )  l=(xl -1 0) u=(xu +1 inf)
    //    (  0 C I )    ( -d )
    //    (  0 E 0 )    ( -f )
    // Spacial dimension nx; parameter dimension ne; number of constraints nc
    const uint nx=this->dimension();
    const uint np=this->number_of_parameters();
    const uint nc=this->_constraints.size();
    const uint ne=this->_equations.size();

    lp.A.resize(nx+nc+ne,nx+np+nc);
    lp.b.resize(nx+nc+ne);
    lp.c.resize(nx+np+nc);
    lp.l.resize(nx+np+nc);
    lp.u.resize(nx+np+nc);

    // Make part of linear program only dependent on set
    // Need to set all values since matrix is uninitialised
    for(uint i=0; i!=nx; ++i) {
        for(uint j=0; j!=nx; ++j) {
            lp.A[i][j]=0;
        }
        lp.A[i][i] = -1;
        for(uint j=0; j!=np; ++j) {
            lp.A[i][nx+j] = +this->_function[i].gradient(j);
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[i][nx+np+j]=0;
        }
        lp.b[i] = -this->_function[i].value();
    }
    for(uint i=0; i!=nc; ++i) {
        for(uint j=0; j!=nx; ++j) {
            lp.A[nx+i][j]=0;
        }
        for(uint j=0; j!=np; ++j) {
            lp.A[nx+i][nx+j] = +this->_constraints[i].gradient(j);
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[nx+i][nx+np+j]=0;
        }
        lp.A[nx+i][nx+np+i] = +1;
        lp.b[nx+i] = -this->_constraints[i].value();
    }
    for(uint i=0; i!=ne; ++i) {
        for(uint j=0; j!=nx; ++j) {
            lp.A[nx+nc+i][j]=0;
        }
        for(uint j=0; j!=np; ++j) {
            lp.A[nx+nc+i][nx+j] = +this->_equations[i].gradient(j);
        }
        for(uint j=0; j!=ne; ++j) {
            lp.A[nx+nc+i][nx+np+j]=0;
        }
        lp.b[nx+nc+i] = -this->_equations[i].value();
    }

    for(uint i=0; i!=np; ++i) {
        lp.l[nx+i]=this->_domain[i].lower();
        lp.u[nx+i]=this->_domain[i].upper();
    }
    for(uint i=0; i!=nc; ++i) {
        lp.l[nx+np+i]=0;
        lp.u[nx+np+i]=inf<Float>();
    }

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




void AffineSet::_robust_adjoin_outer_approximation_to(GridTreeSet& paving, LinearProgram<Float>& lp, GridCell& cell, int depth)
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
        lp.l[i]=bx[i].lower();
        lp.u[i]=bx[i].upper();
    }

    int cell_tree_depth=(cell.depth()-cell.height());
    int maximum_tree_depth=depth*cell.dimension();

    // Check for disjointness using linear program
    tribool feasible=lpsolver.hotstarted_feasible(lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x,lp.y);

    bool done=false;
    while(!done && lp.x[ne+nx+nc]<0.0) {
        done=lpsolver.lpstep(lp.c,lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x);
    }
    Vector<Interval> x=lpsolver.compute_x(lp.l,lp.u,lp.A,lp.b,lp.vt);
    if(x[ne+nx+nc].upper()<0.0) { return; } // No feasible solution

    //feasible=verify_feasibility(lp.A,lp.b,lp.l,lp.u,lp.vt);
    //if(!feasible) { return; }

    // If not disjoint, and we are at maximum depth, adjoin the set.
    // Otherwise, split and continue.

    // FIXME: The cell depth should be an absolute depth, so this conversion should not be needed
    if(cell_tree_depth>=maximum_tree_depth) {
        paving.adjoin(cell);
    } else {
        GridCell subcell1 = cell.split(0);
        GridCell subcell2 = cell.split(1);
        AffineSet::_adjoin_outer_approximation_to(paving,lp,subcell1,depth);
        AffineSet::_adjoin_outer_approximation_to(paving,lp,subcell2,depth);
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
            lp.A[nx+i][j]=this->_constraints[i].gradient(j);
        }
        for(uint j=0; j!=nx; ++j) {
            lp.A[nx+i][ne+j]=0;
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[nx+i][ne+nx+j]=0;
        }
        lp.A[nx+i][nx+ne+i]=+1;
        lp.b[nx+i]=-this->_constraints[i].value();
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
    tribool feasible=lpsolver.hotstarted_feasible(lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x,lp.y);
    ARIADNE_LOG(9,"  vt="<<lp.vt<<"\nx="<<lp.x<<"\n");
    if(!feasible) { return; } // no intersection

    _adjoin_outer_approximation_to(paving,lp,bounding_cell,depth);
}

List<Point2d>
AffineSet::boundary(uint xind, uint yind) const
{
    //ARIADNE_LOG(5,"AffineSet::boundary(xi,yi): self="<<*this<<"\n");

    SimplexSolver<Float> lpsolver;

    static const int MAX_STEPS=1000;
    static const double ERROR_TOLERANCE = std::numeric_limits<float>::epsilon();
    static const Float inf = Ariadne::inf<Float>();

    Box const& domain=this->_domain;
    List< Affine<Float> > const& function=this->_function;
    List< Affine<Float> > const& negative_constraints=this->_constraints;
    List< Affine<Float> > const& zero_constraints=this->_equations;

    const size_t nx=domain.size();
    const size_t nc=negative_constraints.size();
    const size_t ne=zero_constraints.size();

    Affine<Float> xa=function[xind];
    Affine<Float> ya=function[yind];

    List<Point2d> vertices;

    // Set up matrix of function values
    Matrix<Float> G=Matrix<Float>::zero(2,nx+nc);
    for(uint j=0; j!=nx; ++j) { G[0][j]=numeric_cast<float>(xa.gradient(j)); G[1][j]=numeric_cast<float>(ya.gradient(j)); }
    Vector<Float> h(2);
    h[0]=numeric_cast<float>(xa.value()); h[1]=numeric_cast<float>(ya.value());
    ARIADNE_LOG(5,"G="<<G<<" h="<<h<<"\n");

    // Set up linear programming problem over domain variables (x,w)
    // We have constraints l<=(x,w)<=u; Ac*x-w=0; Ae*s=be

    Matrix<Float> A(nc+ne,nx+nc);
    Vector<Float> b(nc+ne);
    Vector<Float> l(nx+nc);
    Vector<Float> u(nx+nc);
    for(uint i=0; i!=nc; ++i) {
        for(uint j=0; j!=nx; ++j) {
            A[i][j]=numeric_cast<float>(negative_constraints[i].gradient(j));
        }
        A[i][i+nx]=1.0;
        b[i]=numeric_cast<float>(-negative_constraints[i].value());
    }
    for(uint i=0; i!=ne; ++i) {
        for(uint j=0; j!=nx; ++j) {
            A[i+nc][j]=numeric_cast<float>(zero_constraints[i].gradient(j));
        }
        b[i+nc]=numeric_cast<float>(-zero_constraints[i].value());
    }
    for(uint j=0; j!=nx; ++j) {
        l[j]=numeric_cast<float>(domain[j].lower());
        u[j]=numeric_cast<float>(domain[j].upper());
    }
    for(uint i=0; i!=nc; ++i) {
        l[nx+i]=0.0;
        u[nx+i]=inf;
        //u[n+i]=+std::numeric_limits<double>::max();
    }
    ARIADNE_LOG(3," A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n");

    // Set up simplex algorithm working variables
    Array<Slackness> vt(0); Array<size_t> p(nc+ne); Matrix<Float> B(nx+nc,nx+nc);
    Vector<Float> x(nx+nc); Vector<Float> y(nc+ne);

    // Find an initial feasible point
    tribool feasible = lpsolver.hotstarted_feasible(l,u,A,b, vt, p,B, x,y);
    ARIADNE_LOG(3," A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<"\n  x="<<x<<" Ax="<< A*x <<"\n");
    lpsolver.consistency_check(l,u,A,b, vt,p,B,x);

    // If problem not feasible, then set is empty; return empty list
    if(!possibly(feasible)) { return vertices; }

    Vector<Float> c(nx+nc,0.0);
    for(uint j=0; j!=nx; ++j) { c[j]=xa[j]; }

    // Find a point on the boundary; choost the point minimising the spacial x-coordinate
    x=lpsolver.hotstarted_minimise(c,l,u,A,b, vt,p,B);
    ARIADNE_LOG(3,"\nA="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<"\n")
    ARIADNE_LOG(3,"  x="<<x<<" Ax="<<(A*x)<<" c="<<c<<" cx="<<dot(c,x)<<" pt="<<G*x+h<<"\n\n");

    lpsolver.consistency_check(l,u,A,b, vt,p,B,x);

    // We always want to turn left, so we need to choose a direction d such that
    // (v0*d0+v1*d1) is positive and the ratio (v0*d1-v1*d0)/(v0*d0+v1*d1) is maximised.
    // Note that in in all cases, v0*d1-v1*d0 must be positive.

    int STEPS=0;
    Array<Slackness> initial_variable_type=vt;

    Vector<Float> pt=G*x+h; // The current point in space
    Vector<Float> last_vec(2,0.0,-1.0); // The direction in space of the last step along the boundary
                                        // Should be set orthogonal to the direction (+1,0) which is minimised in finding first boundary point
    Vector<Float> best_next_vec(2);
    Vector<Float> trial_vec(2);
    Vector<Float> Aj(nc+ne); // The jth column of A
    Vector<Float> BAj(nc+ne); // The jth column of A

    uint last_exiting_variable=nc+nx; // Last variable to exit the basis

    if(nx==ne) { vertices.push_back(Point2d(pt[0],pt[1])); return vertices; }

    do {
        ++STEPS;
        Float cot_theta_max=-inf;
        uint s=nx+nc; // The index giving the variable x[p[s]] to enter the basis

        // Compute direction the point Gx moves in when variable x[j]=x[p[k]] enters the basis
        // This direction is given by Gd where d_B=-A_B^{-1} A_N e_j and
        //   di=+1/-1 depending on whether variable x[j] is at lower or upper bound
        for(uint k=nc+ne; k!=nx+nc; ++k) {
            // Test variable to enter basis; there are m to test, one for each dimension of the domain
            uint j=p[k];
            if(j!=last_exiting_variable || true) {
                uint j=p[k];
                Aj=column(A,j);
                BAj=B*Aj;

                // Compute the direction the point in space moves for x[j] entering the basis
                Vector<Float> d(nx+nc,0.0); // The direction x moves in in changing the basis
                d[j] = (vt[j]==LOWER ? +1 : -1);
                for(uint i=0; i!=nc+ne; ++i) {
                    d[p[i]]=-BAj[i]*d[j];
                }
                trial_vec=G*d;

                // Compare the direction moved in this step with the last step.
                // dot compares the direction of the two steps;
                //   positive means an obtuse angle i.e. the same general direction.
                // cross compares whether we turn to the left or the right.
                // Since we traverse the boundary anticlockwise, cross should be positive.
                Float dot=last_vec[0]*trial_vec[0]+last_vec[1]*trial_vec[1];
                Float cross=last_vec[0]*trial_vec[1]-last_vec[1]*trial_vec[0];
                Float cot_theta=dot/cross;
                ARIADNE_LOG(5,"  k="<<k<<" p[k]="<<p[k]<<" d="<<d<<" trial_vec="<<trial_vec<<" last_vec="<<last_vec<<
                              " dot="<<dot<<" cross="<<cross<<" cot="<<cot_theta<<"\n");


                // Due to roundoff error, the computed cross-product may be negative.
                // If this lies within a reasonable tolerance, set to zero,
                // otherwise abort.
                ARIADNE_ASSERT_MSG(cross >= -ERROR_TOLERANCE,
                                   "AffineSet::boundary(...): cross product is="<<cross<<"; should be positive.");

                if(cross<=0.0 ) {
                    cross=+0.0;
                    cot_theta=(dot>0.0) ? +inf : -inf;
                }

                // Allow for equality; in particular, if cot_theta=-infty,
                // then may turn round and go backwards.
                if(cot_theta>=cot_theta_max) {
                    cot_theta_max=cot_theta;
                    best_next_vec=trial_vec;
                    s=k;
                }
            } // k!=r
        }

        ARIADNE_ASSERT_MSG(s<nc+nx,"Could not find direction to move along boundary of AffineSet.");
        ARIADNE_DEBUG_ASSERT(vt[p[s]]!=BASIS);
        ARIADNE_LOG(3,"  Choosing variable x["<<p[s]<<"]=x[p["<<s<<"]] to enter basis\n");
        lpsolver.lpstep(l,u,A,b, vt,p,B,x,s);
        last_exiting_variable=p[s];
        pt=G*x+h;
        last_vec=best_next_vec;
        ARIADNE_LOG(5,"G="<<G<<" h="<<h<<" x="<<x<<" pt="<<pt<<"\n");
        ARIADNE_LOG(3,"A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" x="<<x<<" Ax="<<A*x<<" pt="<<pt<<" vec="<<best_next_vec<<"\n");

        vertices.push_back(Point2d(pt[0],pt[1]));

    } while(STEPS<MAX_STEPS && vt!=initial_variable_type);

    ARIADNE_LOG(3,"vertices="<<vertices<<"\n");
    return vertices;

}

void AffineSet::draw(CanvasInterface& canvas, const Projection2d& projection) const {
    //std::cerr<<"AffineSet::draw: "<<*this<<"\n";
    List<Point2d> boundary=this->boundary(projection.x_coordinate(),projection.y_coordinate());

    if(boundary.empty()) { return; }
    if(boundary.size()==1) { canvas.dot(boundary[0].x,boundary[0].y); }

    // Trace boundary
    canvas.move_to(boundary[0].x,boundary[0].y);
    for(uint i=1; i!=boundary.size(); ++i) {
        canvas.line_to(boundary[i].x,boundary[i].y);
    }
    canvas.line_to(boundary[0].x,boundary[0].y);
    canvas.fill();
}


std::ostream& AffineSet::write(std::ostream& os) const {
    return os << "AffineSet( domain=" << this->_domain << ", function=" << this->_function << ", negative_constraints=" << this->_constraints << ", zero_constraints=" << this->_equations <<" )";
}


} // namespace Ariadne
