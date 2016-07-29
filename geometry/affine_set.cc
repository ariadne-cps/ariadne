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

#include "function/functional.h"
#include "config.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "solvers/linear_programming.h"
#include "function/function.h"
#include "function/affine.h"
#include "function/affine_model.h"

#include "geometry/box.h"
#include "geometry/grid_set.h"
#include "geometry/affine_set.h"

#include "output/graphics_interface.h"
#include "output/geometry2d.h"


namespace Ariadne {

typedef Vector<Float64> RawFloatVector;
typedef Vector<ExactIntervalType> ExactIntervalVectorType;




template<class X>
struct LinearProgram {
    Matrix<X> A;
    Vector<X> b;
    Vector<X> c;
    Vector<X> l;
    Vector<X> u;
    Array<Slackness> vt;
    Array<SizeType> p;
    Matrix<X> B;
    Vector<X> x;
    Vector<X> y;
    Vector<X> z;
};


ValidatedAffineConstraint operator<=(const Float64Bounds& l, const ValidatedAffine& a) { return ValidatedAffineConstraint(l,a,+infty); }
ValidatedAffineConstraint operator<=(const ValidatedAffine& a, const Float64Bounds& u) { return ValidatedAffineConstraint(-infty,a,u); }
ValidatedAffineConstraint operator==(const ValidatedAffine& a, const Float64Bounds& b) { return ValidatedAffineConstraint(b,a,b); }

ValidatedAffineConstraint operator<=(const ValidatedAffineConstraint& c, const Float64Bounds& u) {
    ARIADNE_ASSERT(c.upper_bound()==infty);
    return ValidatedAffineConstraint(c.lower_bound(),c.function(),u);
}

ValidatedAffineModelConstraint operator<=(const Float64Bounds& l, const ValidatedAffineModel& a) { return ValidatedAffineModelConstraint(l,a,+infty); }
ValidatedAffineModelConstraint operator<=(const ValidatedAffineModel& a, const Float64Bounds& u) { return ValidatedAffineModelConstraint(-infty,a,u); }
ValidatedAffineModelConstraint operator==(const ValidatedAffineModel& a, const Float64Bounds& b) { return ValidatedAffineModelConstraint(b,a,b); }

ValidatedAffineModelConstraint operator<=(const ValidatedAffineModelConstraint& c, const Float64Bounds& u) {
    ARIADNE_ASSERT(c.upper_bound()==infty);
    return ValidatedAffineModelConstraint(c.lower_bound(),c.function(),u);
}

ValidatedAffineModel affine_model(const ValidatedAffine& a);

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBoxType& d,
                     const Vector<ValidatedAffine>& f)
    : _domain(d), _space_models(f.size(),ValidatedAffineModel(d.size(),Precision64()))
{
    if(d==ExactBoxType::unit_box(d.size())) {
        for(Nat i=0; i!=f.size(); ++i) {
            _space_models[i] = affine_model(f[i]);
        }
    }
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBoxType& d,
                     const Vector<ValidatedAffine>& f,
                     const List<ValidatedAffineConstraint>& c)
    : _domain(d), _space_models(f.size(),ValidatedAffineModel(d.size(),Precision64())), _constraint_models()
{
    if(d==ExactBoxType::unit_box(d.size())) {
        for(Nat i=0; i!=f.size(); ++i) {
            _space_models[i] = affine_model(f[i]);
        }
        for(Nat i=0; i!=c.size(); ++i) {
            _constraint_models.append(ValidatedAffineModelConstraint(c[i].lower_bound(),affine_model(c[i].function()),c[i].upper_bound()));
        }
    } else {
        ARIADNE_NOT_IMPLEMENTED;
    }
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBoxType& d,
                     const Vector<ValidatedAffineModel>& f,
                     const List<ValidatedAffineModelConstraint>& c)
    : _domain(d), _space_models(f), _constraint_models(c)
{
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f,
                     const List<ValidatedAffineModelConstraint>& c)
    : _domain(ExactBoxType::unit_box(f[0].argument_size())), _space_models(f), _constraint_models(c)
{
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f)
    : _domain(ExactBoxType::unit_box(f[0].argument_size())), _space_models(f)
{
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBoxType& D, const Matrix<Float64Value>& G, const Vector<Float64Value>& h)
{
    this->construct(D,G,h);
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const Matrix<Float64Value>& G, const Vector<Float64Value>& h)
{
    this->construct(Vector<ExactIntervalType>(G.column_size(),ExactIntervalType(-1,+1)),G,h);
}

Void ValidatedAffineConstrainedImageSet::construct(const ExactBoxType& D, const Matrix<Float64Value>& G, const Vector<Float64Value>& h)
{
    ARIADNE_ASSERT_MSG(G.row_size()==h.size() && G.row_size()>0,"G="<<G<<", h="<<h);
    this->_domain=D;
    this->_space_models=Vector<ValidatedAffineModel>(G.row_size(),ValidatedAffineModel(G.column_size(),Precision64()));
    for(Nat i=0; i!=G.row_size(); ++i) {
        ValidatedAffineModel x(G.column_size(),Precision64());
        x=h[i];
        for(Nat j=0; j!=G.column_size(); ++j) {
            x[j]=G[i][j];
        }
        this->_space_models[i]=x;
    }
}

Void
ValidatedAffineConstrainedImageSet::new_constraint(const ValidatedAffineModelConstraint& c)
{
    ARIADNE_ASSERT(this->_space_models.size()>0);
    ARIADNE_ASSERT_MSG(this->number_of_parameters()==c.argument_size(),"c["<<c.argument_size()<<"]="<<c<<" f["<<this->number_of_parameters()<<"]="<<this->_space_models);
    _constraint_models.append(c);
}

Void
ValidatedAffineConstrainedImageSet::new_parameter_constraint(const ValidatedAffineConstraint& c)
{
    ARIADNE_ASSERT_MSG(this->_space_models.size()>0,"f="<<this->_space_models);
    ARIADNE_ASSERT_MSG(this->number_of_parameters()==c.argument_size(),"c="<<c<<" f="<<this->_space_models);
    this->new_constraint(ValidatedAffineModelConstraint(c.lower_bound(),affine_model(c.function()),c.upper_bound()));
}





ValidatedAffineConstrainedImageSet*
ValidatedAffineConstrainedImageSet::clone() const
{
    return new ValidatedAffineConstrainedImageSet(*this);
}


DimensionType
ValidatedAffineConstrainedImageSet::dimension() const
{
    return this->_space_models.size();
}

SizeType
ValidatedAffineConstrainedImageSet::number_of_parameters() const
{
    ARIADNE_ASSERT(this->_space_models.size()>0);
    return this->_space_models[0].argument_size();
}

SizeType
ValidatedAffineConstrainedImageSet::number_of_constraints() const
{
    return this->_constraint_models.size();
}

ExactBoxType
ValidatedAffineConstrainedImageSet::domain() const
{
    return this->_domain;
}

ValidatedKleenean ValidatedAffineConstrainedImageSet::is_bounded() const {
    return ValidatedKleenean(ExactBoxType(this->domain()).is_bounded()) || ValidatedKleenean(indeterminate);
}

UpperBoxType ValidatedAffineConstrainedImageSet::bounding_box() const {
    UpperBoxType result(this->dimension());
    ExactBoxType domain=this->domain();
    for(Nat i=0; i!=this->dimension(); ++i) {
        //result[i]=evaluate(this->_space_models[i],domain);
        result[i]=this->_space_models[i].evaluate(static_cast<Vector<UpperIntervalType>>(domain));
    }
    return result;
}





ValidatedSierpinskian ValidatedAffineConstrainedImageSet::separated(const ExactBoxType& bx) const {
    ARIADNE_PRECONDITION_MSG(this->dimension()==bx.dimension(),"set="<<*this<<", box="<<bx);
    ExactBoxType wbx=cast_exact_box(widen(bx));
    LinearProgram<Float64> lp;
    this->construct_linear_program(lp);
    for(Nat i=0; i!=bx.size(); ++i) {
        lp.l[i]=sub_down(wbx[i].lower().raw(),this->_space_models[i].error().raw());
        lp.u[i]=add_up(wbx[i].upper().raw(),this->_space_models[i].error().raw());
    }
    //std::cerr<<"\ns="<<*this<<"\nbx="<<bx<<"\n\nA="<<lp.A<<"\nb="<<lp.b<<"\nl="<<lp.l<<"\nu="<<lp.u<<"\n\n";
    ValidatedKleenean feasible=indeterminate;
    try {
        InteriorPointSolver optimiser;
        feasible=optimiser.feasible(lp.l,lp.u,lp.A,lp.b);
        //feasible=SimplexSolver<Float64>().hotstarted_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);
    }
    catch(DegenerateFeasibilityProblemException e) {
        feasible=indeterminate;
    }
    return !feasible;
}

ValidatedSierpinskian ValidatedAffineConstrainedImageSet::is_empty() const {
    return ValidatedSierpinskian(this->separated(cast_exact_box(this->bounding_box()))) || ValidatedKleenean(indeterminate);
}

ValidatedSierpinskian ValidatedAffineConstrainedImageSet::inside(const ExactBoxType& bx) const {
    ARIADNE_PRECONDITION_MSG(this->dimension()==bx.dimension(),"set="<<*this<<", box="<<bx);
    return widen(this->bounding_box()).inside(bx);
}


GridTreeSet
ValidatedAffineConstrainedImageSet::outer_approximation(const Grid& g, Int d) const {
    GridTreeSet r(g);
    this->adjoin_outer_approximation_to(r,d);
    return r;
}


Void ValidatedAffineConstrainedImageSet::_adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float64>& lp, const Vector<Float64>& errors, GridCell& cell, Int depth)
{

    // No need to check if cell is already part of the set
    if(paving.superset(cell)) {
        return;
    }

    // Find concrete cell box
    const ExactBoxType bx=cell.box();

    // Make part of linear program dependent on cell
    for(Nat i=0; i!=cell.dimension(); ++i) {
        //lp.l[i]=bx[i].lower();
        //lp.u[i]=bx[i].upper();
        lp.l[i]=sub_down(bx[i].lower().raw(),errors[i].raw());
        lp.u[i]=add_up(bx[i].upper().raw(),errors[i].raw());
    }

    Int cell_tree_depth=(cell.depth()-cell.height());
    Int maximum_tree_depth=depth*cell.dimension();

    // Check for disjointness using linear program
    //ValidatedKleenean feasible=SimplexSolver<Float64>().hotstarted_feasible(lp.A,lp.b,lp.l,lp.u,lp.vt,lp.p,lp.B,lp.x,lp.y);
    InteriorPointSolver optimiser;
    ValidatedKleenean feasible=optimiser.feasible(lp.l,lp.u,lp.A,lp.b);
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
        _adjoin_outer_approximation_to(paving,lp,errors,subcell1,depth);
        _adjoin_outer_approximation_to(paving,lp,errors,subcell2,depth);
    }

}



Void
ValidatedAffineConstrainedImageSet::construct_linear_program(LinearProgram<Float64>& lp) const
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

	// Warning: Uniform part of space function is not included!

    const Nat nx=this->dimension();
    const Nat np=this->number_of_parameters();
    const Nat nc=this->_constraint_models.size();

    lp.A.resize(nx+nc,nx+np+nc);
    lp.b.resize(nx+nc);
    lp.c.resize(nx+np+nc);
    lp.l.resize(nx+np+nc);
    lp.u.resize(nx+np+nc);

    // Make part of linear program only dependent on set
    // Need to set all values since matrix is uninitialised
    for(Nat i=0; i!=nx; ++i) {
        for(Nat j=0; j!=nx; ++j) {
            lp.A[i][j]=0;
        }
        lp.A[i][i] = -1;
        for(Nat j=0; j!=np; ++j) {
            lp.A[i][nx+j] = +this->_space_models[i].gradient(j).raw();
        }
        for(Nat j=0; j!=nc; ++j) {
            lp.A[i][nx+np+j]=0;
        }
        lp.b[i] = -this->_space_models[i].value().raw();
    }
    for(Nat i=0; i!=nc; ++i) {
        for(Nat j=0; j!=nx; ++j) {
            lp.A[nx+i][j]=0;
        }
        for(Nat j=0; j!=np; ++j) {
            lp.A[nx+i][nx+j] = +this->_constraint_models[i].function().gradient(j).raw();
        }
        for(Nat j=0; j!=nc; ++j) {
            lp.A[nx+i][nx+np+j]=0;
        }
        lp.A[nx+i][nx+np+i] = +1;
        lp.b[nx+i] = -this->_constraint_models[i].function().value().raw();
    }

    for(Nat i=0; i!=np; ++i) {
        //lp.l[nx+i]=this->_domain[i].lower();
        //lp.u[nx+i]=this->_domain[i].upper();
        lp.l[nx+i]=-1.0;
        lp.u[nx+i]=+1.0;
    }
    for(Nat i=0; i!=nc; ++i) {
        lp.l[nx+np+i]=-this->_constraint_models[i].upper_bound().value().raw();
        lp.u[nx+np+i]=-this->_constraint_models[i].lower_bound().value().raw();
    }

    // Make part of linear program dependent on cell be +/-infinity
    for(Nat i=0; i!=nx; ++i) {
        lp.l[i]=-inf;
        lp.u[i]=+inf;
    }

    ARIADNE_LOG(7,"set="<<*this<<", A="<<lp.A<<", b="<<lp.b<<", l="<<lp.l<<", u="<<lp.u<<"\n");
}


Void
ValidatedAffineConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const
{
    ARIADNE_ASSERT(this->dimension()==paving.dimension());

    GridCell bounding_cell=GridCell::smallest_enclosing_primary_cell(this->bounding_box(),paving.grid());

    // Create linear program
    LinearProgram<Float64> lp;
    this->construct_linear_program(lp);
    Vector<Float64> errors(this->dimension());
    for(Nat i=0; i!=this->dimension(); ++i) {
        errors[i]=this->_space_models[i].error().raw();
    }
    _adjoin_outer_approximation_to(paving,lp,errors,bounding_cell,depth);
}




Void ValidatedAffineConstrainedImageSet::_robust_adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float64>& lp, const Vector<Float64>& errors, GridCell& cell, Int depth)
{
    SimplexSolver<Float64> lpsolver;

    const Nat nx=cell.dimension();
    const Nat nc=lp.A.row_size()-nx;
    const Nat ne=lp.A.column_size()-lp.A.row_size()-1u;


    // No need to check if cell is already part of the set
    if(paving.superset(cell)) {
        return;
    }

    // Find concrete cell box
    const ExactBoxType& bx=cell.box();

    // Make part of linear program dependent on cell
    for(Nat i=0; i!=nx; ++i) {
        lp.l[i]=bx[i].lower().raw();
        lp.u[i]=bx[i].upper().raw();
    }

    Int cell_tree_depth=(cell.depth()-cell.height());
    Int maximum_tree_depth=depth*cell.dimension();

    // Check for disjointness using linear program
    ValidatedKleenean feasible=lpsolver.hotstarted_feasible(lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x,lp.y);

    Bool done=false;
    while(!done && lp.x[ne+nx+nc]<0.0) {
        done=lpsolver.lpstep(lp.c,lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x);
    }
    // FIXME: Should be rigorous!
    Vector<Float64> x=lpsolver.compute_x(lp.l,lp.u,lp.A,lp.b,lp.vt);
    if(x[ne+nx+nc]<0.0) { return; } // No feasible solution

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
        ValidatedAffineConstrainedImageSet::_adjoin_outer_approximation_to(paving,lp,errors,subcell1,depth);
        ValidatedAffineConstrainedImageSet::_adjoin_outer_approximation_to(paving,lp,errors,subcell2,depth);
    }

}



Void
ValidatedAffineConstrainedImageSet::robust_adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const {
    ARIADNE_ASSERT(this->dimension()==paving.dimension());

    SimplexSolver<Float64> lpsolver;

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

    // Spacial dimension nx; parameter dimension np; number of constraints nc
    const Nat nx=this->dimension();
    const Nat ne=this->number_of_parameters();
    const Nat nc=this->number_of_constraints();

    // Create linear program
    LinearProgram<Float64> lp;
    lp.A.resize(nx+nc,nx+ne+nc+1u);
    lp.b.resize(nx+nc);
    lp.c.resize(nx+ne+nc+1u);
    lp.l.resize(nx+ne+nc+1u);
    lp.u.resize(nx+ne+nc+1u);
    lp.vt.resize(ne+nx+nc+1u);
    lp.p.resize(ne+nx+nc+1u);

    // Make part of linear program only dependent on set
    // Need to set all values since matrix is uninitialised
    for(Nat i=0; i!=nx; ++i) {
        for(Nat j=0; j!=ne; ++j) {
            lp.A[i][j]=-this->_space_models[i].gradient(j).raw();
        }
        for(Nat j=0; j!=nx; ++j) {
            lp.A[i][ne+j]=0;
        }
        lp.A[i][ne+i]=+1;
        for(Nat j=0; j!=nc; ++j) {
            lp.A[i][ne+nx+j]=0;
        }
        lp.b[i]=this->_space_models[i].value().raw();
    }
    for(Nat i=0; i!=nc; ++i) {
        for(Nat j=0; j!=ne; ++j) {
            lp.A[nx+i][j]=this->_constraint_models[i].function().gradient(j).raw();
        }
        for(Nat j=0; j!=nx; ++j) {
            lp.A[nx+i][ne+j]=0;
        }
        for(Nat j=0; j!=nc; ++j) {
            lp.A[nx+i][ne+nx+j]=0;
        }
        lp.A[nx+i][nx+ne+i]=+1;
        lp.b[nx+i]=-this->_constraint_models[i].function().value().raw();
    }
    for(Nat i=0; i!=ne; ++i) {
        lp.l[i]=-1;
        lp.u[i]=+1;
    }
    for(Nat i=0; i!=nc; ++i) {
        lp.l[ne+nx+i]=0;
        lp.u[ne+nx+i]=inf;
    }


    // Make part of linear program used for robustness
    for(Nat i=0; i!=ne+nx+nc; ++i) {
        lp.c[i]=0;
    }
    for(Nat i=0; i!=nx+nc; ++i) {
        lp.A[i][ne+nx+nc]=1;
    }
    lp.c[ne+nx+nc]=-1;
    //lp.l[ne+nx+nc]=0;
    lp.l[ne+nx+nc]=-1;
    lp.u[ne+nx+nc]=inf;
    lp.vt[ne+nx+nc]=LOWER;
    lp.p[ne+nx+nc]=ne+nx+nc;


    // Make part of linear program dependent on cell
    const ExactBoxType bx=bounding_cell.box();
    for(Nat i=0; i!=nx; ++i) {
        lp.l[ne+i]=bx[i].lower().raw();
        lp.u[ne+i]=bx[i].upper().raw();
    }

    // Take x and s variables to be basic, so the initial basis matrix is the
    // identity
    for(Nat i=0; i!=ne; ++i) {
        lp.vt[i]=LOWER;
        lp.p[nx+nc+i]=i;
    }
    for(Nat i=0; i!=nx+nc; ++i) {
        lp.vt[ne+i]=BASIS;
        lp.p[i]=ne+i;
    }
    lp.B=Matrix<Float64>::identity(nx+nc);

    Vector<Float64> errors(this->dimension());
    for(Nat i=0; i!=this->dimension(); ++i) {
        errors[i]=this->_space_models[i].error().raw();
    }

    ARIADNE_LOG(9,"A="<<lp.A<<"\nb="<<lp.b<<"\nl="<<lp.l<<"\nu="<<lp.u<<"\n");
    ValidatedKleenean feasible=lpsolver.hotstarted_feasible(lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x,lp.y);
    ARIADNE_LOG(9,"  vt="<<lp.vt<<"\nx="<<lp.x<<"\n");
    if(definitely(not feasible)) { return; } // no intersection

    _adjoin_outer_approximation_to(paving,lp,errors,bounding_cell,depth);
}


class PerturbationGenerator {
  public:
    PerturbationGenerator() : _n(0) { }
    double operator()() { _n=_n+3; if(_n>19) { _n=_n-37; } return _n*1e-14; }
  private:
    Int _n;
};

List<Point2d>
ValidatedAffineConstrainedImageSet::boundary(Nat xind, Nat yind) const
{
    ARIADNE_LOG(3,"ValidatedAffineConstrainedImageSet::boundary("<<xind<<","<<yind<<"): self="<<*this<<"\n");

    SimplexSolver<Float64> lpsolver;
    PerturbationGenerator eps;

    static const Int MAX_STEPS=1000;
    static const double ERROR_TOLERANCE = std::numeric_limits<float>::epsilon();
    static const Float64 inf = Ariadne::inf;

    const SizeType nx=_domain.size();
    const SizeType ne=2u;
    const SizeType nc=_constraint_models.size();
    const SizeType np=nx+ne+nc;

    // Set up matrix of function values
    ValidatedAffineModel const& xa=this->_space_models[xind];
    ValidatedAffineModel const& ya=this->_space_models[yind];

    // The set is given by (x,y)=Gs+h, where As=b and l<=s<=u
    Matrix<Float64> G=Matrix<Float64>::zero(2,np);
    for(Nat j=0; j!=nx; ++j) { G[0][j]=numeric_cast<RawFloat64>(xa.gradient(j))+eps(); G[1][j]=numeric_cast<RawFloat64>(ya.gradient(j))+eps(); }
    G[0][nx+0]=xa.error().raw()+std::abs(eps()); G[1][nx+1]=ya.error().raw()+std::abs(eps());
    Vector<Float64> h(2);
    h[0]=numeric_cast<RawFloat64>(xa.value()); h[1]=numeric_cast<RawFloat64>(ya.value());
    ARIADNE_LOG(5,"G="<<G<<" h="<<h<<"\n");

    // Set up linear programming problem Ax=b; l<=x<=u
    // Since the parameter domain is given by cl<=Ay+b+/-e<=cu, -1<=y<=+1, introduce slack variables z such that z-Ay=b, with cl-e<=z<=cu+e

    Matrix<Float64> A(nc,np);
    Vector<Float64> b(nc);
    Vector<Float64> l(np);
    Vector<Float64> u(np);

    for(Nat j=0; j!=nx; ++j) {
        l[j]=-1.0+eps();
        u[j]=+1.0+eps();
    }
    for(Nat i=0; i!=nc; ++i) {
        for(Nat j=0; j!=nx; ++j) {
            A[i][j] = neg( this->_constraint_models[i].function().gradient(j).raw() );
        }
        for(Nat j=nx; j!=nx+nc; ++j) {
            A[i][j] = 0.0;
        }
        A[i][nx+i]=1.0;
        Float64 fb=this->_constraint_models[i].function().value().raw();
        Float64 fe=this->_constraint_models[i].function().error().raw();
        Float64 cl=this->_constraint_models[i].lower_bound().value().raw();
        Float64 cu=this->_constraint_models[i].upper_bound().value().raw();
        b[i]=fb;
        l[nx+i]=cl-fe;
        u[nx+i]=cu+fe;
    }

    if(xa.error().raw()==0.0) {
        l[nx+nc]=0.0; u[nx+nc]=0.0;
    } else {
        l[nx+nc]=-1.0; u[nx+nc]=+1.0;
    }
    if(ya.error().raw()==0.0) {
        l[nx+nc+1]=0.0; u[nx+nc+1]=0.0;
    } else {
        l[nx+nc+1]=-1.0; u[nx+nc+1]=+1.0;
    }
    ARIADNE_LOG(3," A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n");

    List<Point2d> vertices;

    // Set up simplex algorithm working variables
    Array<Slackness> vt(0); Array<SizeType> p(nc); Matrix<Float64> B(nc,nc);
    Vector<Float64> x(np); Vector<Float64> y(nc);

    // Find an initial feasible point
    ValidatedKleenean feasible = lpsolver.hotstarted_feasible(l,u,A,b, vt, p,B, x,y);
    ARIADNE_LOG(3," A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<"\n  x="<<x<<" Ax="<< A*x <<"\n");
    lpsolver.consistency_check(l,u,A,b, vt,p,B,x);

    // If problem not feasible, then set is empty; return empty list
    if(!possibly(feasible)) { return vertices; }

    Vector<Float64> c(np,0.0);
    for(Nat j=0; j!=np; ++j) { c[j]=G[0][j]; }

    // Find a point on the boundary; choose the point minimising the spacial x-coordinate
    x=lpsolver.hotstarted_minimise(c,l,u,A,b, vt,p,B);
    ARIADNE_LOG(3,"\nA="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<"\n")
    ARIADNE_LOG(3,"  x="<<x<<" Ax="<<(A*x)<<" c="<<c<<" cx="<<dot(c,x)<<" pt="<<G*x+h<<"\n\n");

    lpsolver.consistency_check(l,u,A,b, vt,p,B,x);

    // We always want to turn left, so we need to choose a direction d such that
    // (v0*d0+v1*d1) is positive and the ratio (v0*d1-v1*d0)/(v0*d0+v1*d1) is maximised.
    // Note that in in all cases, v0*d1-v1*d0 must be positive.

    Int STEPS=0;
    Array<Slackness> initial_variable_type=vt;

    Vector<Float64> pt=G*x+h; // The current point in space
    Vector<Float64> last_vec({0.0,-1.0}); // The direction in space of the last step along the boundary
                                        // Should be set orthogonal to the direction (+1,0) which is minimised in finding first boundary point
    Vector<Float64> best_next_vec(2);
    Vector<Float64> trial_vec(2);
    Vector<Float64> Aj(nc); // The jth column of A
    Vector<Float64> BAj(nc); // The jth column of A

    Nat last_exiting_variable=np; // Last variable to exit the basis

    // FIXME: Do we need check below?
    // if(nx==ne) { vertices.push_back(Point2d(pt[0],pt[1])); return vertices; }

    do {
        ++STEPS;
        Float64 cot_theta_max=-inf;
        Nat s=np; // The index giving the variable x[p[s]] to enter the basis

        // Compute direction the point Gx moves in when variable x[j]=x[p[k]] enters the basis
        // This direction is given by Gd where d_B=-A_B^{-1} A_N e_j and
        //   di=+1/-1 depending on whether variable x[j] is at lower or upper bound
        for(Nat k=nc; k!=np; ++k) {
            // Test variable to enter basis; there are m to test, one for each dimension of the domain
            Nat j=p[k];
            if(j!=last_exiting_variable || true) {
                Nat j=p[k];
                Aj=column(A,j);
                BAj=B*Aj;

                // Compute the direction the point in space moves for x[j] entering the basis
                Vector<Float64> d(np,0.0); // The direction x moves in in changing the basis
                d[j] = (vt[j]==LOWER ? +1 : -1);
                for(Nat i=0; i!=nc; ++i) {
                    d[p[i]]=-BAj[i]*d[j];
                }
                trial_vec=G*d;

                // Compare the direction moved in this step with the last step.
                // dot compares the direction of the two steps;
                //   positive means an obtuse angle i.e. the same general direction.
                // cross compares whether we turn to the left or the right.
                // Since we traverse the boundary anticlockwise, cross should be positive.
                Float64 dot=last_vec[0]*trial_vec[0]+last_vec[1]*trial_vec[1];
                Float64 cross=last_vec[0]*trial_vec[1]-last_vec[1]*trial_vec[0];
                Float64 cot_theta=dot/cross;
                ARIADNE_LOG(5,"  k="<<k<<" p[k]="<<p[k]<<" d="<<d<<" trial_vec="<<trial_vec<<" last_vec="<<last_vec<<
                              " dot="<<dot<<" cross="<<cross<<" cot="<<cot_theta<<"\n");


                // Due to roundoff error, the computed cross-product may be negative.
                // If this lies within a reasonable tolerance, set to zero,
                // otherwise abort.
                ARIADNE_ASSERT_MSG(cross >= -ERROR_TOLERANCE,
                                   "ValidatedAffineConstrainedImageSet::boundary(...): cross product is="<<cross<<"; should be positive.");

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

        ARIADNE_ASSERT_MSG(s<np,"Could not find direction to move along boundary of ValidatedAffineConstrainedImageSet.");
        ARIADNE_DEBUG_ASSERT(vt[p[s]]!=BASIS);
        ARIADNE_LOG(5,"  Choosing variable x["<<p[s]<<"]=x[p["<<s<<"]] to enter basis\n");
        lpsolver.lpstep(l,u,A,b, vt,p,B,x,s);
        last_exiting_variable=p[s];
        pt=G*x+h;
        last_vec=best_next_vec;
        ARIADNE_LOG(7,"G="<<G<<" h="<<h<<" x="<<x<<" pt="<<pt<<"\n");
        ARIADNE_LOG(5,"A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" x="<<x<<" Ax="<<A*x<<" pt="<<pt<<" vec="<<best_next_vec<<"\n");

        vertices.push_back(Point2d(pt[0],pt[1]));

    } while(STEPS<MAX_STEPS && vt!=initial_variable_type);

    ARIADNE_LOG(3,"vertices="<<vertices<<"\n");
    return vertices;

}

Void ValidatedAffineConstrainedImageSet::draw(CanvasInterface& canvas, const Projection2d& projection) const {
    ARIADNE_LOG(2,"ValidatedAffineConstrainedImageSet::draw(Canvas canvas, Projection2d& projection)\n");
    ARIADNE_LOG(3,"set="<<*this<<"\n");
    ARIADNE_LOG(3,"projection="<<projection<<"\n");

    List<Point2d> boundary;

    try {
        boundary=this->boundary(projection.x_coordinate(),projection.y_coordinate());
        ARIADNE_LOG(3,"boundary="<<boundary<<"\n");
    } catch(std::runtime_error& e) {
        throw e;
    }

    if(boundary.empty()) { return; }
    if(boundary.size()==1) { canvas.dot(boundary[0].x,boundary[0].y); }

    // Trace boundary
    canvas.move_to(boundary[0].x,boundary[0].y);
    for(Nat i=1; i!=boundary.size(); ++i) {
        canvas.line_to(boundary[i].x,boundary[i].y);
    }
    canvas.line_to(boundary[0].x,boundary[0].y);
    canvas.fill();
}


OutputStream& ValidatedAffineConstrainedImageSet::write(OutputStream& os) const {
    return os << "ValidatedAffineConstrainedImageSet( domain=" << this->_domain << ", function=" << this->_space_models << ", constraints=" << this->_constraint_models <<" )";
}


} // namespace Ariadne
