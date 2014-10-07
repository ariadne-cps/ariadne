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

#include "functional.h"
#include "config.h"

#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "linear_programming.h"
#include "function.h"
#include "affine.h"
#include "affine_model.h"

#include "box.h"
#include "grid_set.h"
#include "affine_set.h"

#include "graphics_interface.h"
#include "geometry2d.h"


namespace Ariadne {

typedef Vector<Float> RawFloatVector;
typedef Vector<ExactInterval> ExactIntervalVector;




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


ValidatedAffineConstraint operator<=(const ValidatedFloat& l, const ValidatedAffine& a) { return ValidatedAffineConstraint(l,a,+inf); }
ValidatedAffineConstraint operator<=(const ValidatedAffine& a, const ValidatedFloat& u) { return ValidatedAffineConstraint(-inf,a,u); }
ValidatedAffineConstraint operator==(const ValidatedAffine& a, const ValidatedFloat& b) { return ValidatedAffineConstraint(b,a,b); }

ValidatedAffineConstraint operator<=(const ValidatedAffineConstraint& c, const ValidatedFloat& u) {
    ARIADNE_ASSERT(c.upper_bound()==infty);
    return ValidatedAffineConstraint(c.lower_bound(),c.function(),u);
}

ValidatedAffineModelConstraint operator<=(const ValidatedFloat& l, const ValidatedAffineModel& a) { return ValidatedAffineModelConstraint(l,a,+inf); }
ValidatedAffineModelConstraint operator<=(const ValidatedAffineModel& a, const ValidatedFloat& u) { return ValidatedAffineModelConstraint(-inf,a,u); }
ValidatedAffineModelConstraint operator==(const ValidatedAffineModel& a, const ValidatedFloat& b) { return ValidatedAffineModelConstraint(b,a,b); }

ValidatedAffineModelConstraint operator<=(const ValidatedAffineModelConstraint& c, const ValidatedFloat& u) {
    ARIADNE_ASSERT(c.upper_bound()==infty);
    return ValidatedAffineModelConstraint(c.lower_bound(),c.function(),u);
}

ValidatedAffineModel affine_model(const ValidatedAffine& a);

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBox& d,
                     const Vector<ValidatedAffine>& f)
    : _domain(d), _space_models(f.size(),ValidatedAffineModel(d.size()))
{
    if(d==ExactBox::unit_box(d.size())) {
        for(uint i=0; i!=f.size(); ++i) {
            _space_models[i] = affine_model(f[i]);
        }
    }
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBox& d,
                     const Vector<ValidatedAffine>& f,
                     const List<ValidatedAffineConstraint>& c)
    : _domain(d), _space_models(f.size(),ValidatedAffineModel(d.size())), _constraint_models()
{
    if(d==ExactBox::unit_box(d.size())) {
        for(uint i=0; i!=f.size(); ++i) {
            _space_models[i] = affine_model(f[i]);
        }
        for(uint i=0; i!=c.size(); ++i) {
            _constraint_models.append(ValidatedAffineModelConstraint(c[i].lower_bound(),affine_model(c[i].function()),c[i].upper_bound()));
        }
    } else {
        ARIADNE_NOT_IMPLEMENTED;
    }
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBox& d,
                     const Vector<ValidatedAffineModel>& f,
                     const List<ValidatedAffineModelConstraint>& c)
    : _domain(d), _space_models(f), _constraint_models(c)
{
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f,
                     const List<ValidatedAffineModelConstraint>& c)
    : _domain(ExactBox::unit_box(f[0].argument_size())), _space_models(f), _constraint_models(c)
{
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const Vector<ValidatedAffineModel>& f)
    : _domain(ExactBox::unit_box(f[0].argument_size())), _space_models(f)
{
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const ExactBox& D, const Matrix<ExactFloat>& G, const Vector<ExactFloat>& h)
{
    this->construct(D,G,h);
}

ValidatedAffineConstrainedImageSet::ValidatedAffineConstrainedImageSet(const Matrix<ExactFloat>& G, const Vector<ExactFloat>& h)
{
    this->construct(Vector<ExactInterval>(G.column_size(),ExactInterval(-1,+1)),G,h);
}

void ValidatedAffineConstrainedImageSet::construct(const ExactBox& D, const Matrix<ExactFloat>& G, const Vector<ExactFloat>& h)
{
    ARIADNE_ASSERT_MSG(G.row_size()==h.size() && G.row_size()>0,"G="<<G<<", h="<<h);
    this->_domain=D;
    this->_space_models=Vector<ValidatedAffineModel>(G.row_size(),ValidatedAffineModel(G.column_size()));
    for(uint i=0; i!=G.row_size(); ++i) {
        AffineModel<ValidatedTag> x(G.column_size());
        x=h[i];
        for(uint j=0; j!=G.column_size(); ++j) {
            x[j]=G[i][j];
        }
        this->_space_models[i]=x;
    }
}

void
ValidatedAffineConstrainedImageSet::new_constraint(const ValidatedAffineModelConstraint& c)
{
    ARIADNE_ASSERT(this->_space_models.size()>0);
    ARIADNE_ASSERT_MSG(this->number_of_parameters()==c.argument_size(),"c["<<c.argument_size()<<"]="<<c<<" f["<<this->number_of_parameters()<<"]="<<this->_space_models);
    _constraint_models.append(c);
}

void
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


uint
ValidatedAffineConstrainedImageSet::dimension() const
{
    return this->_space_models.size();
}

uint
ValidatedAffineConstrainedImageSet::number_of_parameters() const
{
    ARIADNE_ASSERT(this->_space_models.size()>0);
    return this->_space_models[0].argument_size();
}

uint
ValidatedAffineConstrainedImageSet::number_of_constraints() const
{
    return this->_constraint_models.size();
}

ExactBox
ValidatedAffineConstrainedImageSet::domain() const
{
    return this->_domain;
}

tribool ValidatedAffineConstrainedImageSet::bounded() const {
    return ExactBox(this->domain()).bounded() || indeterminate;
}

UpperBox ValidatedAffineConstrainedImageSet::bounding_box() const {
    UpperBox result(this->dimension());
    ExactBox domain=this->domain();
    for(uint i=0; i!=this->dimension(); ++i) {
        //result[i]=evaluate(this->_space_models[i],domain);
        result[i]=this->_space_models[i].evaluate(static_cast<Vector<UpperInterval>>(domain));
    }
    return result;
}





tribool ValidatedAffineConstrainedImageSet::separated(const ExactBox& bx) const {
    ARIADNE_PRECONDITION_MSG(this->dimension()==bx.dimension(),"set="<<*this<<", box="<<bx);
    ExactBox wbx=widen(bx);
    LinearProgram<Float> lp;
    this->construct_linear_program(lp);
    for(uint i=0; i!=bx.size(); ++i) {
        lp.l[i]=sub_down(wbx[i].lower().raw(),this->_space_models[i].error().raw());
        lp.u[i]=add_up(wbx[i].upper().raw(),this->_space_models[i].error().raw());
    }
    //std::cerr<<"\ns="<<*this<<"\nbx="<<bx<<"\n\nA="<<lp.A<<"\nb="<<lp.b<<"\nl="<<lp.l<<"\nu="<<lp.u<<"\n\n";
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

tribool ValidatedAffineConstrainedImageSet::empty() const {
    return this->separated(make_exact_box(this->bounding_box()));
}

tribool ValidatedAffineConstrainedImageSet::inside(const ExactBox& bx) const {
    ARIADNE_PRECONDITION_MSG(this->dimension()==bx.dimension(),"set="<<*this<<", box="<<bx);
    return widen(this->bounding_box()).inside(bx) || indeterminate;
}


GridTreeSet
ValidatedAffineConstrainedImageSet::outer_approximation(const Grid& g, int d) const {
    GridTreeSet r(g);
    this->adjoin_outer_approximation_to(r,d);
    return r;
}


void ValidatedAffineConstrainedImageSet::_adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float>& lp, const Vector<Float>& errors, GridCell& cell, int depth)
{

    // No need to check if cell is already part of the set
    if(paving.superset(cell)) {
        return;
    }

    // Find concrete cell box
    const ExactBox bx=cell.box();

    // Make part of linear program dependent on cell
    for(uint i=0; i!=cell.dimension(); ++i) {
        //lp.l[i]=bx[i].lower();
        //lp.u[i]=bx[i].upper();
        lp.l[i]=sub_down(bx[i].lower().raw(),errors[i].raw());
        lp.u[i]=add_up(bx[i].upper().raw(),errors[i].raw());
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
        _adjoin_outer_approximation_to(paving,lp,errors,subcell1,depth);
        _adjoin_outer_approximation_to(paving,lp,errors,subcell2,depth);
    }

}



void
ValidatedAffineConstrainedImageSet::construct_linear_program(LinearProgram<Float>& lp) const
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

    const uint nx=this->dimension();
    const uint np=this->number_of_parameters();
    const uint nc=this->_constraint_models.size();

    lp.A.resize(nx+nc,nx+np+nc);
    lp.b.resize(nx+nc);
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
            lp.A[i][nx+j] = +this->_space_models[i].gradient(j).raw();
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[i][nx+np+j]=0;
        }
        lp.b[i] = -this->_space_models[i].value().raw();
    }
    for(uint i=0; i!=nc; ++i) {
        for(uint j=0; j!=nx; ++j) {
            lp.A[nx+i][j]=0;
        }
        for(uint j=0; j!=np; ++j) {
            lp.A[nx+i][nx+j] = +this->_constraint_models[i].function().gradient(j).raw();
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[nx+i][nx+np+j]=0;
        }
        lp.A[nx+i][nx+np+i] = +1;
        lp.b[nx+i] = -this->_constraint_models[i].function().value().raw();
    }

    for(uint i=0; i!=np; ++i) {
        //lp.l[nx+i]=this->_domain[i].lower();
        //lp.u[nx+i]=this->_domain[i].upper();
        lp.l[nx+i]=-1.0;
        lp.u[nx+i]=+1.0;
    }
    for(uint i=0; i!=nc; ++i) {
        lp.l[nx+np+i]=-this->_constraint_models[i].upper_bound().midpoint().raw();
        lp.u[nx+np+i]=-this->_constraint_models[i].lower_bound().midpoint().raw();
    }

    // Make part of linear program dependent on cell be +/-infinity
    for(uint i=0; i!=nx; ++i) {
        lp.l[i]=-inf;
        lp.u[i]=+inf;
    }

    ARIADNE_LOG(7,"set="<<*this<<", A="<<lp.A<<", b="<<lp.b<<", l="<<lp.l<<", u="<<lp.u<<"\n");
}


void
ValidatedAffineConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    ARIADNE_ASSERT(this->dimension()==paving.dimension());

    GridCell bounding_cell=GridCell::smallest_enclosing_primary_cell(this->bounding_box(),paving.grid());

    // Create linear program
    LinearProgram<Float> lp;
    this->construct_linear_program(lp);
    Vector<Float> errors(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        errors[i]=this->_space_models[i].error().raw();
    }
    _adjoin_outer_approximation_to(paving,lp,errors,bounding_cell,depth);
}




void ValidatedAffineConstrainedImageSet::_robust_adjoin_outer_approximation_to(PavingInterface& paving, LinearProgram<Float>& lp, const Vector<Float>& errors, GridCell& cell, int depth)
{
    SimplexSolver<Float> lpsolver;

    const uint nx=cell.dimension();
    const uint nc=lp.A.row_size()-nx;
    const uint ne=lp.A.column_size()-lp.A.row_size()-1u;


    // No need to check if cell is already part of the set
    if(paving.superset(cell)) {
        return;
    }

    // Find concrete cell box
    const ExactBox& bx=cell.box();

    // Make part of linear program dependent on cell
    for(uint i=0; i!=nx; ++i) {
        lp.l[i]=bx[i].lower().raw();
        lp.u[i]=bx[i].upper().raw();
    }

    int cell_tree_depth=(cell.depth()-cell.height());
    int maximum_tree_depth=depth*cell.dimension();

    // Check for disjointness using linear program
    tribool feasible=lpsolver.hotstarted_feasible(lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x,lp.y);

    bool done=false;
    while(!done && lp.x[ne+nx+nc]<0.0) {
        done=lpsolver.lpstep(lp.c,lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x);
    }
    Vector<UpperInterval> x=lpsolver.compute_x(lp.l,lp.u,lp.A,lp.b,lp.vt);
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
        ValidatedAffineConstrainedImageSet::_adjoin_outer_approximation_to(paving,lp,errors,subcell1,depth);
        ValidatedAffineConstrainedImageSet::_adjoin_outer_approximation_to(paving,lp,errors,subcell2,depth);
    }

}



void
ValidatedAffineConstrainedImageSet::robust_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const {
    ARIADNE_ASSERT(this->dimension()==paving.dimension());

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

    // Spacial dimension nx; parameter dimension np; number of constraints nc
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
            lp.A[i][j]=-this->_space_models[i].gradient(j).raw();
        }
        for(uint j=0; j!=nx; ++j) {
            lp.A[i][ne+j]=0;
        }
        lp.A[i][ne+i]=+1;
        for(uint j=0; j!=nc; ++j) {
            lp.A[i][ne+nx+j]=0;
        }
        lp.b[i]=this->_space_models[i].value().raw();
    }
    for(uint i=0; i!=nc; ++i) {
        for(uint j=0; j!=ne; ++j) {
            lp.A[nx+i][j]=this->_constraint_models[i].function().gradient(j).raw();
        }
        for(uint j=0; j!=nx; ++j) {
            lp.A[nx+i][ne+j]=0;
        }
        for(uint j=0; j!=nc; ++j) {
            lp.A[nx+i][ne+nx+j]=0;
        }
        lp.A[nx+i][nx+ne+i]=+1;
        lp.b[nx+i]=-this->_constraint_models[i].function().value().raw();
    }
    for(uint i=0; i!=ne; ++i) {
        lp.l[i]=-1;
        lp.u[i]=+1;
    }
    for(uint i=0; i!=nc; ++i) {
        lp.l[ne+nx+i]=0;
        lp.u[ne+nx+i]=inf;
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
    lp.u[ne+nx+nc]=inf;
    lp.vt[ne+nx+nc]=LOWER;
    lp.p[ne+nx+nc]=ne+nx+nc;


    // Make part of linear program dependent on cell
    const ExactBox bx=bounding_cell.box();
    for(uint i=0; i!=nx; ++i) {
        lp.l[ne+i]=bx[i].lower().raw();
        lp.u[ne+i]=bx[i].upper().raw();
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

    Vector<Float> errors(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        errors[i]=this->_space_models[i].error().raw();
    }

    ARIADNE_LOG(9,"A="<<lp.A<<"\nb="<<lp.b<<"\nl="<<lp.l<<"\nu="<<lp.u<<"\n");
    tribool feasible=lpsolver.hotstarted_feasible(lp.l,lp.u,lp.A,lp.b,lp.vt,lp.p,lp.B,lp.x,lp.y);
    ARIADNE_LOG(9,"  vt="<<lp.vt<<"\nx="<<lp.x<<"\n");
    if(!feasible) { return; } // no intersection

    _adjoin_outer_approximation_to(paving,lp,errors,bounding_cell,depth);
}


class PerturbationGenerator {
  public:
    PerturbationGenerator() : _n(0) { }
    double operator()() { _n=_n+3; if(_n>19) { _n=_n-37; } return _n*1e-14; }
  private:
    int _n;
};

List<Point2d>
ValidatedAffineConstrainedImageSet::boundary(uint xind, uint yind) const
{
    ARIADNE_LOG(3,"ValidatedAffineConstrainedImageSet::boundary("<<xind<<","<<yind<<"): self="<<*this<<"\n");

    SimplexSolver<Float> lpsolver;
    PerturbationGenerator eps;

    static const int MAX_STEPS=1000;
    static const double ERROR_TOLERANCE = std::numeric_limits<float>::epsilon();
    static const Float inf = Ariadne::inf;

    const size_t nx=_domain.size();
    const size_t ne=2u;
    const size_t nc=_constraint_models.size();
    const size_t np=nx+ne+nc;

    // Set up matrix of function values
    AffineModel<ValidatedTag> const& xa=this->_space_models[xind];
    AffineModel<ValidatedTag> const& ya=this->_space_models[yind];

    // The set is given by (x,y)=Gs+h, where As=b and l<=s<=u
    Matrix<Float> G=Matrix<Float>::zero(2,np);
    for(uint j=0; j!=nx; ++j) { G[0][j]=numeric_cast<RawFloat>(xa.gradient(j))+eps(); G[1][j]=numeric_cast<RawFloat>(ya.gradient(j))+eps(); }
    G[0][nx+0]=xa.error().raw()+std::abs(eps()); G[1][nx+1]=ya.error().raw()+std::abs(eps());
    Vector<Float> h(2);
    h[0]=numeric_cast<RawFloat>(xa.value()); h[1]=numeric_cast<RawFloat>(ya.value());
    ARIADNE_LOG(5,"G="<<G<<" h="<<h<<"\n");

    // Set up linear programming problem Ax=b; l<=x<=u
    // Since the parameter domain is given by cl<=Ay+b+/-e<=cu, -1<=y<=+1, introduce slack variables z such that z-Ay=b, with cl-e<=z<=cu+e

    Matrix<Float> A(nc,np);
    Vector<Float> b(nc);
    Vector<Float> l(np);
    Vector<Float> u(np);

    for(uint j=0; j!=nx; ++j) {
        l[j]=-1.0+eps();
        u[j]=+1.0+eps();
    }
    for(uint i=0; i!=nc; ++i) {
        for(uint j=0; j!=nx; ++j) {
            A[i][j] = neg( this->_constraint_models[i].function().gradient(j).raw() );
        }
        for(uint j=nx; j!=nx+nc; ++j) {
            A[i][j] = 0.0;
        }
        A[i][nx+i]=1.0;
        Float fb=this->_constraint_models[i].function().value().raw();
        Float fe=this->_constraint_models[i].function().error().raw();
        Float cl=this->_constraint_models[i].lower_bound().midpoint().raw();
        Float cu=this->_constraint_models[i].upper_bound().midpoint().raw();
        b[i]=fb;
        l[nx+i]=cl-fe;
        u[nx+i]=cu+fe;
    }

    if(xa.error()==0.0) {
        l[nx+nc]=0.0; u[nx+nc]=0.0;
    } else {
        l[nx+nc]=-1.0; u[nx+nc]=+1.0;
    }
    if(ya.error()==0.0) {
        l[nx+nc+1]=0.0; u[nx+nc+1]=0.0;
    } else {
        l[nx+nc+1]=-1.0; u[nx+nc+1]=+1.0;
    }
    ARIADNE_LOG(3," A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n");

    List<Point2d> vertices;

    // Set up simplex algorithm working variables
    Array<Slackness> vt(0); Array<size_t> p(nc); Matrix<Float> B(nc,nc);
    Vector<Float> x(np); Vector<Float> y(nc);

    // Find an initial feasible point
    tribool feasible = lpsolver.hotstarted_feasible(l,u,A,b, vt, p,B, x,y);
    ARIADNE_LOG(3," A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<"\n  x="<<x<<" Ax="<< A*x <<"\n");
    lpsolver.consistency_check(l,u,A,b, vt,p,B,x);

    // If problem not feasible, then set is empty; return empty list
    if(!possibly(feasible)) { return vertices; }

    Vector<Float> c(np,0.0);
    for(uint j=0; j!=np; ++j) { c[j]=G[0][j]; }

    // Find a point on the boundary; choose the point minimising the spacial x-coordinate
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
    Vector<Float> last_vec({0.0,-1.0}); // The direction in space of the last step along the boundary
                                        // Should be set orthogonal to the direction (+1,0) which is minimised in finding first boundary point
    Vector<Float> best_next_vec(2);
    Vector<Float> trial_vec(2);
    Vector<Float> Aj(nc); // The jth column of A
    Vector<Float> BAj(nc); // The jth column of A

    uint last_exiting_variable=np; // Last variable to exit the basis

    // FIXME: Do we need check below?
    // if(nx==ne) { vertices.push_back(Point2d(pt[0],pt[1])); return vertices; }

    do {
        ++STEPS;
        Float cot_theta_max=-inf;
        uint s=np; // The index giving the variable x[p[s]] to enter the basis

        // Compute direction the point Gx moves in when variable x[j]=x[p[k]] enters the basis
        // This direction is given by Gd where d_B=-A_B^{-1} A_N e_j and
        //   di=+1/-1 depending on whether variable x[j] is at lower or upper bound
        for(uint k=nc; k!=np; ++k) {
            // Test variable to enter basis; there are m to test, one for each dimension of the domain
            uint j=p[k];
            if(j!=last_exiting_variable || true) {
                uint j=p[k];
                Aj=column(A,j);
                BAj=B*Aj;

                // Compute the direction the point in space moves for x[j] entering the basis
                Vector<Float> d(np,0.0); // The direction x moves in in changing the basis
                d[j] = (vt[j]==LOWER ? +1 : -1);
                for(uint i=0; i!=nc; ++i) {
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

void ValidatedAffineConstrainedImageSet::draw(CanvasInterface& canvas, const Projection2d& projection) const {
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
    for(uint i=1; i!=boundary.size(); ++i) {
        canvas.line_to(boundary[i].x,boundary[i].y);
    }
    canvas.line_to(boundary[0].x,boundary[0].y);
    canvas.fill();
}


std::ostream& ValidatedAffineConstrainedImageSet::write(std::ostream& os) const {
    return os << "ValidatedAffineConstrainedImageSet( domain=" << this->_domain << ", function=" << this->_space_models << ", constraints=" << this->_constraint_models <<" )";
}


} // namespace Ariadne
