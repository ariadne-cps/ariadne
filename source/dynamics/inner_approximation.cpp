/***************************************************************************
 *            inner_approximation.cpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "inner_approximation.hpp"

using namespace ConcLog;

using namespace std;

namespace Ariadne {

FloatDP NativeLinearSolver::minimise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {
    auto nv = A.column_size();
    auto c = Vector<FloatDP>::unit(nv,k,DoublePrecision());
    return solve(k,c,A,b,xl,xu);
}

FloatDP NativeLinearSolver::maximise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {
    auto nv = A.column_size();
    auto c = -Vector<FloatDP>::unit(nv,k,DoublePrecision());
    return solve(k,c,A,b,xl,xu);
}

FloatDP NativeSimplex::solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {
    CONCLOG_SCOPE_CREATE
    SimplexSolver<FloatDP> solver;
    auto solution = solver.minimise(c,xl,xu,A,b);
    return solution.at(k).value();
}

LinearSolverInterface* NativeSimplex::clone() const {
    return new NativeSimplex();
}

FloatDP NativeIPM::solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {
    CONCLOG_SCOPE_CREATE
    InteriorPointSolver solver;
    auto solution = solver.minimise(c,xl,xu,A,b);
    return get<1>(solution).at(k).raw();
}
LinearSolverInterface* NativeIPM::clone() const {
    return new NativeIPM();
}

FloatDP GLPKSolver::minimise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {
    return _solve(GLP_MIN,k,A,b,xl,xu);
}

FloatDP GLPKSolver::maximise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {
    auto nv = A.column_size();
    auto c = -Vector<FloatDP>::unit(nv,k,DoublePrecision());
    return _solve(GLP_MAX,k,A,b,xl,xu);
}

FloatDP GLPKSolver::_solve(int dir, SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {

    SizeType num_auxiliary = A.row_size();
    SizeType num_structural = A.column_size()-num_auxiliary;

    glp_term_out(GLP_OFF);
    glp_prob *lp;
    int ia[1+1000], ja[1+1000];
    double ar[1+1000];
    lp = glp_create_prob();

    glp_set_obj_dir(lp, dir);

    glp_add_rows(lp, num_auxiliary);
    for (SizeType i=1; i<=num_auxiliary; ++i) {
        char str[3] = "";
        snprintf(str,3,"a%lu",i);
        glp_set_row_name(lp, i, str);
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, b.at(i-1).get_d());
    }

    glp_add_cols(lp, num_structural);
    for (SizeType i=1; i<=num_structural; ++i) {
        char str[3] = "";
        snprintf(str,3,"x%lu",i);
        glp_set_col_name(lp, i, str);
        glp_set_col_bnds(lp, i, GLP_DB, xl.at(i-1).get_d(), xu.at(i-1).get_d());
        glp_set_obj_coef(lp,i, 0);
    }
    glp_set_obj_coef(lp,k+1, 1);

    SizeType pr = 1;
    for (SizeType i=0; i<num_auxiliary; ++i) {
        for (SizeType j=0; j<num_structural; ++j) {
            ia[pr] = i+1;
            ja[pr] = j+1;
            ar[pr] = A.at(i,j).get_d();
            ++pr;
        }
    }

    glp_load_matrix(lp, num_auxiliary*num_structural, ia, ja, ar);
    optimisation_method(lp);
    double result = glp_get_col_prim(lp, k+1);
    glp_delete_prob(lp);

    return FloatDP(cast_exact(result),DoublePrecision());
}

void GLPKSimplex::optimisation_method(glp_prob* lp) const {
    glp_simplex(lp,NULL);
    auto status = glp_get_status(lp);
    switch (status) {
        case GLP_UNDEF :
            throw std::runtime_error("Undefined solution to linear problem.");
        case GLP_INFEAS :
            throw std::runtime_error("Infeasible solution to linear problem.");
        case GLP_NOFEAS :
            throw std::runtime_error("No feasible solution to linear problem.");
        default :
            throw std::runtime_error("Unhandled GLP status value " + to_string(status) + ".");
    }
}

LinearSolverInterface* GLPKSimplex::clone() const {
    return new GLPKSimplex();
}

void GLPKIPM::optimisation_method(glp_prob* lp) const {
    glp_interior(lp,NULL);
    auto status = glp_ipt_status(lp);
    switch (status) {
        case GLP_UNDEF :
            throw std::runtime_error("Undefined solution to linear problem.");
        case GLP_INFEAS :
            throw std::runtime_error("Infeasible solution to linear problem.");
        case GLP_NOFEAS :
            throw std::runtime_error("No feasible solution to linear problem.");
        case GLP_UNBND :
            throw std::runtime_error("Unbounded solution to linear problem.");
        default :
            throw std::runtime_error("Unhandled GLP status value " + to_string(status) + ".");
    }
}

LinearSolverInterface* GLPKIPM::clone() const {
    return new GLPKSimplex();
}

List<LabelledEnclosure> boundary(LabelledEnclosure const& enclosure) {
    List<LabelledEnclosure> result;
    auto nv = enclosure.state_space().size();
    for (SizeType i=0; i<nv; ++i) {
        {
            auto boundary_piece  = enclosure;
            auto bounding_domain = boundary_piece.domain();
            bounding_domain[i].set_lower_bound(bounding_domain[i].upper_bound());
            boundary_piece.restrict(bounding_domain);
            result.push_back(boundary_piece);
        }
        {
            auto boundary_piece  = enclosure;
            auto bounding_domain = boundary_piece.domain();
            bounding_domain[i].set_upper_bound(bounding_domain[i].lower_bound());
            boundary_piece.restrict(bounding_domain);
            result.push_back(boundary_piece);
        }
    }
    return result;
}

Tuple<Matrix<FloatDP>,Vector<FloatDP>,Vector<FloatDP>,Vector<FloatDP>> construct_parallel_linearisation_problem(ValidatedVectorMultivariateFunction const& f, ExactBoxType const& d) {

    auto x0 = midpoint(d);

    auto Jx0 = f.jacobian(x0);
    auto J_rng = jacobian_range(f,cast_vector(d));
    auto fx0 = f(x0);
    auto b_rng = fx0 -Jx0 * Vector<FloatDPUpperInterval>(x0) + (J_rng-Jx0) * (d-x0);

    auto n = f.result_size();
    auto m = f.argument_size();
    auto nv = m+2*n;

    Matrix<FloatDP> A(2*n,nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<n; ++i) {
        for (SizeType j=0; j<m; ++j)
            A.at(i,j) = -Jx0.at(i,j).value();
        A.at(i,m+i) = 1;
    }
    for (SizeType i=0; i<n; ++i) {
        for (SizeType j=0; j<m; ++j)
            A.at(n+i,j) = Jx0.at(i,j).value();
        A.at(n+i,m+n+i) = 1;
    }

    Vector<FloatDP> b(2*n,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<n; ++i) {
        b.at(i) = b_rng.at(i).upper_bound().raw();
    }
    for (SizeType i=0; i<n; ++i) {
        b.at(n+i) = -b_rng.at(i).lower_bound().raw();
    }

    Vector<FloatDP> xl(nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<m; ++i) {
        xl.at(i) = d[i].lower_bound();
    }

    Vector<FloatDP> xu(nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<m; ++i) {
        xu.at(i) = d[i].upper_bound();
    }
    for (SizeType i=m; i<nv; ++i) {
        xu.at(i) = inf;
    }

    return std::make_tuple(A,b,xl,xu);
}

ApproximateOptimiser::ValidatedNumericType nonlinear_intersection_bound(ValidatedVectorMultivariateFunctionPatch const& h, ExactBoxType const& d, SizeType var_idx, bool maximise) {
    CONCLOG_SCOPE_CREATE
    NonlinearInfeasibleInteriorPointOptimiser nonlinear_solver;
    auto obj = ValidatedScalarMultivariateFunction::coordinate(d.dimension(),var_idx);
    if (maximise) obj = -obj;
    ValidatedVectorMultivariateFunction g(0,d.dimension());

    auto dom = d;
    // Replace the domain for var_idx with that of h, in order to try the maximum possible range anyway
    dom[var_idx] = h.domain()[var_idx];

    try {
        auto point = nonlinear_solver.minimise(obj,d,g,h);
        return point.at(var_idx);
    } catch (std::exception& e) {
        CONCLOG_PRINTLN("Nonlinear minimisation failed, using original domain bound instead")
        return (maximise ? dom[var_idx].lower_bound() : dom[var_idx].upper_bound());
    }
}

ExactBoxType nonlinear_nonintersection_domain(ApproximateOptimiser::ValidatedNumericType const& bound, ExactBoxType const& d, SizeType var_idx, bool maximise, double scaling) {

    auto result = d;
    double value_correction = d[var_idx].width().get_d()*(1.0-scaling);
    if (maximise)
        result[var_idx].set_lower_bound((bound.upper_raw()+FloatDP(cast_exact(value_correction),DoublePrecision())).upper_raw());
    else
        result[var_idx].set_upper_bound((bound.lower_raw()-FloatDP(cast_exact(value_correction),DoublePrecision())).lower_raw());

    return result;
}

ExactBoxType inner_difference(ExactBoxType const& bx1, ExactBoxType const& bx2) {
    if (not bx1.intersects(bx2)) return bx1;

    auto result = bx1;

    auto n = bx1.dimension();

    SizeType change_idx = 0;
    Interval<FloatDPBounds> change_range(0,0);
    FloatDPBounds vmax = FloatDP(-inf,DoublePrecision());

    for (SizeType i=0; i<n; ++i) {
        auto dl = bx2[i].lower_bound()-bx1[i].lower_bound();
        auto du = bx1[i].upper_bound()-bx2[i].upper_bound();
        auto v = max(dl,du);
        if (definitely(v > vmax)) {
            change_idx = i;
            change_range = (definitely(dl >= du) ? Interval<FloatDPBounds>(bx1[i].lower_bound(),bx2[i].lower_bound()) : Interval<FloatDPBounds>(bx2[i].upper_bound(),bx1[i].upper_bound()));
            vmax = v;
        }
    }

    result[change_idx].set_lower_bound(change_range.lower_bound().lower_raw());
    result[change_idx].set_upper_bound(change_range.upper_bound().upper_raw());

    return result;
}


ParallelLinearisationContractor::ParallelLinearisationContractor(ParallelLinearisationContractor const& other) :
    _solver_ptr(other._solver_ptr) , _num_iterations(other._num_iterations), _split_depth(other._split_depth) { }

ParallelLinearisationContractor::ParallelLinearisationContractor(LinearSolverInterface const& solver, SizeType num_iterations, SizeType split_depth) :
_solver_ptr(solver.clone()), _num_iterations(num_iterations), _split_depth(split_depth) { }

ExactBoxType ParallelLinearisationContractor::contract(ValidatedVectorMultivariateFunctionPatch const& f, ExactBoxType const& d) const {
    CONCLOG_SCOPE_CREATE

    auto n = d.dimension();
    auto m = f.result_size();

    List<ExactBoxType> splits;
    splits.push_back(d);
    for (SizeType k=0; k<_split_depth; ++k)
        for (SizeType i=0; i<m; ++i) {
            List<ExactBoxType> new_splits;
            for (auto const& b : splits) {
                auto s1s2 = split(b,i);
                new_splits.push_back(s1s2.first);
                new_splits.push_back(s1s2.second);
            }
            splits =new_splits;
        }

    List<ExactBoxType> contracted;

    for (auto const& s : splits) {
        auto contraction = s;
        try {
            for (SizeType i=0; i<_num_iterations; ++i) {

                auto problem = construct_parallel_linearisation_problem(f, contraction);
                auto const& A = get<0>(problem);
                auto const& b = get<1>(problem);
                auto const& xl = get<2>(problem);
                auto const& xu = get<3>(problem);

                ExactBoxType q(n,ExactIntervalType::empty_interval());

                for (SizeType p=0; p<n; ++p) {
                    auto lb = _solver_ptr->minimise(p,A,b,xl,xu);
                    auto ub = _solver_ptr->maximise(p,A,b,xl,xu);
                    if (lb < ub) {
                        q[p].set_lower_bound(lb);
                        q[p].set_upper_bound(ub);
                    }
                }
                contraction = q;
            }
            contracted.push_back(contraction);
        } catch (std::exception& e) {
        }
    }

    if (contracted.empty()) return ExactBoxType(n,ExactIntervalType::empty_interval());

    ExactBoxType result = contracted.at(0);
    for (SizeType i=1; i<contracted.size(); ++i) {
        result = hull(result,contracted.at(i));
    }

    return result;
}

ContractorInterface* ParallelLinearisationContractor::clone() const {
    return new ParallelLinearisationContractor(*this);
}

InnerApproximatorBase::InnerApproximatorBase(ContractorInterface const& contractor) :
    _contractor_ptr(contractor.clone()) { }

NonlinearCandidateValidationInnerApproximator::NonlinearCandidateValidationInnerApproximator(ContractorInterface const& contractor, SizeType max_rounds) :
    InnerApproximatorBase(contractor), _max_rounds(max_rounds) { }

LabelledEnclosure NonlinearCandidateValidationInnerApproximator::compute_from(LabelledEnclosure const& outer) const {
    CONCLOG_SCOPE_CREATE
    auto reconditioned_outer = outer;
    reconditioned_outer.uniform_error_recondition();
    auto const& outer_function = reconditioned_outer.state_function();
    auto outer_domain = reconditioned_outer.domain();

    auto n = outer_function.result_size();

    List<ValidatedVectorMultivariateFunctionPatch> boundaries;
    for (SizeType i = 0; i < n; ++i) {
        boundaries.push_back(partial_evaluate(outer_function, i, cast_exact(outer_domain[i].upper_bound().get_d())));
        boundaries.push_back(partial_evaluate(outer_function, i, cast_exact(outer_domain[i].lower_bound().get_d())));
    }

    ExactBoxType I = project(outer_domain, Range(0, n));
    CONCLOG_PRINTLN("Starting domain: " << I)
    Vector<Kleenean> verified(boundaries.size());
    Vector<bool> bound_found(boundaries.size());

    for (SizeType i = 0; i < verified.size(); ++i) {
        verified[i] = indeterminate;
        bound_found[i] = false;
    }

    SizeType rnd = 0;
    SizeType max_rnds = std::max(boundaries.size(),_max_rounds);
    bool completed = false;
    bool all_bounds_found = false;
    while (not completed and rnd < max_rnds) {

        auto bnd_idx = rnd % boundaries.size();

        if (possibly(not verified[bnd_idx])) {

            CONCLOG_PRINTLN_AT(1,"Round " << rnd)

            SizeType var_idx = bnd_idx/2;
            bool is_lower_boundary = (bnd_idx % 2 == 1);

            CONCLOG_PRINTLN_AT(1,"Checking boundary " << bnd_idx << " (outer reach evaluated on the " << (is_lower_boundary ? "lower" : "upper") << " bound of x" << var_idx << ")")

            auto const& boundary = boundaries.at(bnd_idx);
            auto outer_extension = embed(outer_function, boundary.domain());
            auto boundary_extension = embed(outer_function.domain(), boundary);

            ValidatedVectorMultivariateFunctionPatch f = outer_extension - boundary_extension;

            auto extended_domain_restriction = product(I,project(outer_domain, Range(n, outer_function.argument_size())),boundary.domain());

            auto intersection_bound = nonlinear_intersection_bound(f, extended_domain_restriction, var_idx, is_lower_boundary);
            CONCLOG_PRINTLN_AT(2,"Candidate solution for intersection: x" << var_idx << (is_lower_boundary ? " <= " : " >= ") << intersection_bound)

            BisectionSearch<double> scaling_search(0.01,0.99,0.01);
            CONCLOG_PRINTLN_AT(2,"Trying with scaling " << scaling_search.current())
            auto non_intersection_dom = nonlinear_nonintersection_domain(intersection_bound, extended_domain_restriction, var_idx, is_lower_boundary, scaling_search.current());
            if (non_intersection_dom.is_empty()) {
                CONCLOG_PRINTLN_AT(1,"Non-intersection domain from nonlinear optimisation is empty.")
            }

            auto scaled_bound_is_an_improvement = true;
            while (true) {

                CONCLOG_PRINTLN_AT(2,"Resulting candidate (shrinked) domain for non-intersection: " << non_intersection_dom)

                scaled_bound_is_an_improvement =  (not bound_found[bnd_idx] or (is_lower_boundary ? non_intersection_dom[var_idx].lower_bound() < extended_domain_restriction[var_idx].lower_bound() : non_intersection_dom[var_idx].upper_bound() > extended_domain_restriction[var_idx].upper_bound()));

                bool current_outcome = false;

                auto feasible_dom = _contractor_ptr->contract(f, non_intersection_dom);
                if (feasible_dom.is_empty()) {
                    CONCLOG_PRINTLN_AT(2,"Nonlinear solution is validated")
                    current_outcome = true;
                } else {
                    CONCLOG_PRINTLN_AT(2,"Nonlinear solution is not validated: found feasible domain " << feasible_dom << ", retrying...")
                }

                scaling_search.move_next(current_outcome);
                if (scaling_search.ended()) break;

                CONCLOG_PRINTLN_AT(2,"Trying with scaling " << scaling_search.current())
                non_intersection_dom = nonlinear_nonintersection_domain(intersection_bound, extended_domain_restriction, var_idx, is_lower_boundary, scaling_search.current());
                CONCLOG_PRINTLN_AT(1,non_intersection_dom << " empty? " << non_intersection_dom.is_empty())
            }

            if (non_intersection_dom.is_empty() or not scaling_search.solution_found()) {
                if (not bound_found[bnd_idx]) {
                    CONCLOG_PRINTLN_AT(1,"No solution found, setting failure for this boundary.")
                    verified[bnd_idx] = false;
                } else {
                    CONCLOG_PRINTLN_AT(1,"No solution found, keeping the original value for the bound, setting this boundary as verified.")
                    verified[bnd_idx] = true;
                }
            } else {
                if (scaled_bound_is_an_improvement) {
                    bound_found[bnd_idx] = true;
                    SizeType nbf = 0;
                    for (SizeType i=0; i<bound_found.size(); ++i) if (bound_found[i]) ++nbf;
                    if (nbf == bound_found.size()) all_bounds_found = true;

                    CONCLOG_PRINTLN_AT(1,"Using optimal valid solution as a restriction to the inner domain, resetting all remaining boundaries to try again.")
                    CONCLOG_PRINTLN_AT(2,"Value identified: x" << var_idx << (is_lower_boundary ? " >= " : " <= ") << (is_lower_boundary ? non_intersection_dom[var_idx].lower_bound() : non_intersection_dom[var_idx].upper_bound()))
                    I = project(non_intersection_dom, Range(0, n));
                    CONCLOG_PRINTLN_AT(1,"Current candidate: " << I << (not all_bounds_found ? " (incomplete) " : ""))
                    for (SizeType h = 0; h < verified.size(); ++h) {
                        verified[h] = indeterminate;
                    }
                } else {
                    CONCLOG_PRINTLN_AT(1,"No improvement, keeping the original value for the bound, setting this boundary as verified.")
                }
                verified[bnd_idx] = true;
            }
            CONCLOG_PRINTLN_AT(1,"Current boundary non-intersection verification status: " << verified)
        }

        completed = true;
        for (SizeType i = 0; i < verified.size(); ++i) {
            if (possibly(not verified.at(i)) and possibly(verified.at(i))) {
                completed = false;
                break;
            }
        }
        ++rnd;
    }

    for (SizeType i = 0; i < bound_found.size(); ++i) {
        if (not bound_found.at(i)) {
            throw std::runtime_error("No inner approximation could be computed.");
        }
    }

    auto full_restricted_domain = product(I,project(outer.domain(),Range(outer_function.result_size(),outer.state_function().argument_size())));

    CONCLOG_PRINTLN("Domain for inner: " << full_restricted_domain)

    auto result = outer;
    result.restrict(full_restricted_domain);

    return result;
}

}