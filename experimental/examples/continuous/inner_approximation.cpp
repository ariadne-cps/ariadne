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

#include "ariadne.hpp"
#include "ariadne_main.hpp"
#include "io/drawer.hpp"
#include "solvers/linear_programming.hpp"
#include "solvers/nonlinear_programming.hpp"
#include "utility/stopwatch.hpp"

#include "glpk.h"

using namespace ConcLog;

using namespace Ariadne;
using namespace std;

struct ParallelLinearisationInterface {
    virtual FloatDP minimise(SizeType i, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const = 0;
    virtual FloatDP maximise(SizeType i, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const = 0;

    virtual ~ParallelLinearisationInterface() = default;
};

struct NativeParallelLinearisation : ParallelLinearisationInterface {

    FloatDP minimise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override {
        auto nv = A.column_size();
        auto c = Vector<FloatDP>::unit(nv,k,DoublePrecision());
        return solve(k,c,A,b,xl,xu);
    }

    FloatDP maximise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override {
        auto nv = A.column_size();
        auto c = -Vector<FloatDP>::unit(nv,k,DoublePrecision());
        return solve(k,c,A,b,xl,xu);
    }

  protected:

    virtual FloatDP solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const = 0;

};

struct NativeSimplexParallelLinearisation : NativeParallelLinearisation {
  protected:
    FloatDP solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override {
        SimplexSolver<FloatDP> solver;
        auto solution = solver.minimise(c,xl,xu,A,b);
        return solution.at(k).value();
    }
};

struct NativeIPMParallelLinearisation : NativeParallelLinearisation {
  protected:
    FloatDP solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override {
        InteriorPointSolver solver;
        auto solution = solver.minimise(c,xl,xu,A,b);
        return get<1>(solution).at(k).raw();
    }
};

struct GLPKParallelLinearisation : ParallelLinearisationInterface {

    FloatDP minimise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override {
        return _solve(GLP_MIN,k,A,b,xl,xu);
    }

    FloatDP maximise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override {
        auto nv = A.column_size();
        auto c = -Vector<FloatDP>::unit(nv,k,DoublePrecision());
        return _solve(GLP_MAX,k,A,b,xl,xu);
    }

    virtual void optimisation_method(glp_prob* lp) const = 0;

  private:

    FloatDP _solve(int dir, SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const {

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
        char filename[30] = "";
        std::string dir_s = (dir == GLP_MIN ? "min" : "max");
        snprintf(filename,30,"problem-x%lu-%s.txt",k+1,dir_s.c_str());
        glp_write_lp(lp,NULL,filename);
        optimisation_method(lp);
        double result = glp_get_col_prim(lp, k+1);
        glp_delete_prob(lp);

        return FloatDP(cast_exact(result),DoublePrecision());
    }
};

struct GLPKSimplexParallelLinearisation : GLPKParallelLinearisation {
    void optimisation_method(glp_prob* lp) const override {
        glp_simplex(lp,NULL);
        auto status = glp_get_status(lp);
        switch (status) {
            case GLP_UNDEF :
                throw std::runtime_error("Undefined solution to linear problem.");
            case GLP_INFEAS :
                throw std::runtime_error("Infeasible solution to linear problem.");
            case GLP_NOFEAS :
                throw std::runtime_error("No feasible solution to linear problem.");
        }
    }
};

struct GLPKIPMParallelLinearisation : GLPKParallelLinearisation {
    void optimisation_method(glp_prob* lp) const override {
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
        }
    }
};

double gamma(LabelledEnclosure const& encl, SizeType idx) {
    auto rng = encl.domain()[idx];
    return rng.width().get_d();
}

double gamma_min(LabelledEnclosure const& inner, LabelledEnclosure const& outer) {
    double result = std::numeric_limits<double>::infinity();
    for (SizeType i=0; i<inner.dimension(); ++i)
        result = min(result,gamma(inner,i)/gamma(outer,i));
    return result;
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

Tuple<Matrix<FloatDP>,Vector<FloatDP>,Vector<FloatDP>,Vector<FloatDP>> construct_problem(ValidatedVectorMultivariateFunction const& f, ExactBoxType const& d) {

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

void print_problem(Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu, SizeType i) {
    CONCLOG_SCOPE_CREATE

    SizeType num_auxiliary = A.row_size();
    SizeType num_structural = A.column_size()-num_auxiliary;

    SizeType coord = 1 + i / (num_auxiliary/2);

    std::ostringstream ss;
    ss << (i % 2 == 0 ? "min" : "max") << " x" << coord << std::endl << std::endl;

    for (SizeType i=0; i<num_auxiliary; ++i) {
        ss << "a" << i+1 << ": ";
        for (SizeType j=num_structural; j>0; j--) {
            if (A.at(i,j-1) != 0)
                ss << (A.at(i,j-1).get_d() > 0 ? "+" : "") << A.at(i,j-1).get_d() << " x" << j << " ";
        }
        ss << "<= " << b.at(i).get_d() << std::endl;
    }

    ss << std::endl;
    for (SizeType i=0; i<num_structural; ++i) {
        ss << xl.at(i).get_d() << " <= x" << i+1 << " <= " << xu.at(i).get_d() << std::endl;
    }
    CONCLOG_PRINTLN(ss.str())
}

ExactBoxType intersection_domain(ValidatedVectorMultivariateFunction const& f, ExactBoxType const& d, SizeType i, std::shared_ptr<ParallelLinearisationInterface> solver) {

    auto problem = construct_problem(f,d);
    auto const& A = get<0>(problem);
    auto const& b = get<1>(problem);
    auto const& xl = get<2>(problem);
    auto const& xu = get<3>(problem);

    //print_problem(A,b,xl,xu,i);

    auto n = d.dimension();

    ExactBoxType q(n,ExactIntervalType::empty_interval());
    for (SizeType p=0;p<n;++p) {
        auto lb = solver->minimise(p,A,b,xl,xu);
        auto ub = solver->maximise(p,A,b,xl,xu);
        if (lb < ub) {
            q[p].set_lower_bound(lb);
            q[p].set_upper_bound(ub);
        }
    }
    return q;
}

ExactBoxType nonlinear_nonintersection_domain(ValidatedVectorMultivariateFunction const& h, ExactBoxType const& d, SizeType i, double scaling) {

    bool minimise = (i % 2 == 0);
    SizeType var_idx = i/2;

    NonlinearInfeasibleInteriorPointOptimiser nonlinear_solver;
    auto obj = ValidatedScalarMultivariateFunction::coordinate(d.dimension(),var_idx);
    if (not minimise) obj = -obj;
    ValidatedVectorMultivariateFunction g(0,d.dimension());

    auto solution = nonlinear_solver.minimise(obj,d,g,h);
    CONCLOG_PRINTLN_AT(1,"Solution for intersection: x" << var_idx << (minimise ? " >= " : " <= ") << solution.at(var_idx))

    auto result = d;
    double value_correction = d[var_idx].width().get_d()*(1.0-scaling);
    if (minimise)
        result[var_idx].set_upper_bound((solution[var_idx].lower_raw()-FloatDP(cast_exact(value_correction),DoublePrecision())).lower_raw());
    else
        result[var_idx].set_lower_bound((solution[var_idx].upper_raw()+FloatDP(cast_exact(value_correction),DoublePrecision())).upper_raw());
    CONCLOG_PRINTLN_AT(1,"Resulting (shrinked) domain for non-intersection: " << result)
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

template<class T> class IntervalSearchInterface {
  public:
    //! \brief Return the current point
    virtual T const& current() const = 0;
    //! \brief Move to the next point according to the current outcome
    //! \return Whether a next point is available
    virtual bool move_next(bool current_outcome) = 0;
    //! \brief Return whether a next point is available
    virtual bool ended() const = 0;
    //! \brief Return whether one solution has been found
    //! \details Even if one solution has been found, the search might not have ended, in case we want to refine it
    virtual bool solution_found() const = 0;
};

template<class T> class DownwardSearch : public IntervalSearchInterface<T> {
  public:
    DownwardSearch(T const& lower, T const& upper, T const& decrement) : _current(upper), _lower(lower), _decrement(decrement), _solution_found(false) { }
    T const& current() const override { return _current; }
    bool move_next(bool current_outcome) override {
        if (current_outcome) {
            _solution_found = true;
            _current = _lower;
            return false;
        } else if (not ended()) {
            _current -= _decrement;
            return true;
        } else {
            return false;
        }
    }

    bool ended() const override { return _current-_lower < _decrement; }
    bool solution_found() const override { return _solution_found; }
  private:
    T _current;
    T _lower;
    T const _decrement;
    bool _solution_found;
};

template<class T> class BisectionSearch : public IntervalSearchInterface<T> {
public:
    BisectionSearch(T const& lower, T const& upper, T const& threshold) :
        _current(upper), _lower(lower), _lower_is_set(false), _upper(upper), _threshold(threshold), _solution_found(false) { }
    T const& current() const override { return _current; }
    bool move_next(bool current_outcome) override {
        if (not ended()) {
            if (current_outcome) {
                _solution_found = true;
                _lower = _current;
                _lower_is_set = true;
                if (_current == _upper) {
                    return false;
                }
            } else {
                _upper = _current;
                if (_current == _lower) {
                    return false;
                }
            }
            _current = (_lower_is_set ? _lower + (_upper-_lower)/2 : _lower);
            _lower_is_set = true;
            return true;
        } else {
            return false;
        }
    }
    bool ended() const override { return _upper-_lower <= _threshold; }
    bool solution_found() const override { return _solution_found; }
private:
    T _current;
    T _lower;
    bool _lower_is_set;
    T _upper;
    T const _threshold;
    bool _solution_found;
};

LabelledEnclosure inner_approximation(LabelledEnclosure const& outer, std::shared_ptr<ParallelLinearisationInterface> solver) {
    auto result = outer;
    result.uniform_error_recondition();
    auto const &outer_function = result.state_function();
    auto outer_domain = result.domain();

    auto n = outer_function.result_size();

    List<ValidatedVectorMultivariateFunctionPatch> boundaries;
    for (SizeType i = 0; i < n; ++i) {
        boundaries.push_back(partial_evaluate(outer_function, i, cast_exact(outer_domain[i].upper_bound().get_d())));
        boundaries.push_back(partial_evaluate(outer_function, i, cast_exact(outer_domain[i].lower_bound().get_d())));
    }

    ExactBoxType I = project(outer_domain, Range(0, n));
    CONCLOG_PRINTLN_VAR(I)
    Vector<Kleenean> verified(boundaries.size());

    for (SizeType i = 0; i < verified.size(); ++i) {
        verified[i] = indeterminate;
    }

    while (true) {
        for (SizeType i = 0; i < boundaries.size(); ++i) {
            if (possibly(not verified[i])) {
                CONCLOG_PRINTLN("Checking boundary " << i << " (outer reach evaluated on the " << (i % 2 == 0 ? "upper" : "lower") << " bound of x" << i/2 << ")")

                auto const& boundary = boundaries.at(i);
                auto outer_extension = embed(outer_function, boundary.domain());
                auto boundary_extension = embed(outer_function.domain(), boundary);

                ValidatedVectorMultivariateFunctionPatch f = outer_extension - boundary_extension;

                auto extended_domain_restriction = product(I,project(outer_domain, Range(n, outer_function.argument_size())),boundary.domain());

                BisectionSearch<double> scaling_search(0.01,0.99,0.01);
                while (not scaling_search.ended()) {
                    CONCLOG_PRINTLN("Trying with scaling " << scaling_search.current())
                    bool current_outcome = false;
                    try {
                        auto non_intersection_dom = nonlinear_nonintersection_domain(f, extended_domain_restriction, i, scaling_search.current());
                        try {
                            auto feasible_dom = intersection_domain(f, non_intersection_dom, i, solver);
                            CONCLOG_PRINTLN("First round feasible domain still not empty: " << feasible_dom)
                            feasible_dom = intersection_domain(f, feasible_dom, i, solver);
                            CONCLOG_PRINTLN("Nonlinear solution is not validated, after 2 rounds found feasible domain " << feasible_dom << ", retrying...")
                        } catch (std::exception& e) {
                            CONCLOG_PRINTLN("Nonlinear solution is ultimately validated, using it as a restriction to the inner domain, resetting all not verified boundaries to try again")
                            I = project(non_intersection_dom, Range(0, n));
                            CONCLOG_PRINTLN_VAR(I)
                            current_outcome = true;
                            verified[i] = true;
                            for (SizeType h = 0; h < verified.size(); ++h) {
                                if (not possibly(verified[h]))
                                    verified[h] = indeterminate;
                            }
                        }
                    } catch (std::exception& e) {
                        CONCLOG_PRINTLN("No feasible nonlinear solution...")
                    }
                    scaling_search.move_next(current_outcome);
                }
                if (not scaling_search.solution_found())
                    verified[i] = false;
                CONCLOG_PRINTLN("Current boundary non-intersection verification status: " << verified)
            }
        }

        bool completed = true;
        for (SizeType i = 0; i < verified.size(); ++i) {
            if (possibly(not verified.at(i)) and possibly(verified.at(i))) {
                completed = false;
                break;
            }
        }
        if (completed) break;
    }

    for (SizeType i = 0; i < verified.size(); ++i) {
        if (definitely(not verified.at(i))) {
            throw std::runtime_error("No inner approximation could be computed");
        }
    }

    auto full_restricted_domain = product(I,project(outer_domain,Range(outer_function.result_size(),outer_function.argument_size())));

    result.restrict(full_restricted_domain);

    return result;
}

LabelledFigure bounded_figure(LabelledEnclosure const& e) {
    auto bx = e.bounding_box().euclidean_set();
    auto vars = List<RealVariable>(e.bounding_box().variables());
    auto x = vars[0];
    auto y = vars[1];
    auto xlb = bx[0].lower_bound().get_d();
    auto xub = bx[0].upper_bound().get_d();
    auto ylb = bx[1].lower_bound().get_d();
    auto yub = bx[1].upper_bound().get_d();

    auto xw = xub-xlb;
    auto yw = yub-ylb;
    xlb = xlb - 0.1*xw;
    xub = xub + 0.1*xw;
    ylb = ylb - 0.1*yw;
    yub = yub + 0.1*yw;

    return LabelledFigure(Axes2d(xlb<=x<=xub,ylb<=y<=yub));
}

LabelledEnclosure brusselator_sample() {
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=-y-1.5_dec*pow(x,2)-0.5_dec*pow(x,3)-0.5_dec,dot(y)=3*x-y});

    Real e1=5/100_q; Real e2=7/100_q;
    RealExpressionBoundedConstraintSet initial_set({1-e1<=x<=1+e1,1-e2<=y<=1+e2});

    StepMaximumError max_err=1e-6;
    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);

    Real evolution_time = 1.0_dec;

    LabelledFigure fig=LabelledFigure({0.8_dec<=x<=1.1_dec,0.9_dec<=y<=1.2_dec});

    return evolver.orbit(initial_set,evolution_time,Semantics::UPPER).final()[0];
}

LabelledEnclosure article_sample() {
    using VFT = ValidatedVectorMultivariateTaylorFunctionModelDP;
    using SFT = ValidatedScalarMultivariateTaylorFunctionModelDP;

    RealVariable x1("x1"), x2("x2");
    RealSpace spc({x1,x2});

    ExactBoxType domain({{-1,1},{-1,1},{-1,1},{-1,1}});

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-9);

    auto p0 = SFT::coordinate(domain,0,sweeper);
    auto p1 = SFT::coordinate(domain,1,sweeper);
    auto p2 = SFT::coordinate(domain,2,sweeper);
    auto p3 = SFT::coordinate(domain,3,sweeper);

    auto f = VFT(2,domain,sweeper);
    f[0] = 6.39_x + 1.06_x*p0 + 0.5_x*p1 - 0.02_x*p0*p0 - 0.01_x*p0*p1 + 0.05_x*p2;
    f[1] = 5.6_x + 0.08_x*p0 + 0.92_x*p1 - 0.07_x*p0*p0 - 0.06_x*p0*p1 + 0.04_x*p3;

    auto factory = TaylorFunctionFactory(sweeper);
    EnclosureConfiguration config(factory);

    return {Enclosure(domain,f,config),spc};
}

LabelledEnclosure basic_linear_sample() {
    using VFT = ValidatedVectorMultivariateTaylorFunctionModelDP;
    using SFT = ValidatedScalarMultivariateTaylorFunctionModelDP;

    RealVariable x1("x1"), x2("x2");
    RealSpace spc({x1,x2});

    ExactBoxType domain({{-1,1},{-1,1}});

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-9);

    auto p0 = SFT::coordinate(domain,0,sweeper);
    auto p1 = SFT::coordinate(domain,1,sweeper);

    auto f = VFT(2,domain,sweeper);
    f[0] = p0;
    f[1] = p1;

    auto factory = TaylorFunctionFactory(sweeper);
    EnclosureConfiguration config(factory);

    return {Enclosure(domain,f,config),spc};
}

LabelledEnclosure basic_linear_noisy_sample() {
    using VFT = ValidatedVectorMultivariateTaylorFunctionModelDP;
    using SFT = ValidatedScalarMultivariateTaylorFunctionModelDP;

    RealVariable x1("x1"), x2("x2");
    RealSpace spc({x1,x2});

    ExactBoxType domain({{-1,1},{-1,1},{-1,1},{-1,1}});

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-9);

    auto p0 = SFT::coordinate(domain,0,sweeper);
    auto p1 = SFT::coordinate(domain,1,sweeper);
    auto p2 = SFT::coordinate(domain,2,sweeper);
    auto p3 = SFT::coordinate(domain,3,sweeper);

    auto f = VFT(2,domain,sweeper);

    ExactDouble noise_level = 1e-3_x;

    f[0] = p0+noise_level*p2;
    f[1] = p1+noise_level*p3;

    auto factory = TaylorFunctionFactory(sweeper);
    EnclosureConfiguration config(factory);

    return {Enclosure(domain,f,config),spc};
}

LabelledEnclosure vanderpol_sample() {
    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    StepMaximumError max_err=1e-6;
    GradedTaylorSeriesIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e2);
    CONCLOG_PRINTLN(evolver.configuration());

    Real x0 = 1.40_dec;
    Real y0 = 2.40_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    CONCLOG_PRINTLN("Initial set: " << initial_set);
    Real evolution_time = 0.3_dec;

    return evolver.orbit(initial_set,evolution_time,Semantics::UPPER).final()[0];
}

void ariadne_main() {

    //std::shared_ptr<ParallelLinearisationInterface> solver(new NativeSimplexParallelLinearisation());
    //std::shared_ptr<ParallelLinearisationInterface> solver(new NativeIPMParallelLinearisation());
    std::shared_ptr<ParallelLinearisationInterface> solver(new GLPKSimplexParallelLinearisation());
    //std::shared_ptr<ParallelLinearisationInterface> solver(new GLPKIPMParallelLinearisation());

    auto outer_final = vanderpol_sample();

    CONCLOG_PRINTLN_AT(1,"enclosure function = " << outer_final.state_function())

    auto fig = bounded_figure(outer_final);

    bool inner_found = false;
    auto inner_final = outer_final;
    try {
        Stopwatch<Milliseconds> sw;
        inner_final = inner_approximation(outer_final, solver);
        sw.click();
        CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.");

        inner_found = true;
        auto gamma = gamma_min(inner_final, outer_final);
        CONCLOG_PRINTLN_VAR(gamma)

    } catch (std::exception& e) {
        CONCLOG_PRINTLN("Inner approximation could not be found")
    }

    GraphicsManager::instance().set_drawer(AffineDrawer(7));
    fig << fill_colour(lightgrey) << outer_final << fill_colour(red) << line_colour(red) << line_width(3.0);

    auto outer_final_boundary = boundary(outer_final);
    for (auto const &encl: outer_final_boundary)
        fig << encl;

    if (inner_found) fig << line_colour(black) << line_width(1.0) << fill_colour(orange) << inner_final;
    CONCLOG_RUN_AT(2,fig.write("inner_approximation"));

}