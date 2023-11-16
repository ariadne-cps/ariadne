/***************************************************************************
 *            optimization_submodule.cpp
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

#include "pybind11.hpp"
#include "utilities.hpp"

#include "config.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "geometry/box.hpp"

#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "solvers/linear_programming.hpp"
#include "solvers/nonlinear_programming.hpp"
#include "solvers/constraint_solver.hpp"

using namespace Ariadne;




Void export_slackness(pybind11::module& module)
{
    pybind11::enum_<Slackness> variable_enum(module,"Slackness");
    variable_enum.value("BASIS", Slackness::BASIS);
    variable_enum.value("LOWER", Slackness::LOWER);
    variable_enum.value("UPPER", Slackness::UPPER);
}

Void export_constraint(pybind11::module& module)
{
    pybind11::class_<EffectiveConstraint> effective_constraint_class(module,"EffectiveConstraint");
    effective_constraint_class.def(pybind11::init<EffectiveConstraint>());
    effective_constraint_class.def(pybind11::init<EffectiveNumber,EffectiveScalarMultivariateFunction,EffectiveNumber>());
    effective_constraint_class.def(pybind11::init([](EffectiveScalarMultivariateFunction const& f, EffectiveNumber const& u){return EffectiveConstraint(EffectiveNumber(-infty),f,u);}));
    effective_constraint_class.def(pybind11::init([](EffectiveNumber const& l, EffectiveScalarMultivariateFunction const& f){return EffectiveConstraint(l,f,EffectiveNumber(infty));}));
    effective_constraint_class.def(pybind11::init<EffectiveScalarMultivariateFunction,EffectiveNumber>());

    effective_constraint_class.def("__str__",&__cstr__<EffectiveConstraint>);

    pybind11::class_<ValidatedConstraint> validated_constraint_class(module,"ValidatedConstraint");
    validated_constraint_class.def(pybind11::init<ValidatedNumber,ValidatedScalarMultivariateFunction,ValidatedNumber
    >());
    validated_constraint_class.def(pybind11::init<ValidatedConstraint>());
    validated_constraint_class.def(pybind11::init<EffectiveConstraint>());
    validated_constraint_class.def("lower_bound", &ValidatedConstraint::lower_bound);
    validated_constraint_class.def("upper_bound", &ValidatedConstraint::upper_bound);
    validated_constraint_class.def("function", (const ValidatedScalarMultivariateFunction&(ValidatedConstraint::*)()const) &ValidatedConstraint::function);
    validated_constraint_class.def("__str__",&__cstr__<ValidatedConstraint>);
}

template<class P> Void export_optimisation_problem(pybind11::module& module)
{
    pybind11::class_<FeasibilityProblem<P>> feasibility_problem_class(module,(class_name<P>()+"FeasibilityProblem").c_str());
    feasibility_problem_class.def(pybind11::init<BoxType<P>,VectorMultivariateFunction<P>,BoxType<P>>());
    feasibility_problem_class.def("__str__",__cstr__<FeasibilityProblem<P>>);
    feasibility_problem_class.attr("D");
    feasibility_problem_class.attr("g");
    feasibility_problem_class.attr("C");

    pybind11::class_<OptimisationProblem<P>,pybind11::bases<FeasibilityProblem<P>>> optimisation_problem_class(module,(class_name<P>()+"OptimisationProblem").c_str());
    optimisation_problem_class.def(pybind11::init<ScalarMultivariateFunction<P>,BoxType<P>,VectorMultivariateFunction<P>,BoxType<P>>());
    optimisation_problem_class.def("__str__",__cstr__<OptimisationProblem<P>>);
    optimisation_problem_class.attr("D");

    if constexpr (Same<P,ApproximateTag>) {
        feasibility_problem_class.def(pybind11::init<FeasibilityProblem<ValidatedTag>>());
        optimisation_problem_class.def(pybind11::init<OptimisationProblem<ValidatedTag>>());
    }
}

Void export_optimisation_problems(pybind11::module& module)
{
    export_optimisation_problem<ValidatedTag>(module);
    export_optimisation_problem<ApproximateTag>(module);
    pybind11::implicitly_convertible<FeasibilityProblem<ValidatedTag>,FeasibilityProblem<ApproximateTag>>();
    pybind11::implicitly_convertible<OptimisationProblem<ValidatedTag>,OptimisationProblem<ApproximateTag>>();

}


Void export_feasibility_checker_interface(pybind11::module& module)
{
    typedef FeasibilityCheckerInterface::ApproximateVectorType ApproximateVectorType;
    typedef FeasibilityCheckerInterface::ValidatedVectorType ValidatedVectorType;
    typedef FeasibilityCheckerInterface::ExactVectorType ExactVectorType;

    pybind11::class_<FeasibilityCheckerInterface> feasibility_checker_interface_class(module,"FeasibilityCheckerInterface");

    feasibility_checker_interface_class.def("almost_feasible_point", (ApproximateKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ApproximateVectorType, ApproximateNumber)const) &FeasibilityCheckerInterface::almost_feasible_point);
    feasibility_checker_interface_class.def("is_feasible_point", (ValidatedKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ExactVectorType)const) &FeasibilityCheckerInterface::is_feasible_point);
    feasibility_checker_interface_class.def("contains_feasible_point", (ValidatedKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, BoxRangeType)const) &FeasibilityCheckerInterface::contains_feasible_point);
    feasibility_checker_interface_class.def("check_feasibility", (ValidatedKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ValidatedVectorType, ExactVectorType)const) &FeasibilityCheckerInterface::check_feasibility);

    feasibility_checker_interface_class.def("validate_feasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ValidatedVectorType)const) &FeasibilityCheckerInterface::validate_feasibility);
    feasibility_checker_interface_class.def("validate_feasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedVectorMultivariateFunction, ValidatedVectorType)const) &FeasibilityCheckerInterface::validate_feasibility);
    feasibility_checker_interface_class.def("validate_infeasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ApproximateVectorType,ExactVectorType)const) &FeasibilityCheckerInterface::validate_infeasibility);
    feasibility_checker_interface_class.def("validate_infeasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, UpperBoxType, ExactVectorType)const) &FeasibilityCheckerInterface::validate_infeasibility);
    feasibility_checker_interface_class.def("validate_infeasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ExactVectorType)const) &FeasibilityCheckerInterface::validate_infeasibility);

    pybind11::class_<FeasibilityChecker,FeasibilityCheckerInterface> feasibility_checker_class(module,"FeasibilityChecker");
    feasibility_checker_class.def(pybind11::init<>());
    feasibility_checker_class.def("__str__", &__cstr__<FeasibilityChecker>);


}

Void export_optimiser_interface(pybind11::module& module)
{
    typedef OptimiserInterface::ValidatedVectorType ValidatedVectorType;
    typedef OptimiserInterface::ApproximateVectorType ApproximateVectorType;

    pybind11::class_<OptimiserInterface> optimiser_interface_class(module,"OptimiserInterface");
    optimiser_interface_class.def("minimise", (ValidatedVectorType(OptimiserInterface::*)(ValidatedScalarMultivariateFunction, ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType)const) &OptimiserInterface::minimise);
    optimiser_interface_class.def("minimise", (ValidatedVectorType(OptimiserInterface::*)(ValidatedOptimisationProblem)const) &OptimiserInterface::minimise);
    optimiser_interface_class.def("feasible", (ValidatedKleenean(OptimiserInterface::*)(ValidatedFeasibilityProblem)const) &OptimiserInterface::feasible);

    optimiser_interface_class.def("minimise", (ApproximateVectorType(OptimiserInterface::*)(ApproximateOptimisationProblem)const) &OptimiserInterface::minimise);
    optimiser_interface_class.def("feasible", (ApproximateKleenean(OptimiserInterface::*)(ApproximateFeasibilityProblem)const) &OptimiserInterface::feasible);
    //NOTE: Not in C++ API
    //optimiser_interface_class.def("__str__", &__cstr__<OptimiserInterface>);
}

Void export_interior_point_solvers(pybind11::module& module)
{
    using AOP = ApproximateOptimisationProblem;
    using AFP = ApproximateFeasibilityProblem;
    using VOP = ValidatedOptimisationProblem;
    using VFP = ValidatedFeasibilityProblem;

    using AK = ApproximateKleenean;
    using VK = ValidatedKleenean;
    using XA = FloatDPApproximation;
    using VXA = Vector<XA>;
    using XB = FloatDPBounds;
    using VXB = Vector<XB>;
    using AB = ApproximateBoxType;

    {
        using Self=InteriorPointOptimiser;

        pybind11::class_<InteriorPointOptimiser,OptimiserInterface> interior_point_solver_class(module,"InteriorPointOptimiser");
        interior_point_solver_class.def(pybind11::init<>());
        interior_point_solver_class.def("__str__", &__cstr__<InteriorPointOptimiser>);

        interior_point_solver_class.def("minimise_hotstarted", (Tuple<XB,VXB,VXB>(Self::*)(const VOP&, const VXA&, const VXA&)const) &Self::minimise_hotstarted);
        interior_point_solver_class.def("feasible_hotstarted", (Tuple<VK,VXB,VXB>(Self::*)(const VFP&, const VXA&, const VXA&)const) &Self::feasible_hotstarted);

        interior_point_solver_class.def("minimise_hotstarted", (Tuple<XA,VXA,VXA>(Self::*)(const AOP&, const VXA&, const VXA&)const) &Self::minimise_hotstarted);
        interior_point_solver_class.def("feasible_hotstarted", (Tuple<AK,VXA,VXA>(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::feasible_hotstarted);

        interior_point_solver_class.def("minimisation_step", (Void(Self::*)(const AOP&, VXA&, VXA&, VXA&, VXA&, XA&)const) &Self::minimisation_step);
        interior_point_solver_class.def("feasibility_step", (Void(Self::*)(const AFP&, VXA&, VXA&, VXA&, VXA&)const) &Self::feasibility_step);

        interior_point_solver_class.def("compute_dual", (VXA(Self::*)(const AB&, const VXA&, const XA&)const) &Self::compute_dual);
        interior_point_solver_class.def("initial_step_data", (Self::StepData(Self::*)(const AFP&)const) &Self::initial_step_data);
        interior_point_solver_class.def("initial_step_data_hotstarted", (Self::StepData(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::initial_step_data_hotstarted);
        interior_point_solver_class.def("compute_w", (VXA(Self::*)(const AFP&, const VXA&, const XA&)const) &Self::compute_w);
        interior_point_solver_class.def("compute_z", (VXA(Self::*)(const AFP&, const VXA&, const XA&)const) &Self::compute_z);
        interior_point_solver_class.def("compute_mu", (XA(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::compute_mu);

        interior_point_solver_class.def("compute_tz", (Void(Self::*)(const AFP&, VXA&, XA&, VXA&)const) &Self::compute_tz);

        pybind11::class_<InteriorPointOptimiser::StepData> interior_point_solver_step_data_class(module,"InteriorPointOptimiserStepData");
        interior_point_solver_step_data_class.def("__str__", &__cstr__<InteriorPointOptimiser::StepData>);

    }

    {
        using Self=PrimalDualOnlyInteriorPointOptimiser;

        pybind11::class_<PrimalDualOnlyInteriorPointOptimiser,InteriorPointOptimiser> primal_dual_only_interior_point_solver_class(module,"PrimalDualOnlyInteriorPointOptimiser");

        primal_dual_only_interior_point_solver_class.def("feasible_hotstarted", (Tuple<AK,VXA,VXA>(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::feasible_hotstarted);
        primal_dual_only_interior_point_solver_class.def("feasibility_step", (Void(Self::*)(const AFP&, VXA&, VXA&, XA&)const) &Self::feasibility_step);
    }

    {
        using Self = InfeasibleInteriorPointOptimiser;

        pybind11::class_<InfeasibleInteriorPointOptimiser,OptimiserInterface> infeasible_interior_point_solver_class(module,"InfeasibleInteriorPointOptimiser");
        infeasible_interior_point_solver_class.def(pybind11::init<>());
        infeasible_interior_point_solver_class.def("__str__", &__cstr__<InfeasibleInteriorPointOptimiser>);

        infeasible_interior_point_solver_class.def("feasible_hotstarted", (Pair<ValidatedKleenean,VXA>(Self::*)(VFP, const Self::PrimalDualData&)const) & Self::feasible_hotstarted);
        infeasible_interior_point_solver_class.def("setup_feasibility", (Void(Self::*)(const AFP&, Self::StepData&)const) &Self::setup_feasibility);
        infeasible_interior_point_solver_class.def("minimisation_step", (Void(Self::*)(const AOP&, Self::StepData&)const) &Self::minimisation_step);
    }
}

Void export_approximate_optimiser(pybind11::module& module)
{
    pybind11::class_<ApproximateOptimiser,InteriorPointOptimiser> approximate_optimiser_class(module,"ApproximateOptimiser");
    approximate_optimiser_class.def("feasible_zero", (ValidatedKleenean(ApproximateOptimiser::*)(ExactBoxType, ValidatedVectorMultivariateFunction)const) &ApproximateOptimiser::feasible_zero);
    approximate_optimiser_class.def("feasibility_step", (Void(ApproximateOptimiser::*)(const ExactBoxType&, const ApproximateVectorMultivariateFunction&, FloatDPApproximationVector&, FloatDPApproximationVector& Lambda) const) &ApproximateOptimiser::feasibility_step);
}

Void export_interval_optimiser(pybind11::module& module)
{
    pybind11::class_<IntervalOptimiser,InteriorPointOptimiser> interval_optimiser_class(module,"IntervalOptimiser");
    interval_optimiser_class.def("feasible_zero", (ValidatedKleenean(IntervalOptimiser::*)(ExactBoxType, ValidatedVectorMultivariateFunction)const) &IntervalOptimiser::feasible_zero);
    interval_optimiser_class.def("feasibility_step", (Void(IntervalOptimiser::*)(const FloatDPVector&, const FloatDPVector&, const ValidatedVectorMultivariateFunction&, FloatDPBoundsVector&, FloatDPBoundsVector&, FloatDPBoundsVector&, FloatDPBoundsVector, FloatDPBounds& mu) const) &IntervalOptimiser::feasibility_step);
}

Void export_linear_interior_point_solver(pybind11::module& module)
{
    pybind11::class_<InteriorPointLinearOptimiser> linear_interior_point_solver_class(module,"InteriorPointLinearOptimiser");
    linear_interior_point_solver_class.def(pybind11::init<>());
    linear_interior_point_solver_class.def("minimise", &InteriorPointLinearOptimiser::minimise);
    linear_interior_point_solver_class.def("feasible", (ValidatedKleenean(InteriorPointLinearOptimiser::*)(const Vector<FloatDP>&,const Vector<FloatDP>&, const Matrix<FloatDP>&,const Vector<FloatDP>&)const) &InteriorPointLinearOptimiser::feasible);
    linear_interior_point_solver_class.def("validate_feasibility", &InteriorPointLinearOptimiser::validate_feasibility);
}


Void export_constraint_solver(pybind11::module& module)
{
    pybind11::class_<ConstraintSolver> constraint_solver_class(module,"ConstraintSolver");
    constraint_solver_class.def(pybind11::init<>());
    constraint_solver_class.def("hull_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&)const) &ConstraintSolver::hull_reduce);
    constraint_solver_class.def("box_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&,SizeType)const) &ConstraintSolver::box_reduce);
    constraint_solver_class.def("monotone_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&,SizeType)const) &ConstraintSolver::monotone_reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const List<ValidatedConstraint>&)const) &ConstraintSolver::reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const) &ConstraintSolver::reduce);
}


template<class X>
Void export_simplex_solver(pybind11::module& module)
{
    typedef Array<SizeType> SizeArray;

    typedef RigorousNumericType<X> XX;

    pybind11::class_< SimplexLinearOptimiser<X> > simplex_solver_class(module,(class_name<X>()+"SimplexLinearOptimiser").c_str());
    simplex_solver_class.def(pybind11::init<>());
    simplex_solver_class.def("lpstep",(Bool(SimplexLinearOptimiser<X>::*)(const Vector<X>&,const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,Array<Slackness>& ,SizeArray&,Matrix<XX>&,Vector<XX>&)const) &SimplexLinearOptimiser<X>::lpstep);


    simplex_solver_class.def("feasible",(ValidatedKleenean(SimplexLinearOptimiser<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&)const) &SimplexLinearOptimiser<X>::feasible);

    simplex_solver_class.def("verify_feasibility",(ValidatedKleenean(SimplexLinearOptimiser<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,const Array<Slackness>&)const) &SimplexLinearOptimiser<X>::verify_feasibility);

    simplex_solver_class.def("compute_basis",(Pair< SizeArray, Matrix<XX> >(SimplexLinearOptimiser<X>::*)(const Matrix<X>&)const) &SimplexLinearOptimiser<X>::compute_basis);

}


Void optimization_submodule(pybind11::module& module) {
    export_slackness(module);
    export_constraint(module);

    export_optimisation_problems(module);
    export_feasibility_checker_interface(module);
    export_optimiser_interface(module);

    export_interior_point_solvers(module);
    export_approximate_optimiser(module);
    export_interval_optimiser(module);

    export_simplex_solver<FloatDPApproximation>(module);
    export_simplex_solver<FloatDP>(module);
    export_simplex_solver<Rational>(module);
    export_linear_interior_point_solver(module);

    export_constraint_solver(module);
}
