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

#include "solvers/nonlinear_programming.tpl.hpp"

using namespace Ariadne;




Void export_slackness(pybind11::module& module)
{
    pybind11::enum_<Slackness> variable_enum(module,"Slackness");
    variable_enum.value("BASIS", Slackness::BASIS);
    variable_enum.value("LOWER", Slackness::LOWER);
    variable_enum.value("UPPER", Slackness::UPPER);
}

template<class P> pybind11::class_<ConstraintType<P>> export_constraint(pybind11::module& module)
{
    pybind11::class_<ConstraintType<P>> constraint_class(module,(class_name<P>()+"Constraint").c_str());
    constraint_class.def(pybind11::init<ConstraintType<P>>());
    constraint_class.def(pybind11::init<Number<P>,ScalarMultivariateFunction<P>,Number<P>>());
    constraint_class.def(pybind11::init([](ScalarMultivariateFunction<P> const& f, Number<P> const& u){return ConstraintType<P>(Number<P>(-infty),f,u);}));
    constraint_class.def(pybind11::init([](Number<P> const& l, ScalarMultivariateFunction<P> const& f){return ConstraintType<P>(l,f,Number<P>(infty));}));
    constraint_class.def(pybind11::init<ScalarMultivariateFunction<P>,Number<P>>());
    constraint_class.def("lower_bound", &ConstraintType<P>::lower_bound);
    constraint_class.def("upper_bound", &ConstraintType<P>::upper_bound);
    constraint_class.def("function", (const ScalarMultivariateFunction<P>&(ConstraintType<P>::*)()const) &ConstraintType<P>::function);
    constraint_class.def("__str__",&__cstr__<ConstraintType<P>>);
    constraint_class.def("__bool__",[](ConstraintType<P> z){throw pybind11::type_error("Cannot use " + class_name<P>() + "Constraint" + " object in boolean context");});

    pybind11::class_<LowerConstraintType<P>,pybind11::bases<ConstraintType<P>>> lower_constraint_class(module,(class_name<P>()+"LowerConstraint").c_str());
    lower_constraint_class.def("__le__", (ConstraintType<P>(*)(LowerConstraintType<P> const&, Number<P>const&)) operator<=);
    lower_constraint_class.def("__str__",&__cstr__<LowerConstraintType<P>>);

    pybind11::class_<UpperConstraintType<P>,pybind11::bases<ConstraintType<P>>> upper_constraint_class(module,(class_name<P>()+"UpperConstraint").c_str());
    upper_constraint_class.def("__ge__", (ConstraintType<P>(*)(UpperConstraintType<P> const&, Number<P>const&)) operator>=);
    upper_constraint_class.def("__str__",&__cstr__<UpperConstraintType<P>>);

    pybind11::class_<EqualityConstraintType<P>,pybind11::bases<ConstraintType<P>>> equality_constraint_class(module,(class_name<P>()+"EqualityConstraint").c_str());
    equality_constraint_class.def("__str__",&__cstr__<EqualityConstraintType<P>>);

    return constraint_class;
}

Void export_constraints(pybind11::module& module)
{
    auto effective_constraint_class=export_constraint<EffectiveTag>(module);
    auto validated_constraint_class=export_constraint<ValidatedTag>(module);
    auto approximate_constraint_class=export_constraint<ApproximateTag>(module);

    validated_constraint_class.def(pybind11::init<EffectiveConstraint>());
    pybind11::implicitly_convertible<EffectiveConstraint,ValidatedConstraint>();

    approximate_constraint_class.def(pybind11::init<EffectiveConstraint>());
    pybind11::implicitly_convertible<EffectiveConstraint,ApproximateConstraint>();
    approximate_constraint_class.def(pybind11::init<ValidatedConstraint>());
    pybind11::implicitly_convertible<ValidatedConstraint,ApproximateConstraint>();
}

template<class P> Void export_optimisation_problem(pybind11::module& module)
{
    pybind11::class_<FeasibilityProblem<P>> feasibility_problem_class(module,(class_name<P>()+"FeasibilityProblem").c_str());
    feasibility_problem_class.def(pybind11::init<BoxType<P>,VectorMultivariateFunction<P>,BoxType<P>>());
    feasibility_problem_class.def(pybind11::init<BoxType<P>,const List<ConstraintType<P>>&>());
    feasibility_problem_class.def("__str__",__cstr__<FeasibilityProblem<P>>);
    feasibility_problem_class.def_readwrite("D", &FeasibilityProblem<P>::D);
    feasibility_problem_class.def_readwrite("g", &FeasibilityProblem<P>::g);
    feasibility_problem_class.def_readwrite("C", &FeasibilityProblem<P>::C);
    feasibility_problem_class.def("number_of_constraints", &FeasibilityProblem<P>::number_of_constraints);
    feasibility_problem_class.def("number_of_variables", &FeasibilityProblem<P>::number_of_variables);

    pybind11::class_<OptimisationProblem<P>,pybind11::bases<FeasibilityProblem<P>>> optimisation_problem_class(module,(class_name<P>()+"OptimisationProblem").c_str());
    optimisation_problem_class.def(pybind11::init<ScalarMultivariateFunction<P>,BoxType<P>,VectorMultivariateFunction<P>,BoxType<P>>());
    optimisation_problem_class.def(pybind11::init<ScalarMultivariateFunction<P>,BoxType<P>,const List<ConstraintType<P>>&>());
    optimisation_problem_class.def(pybind11::init<ScalarMultivariateFunction<P>,BoxType<P>>());
    optimisation_problem_class.def("__str__",__cstr__<OptimisationProblem<P>>);
    optimisation_problem_class.def_readwrite("f", &OptimisationProblem<P>::f);

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
    typedef FeasibilityCheckerInterface::ApproximateVector ApproximateVector;
    typedef FeasibilityCheckerInterface::ValidatedVector ValidatedVector;
    typedef FeasibilityCheckerInterface::ExactVector ExactVector;

    pybind11::class_<FeasibilityCheckerInterface> feasibility_checker_interface_class(module,"FeasibilityCheckerInterface");

    feasibility_checker_interface_class.def("almost_feasible_point", (ApproximateKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ApproximateVector, ApproximateNumber)const) &FeasibilityCheckerInterface::almost_feasible_point);
    feasibility_checker_interface_class.def("is_feasible_point", (ValidatedKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ExactVector)const) &FeasibilityCheckerInterface::is_feasible_point);
    feasibility_checker_interface_class.def("contains_feasible_point", (ValidatedKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, BoxRangeType)const) &FeasibilityCheckerInterface::contains_feasible_point);
    feasibility_checker_interface_class.def("check_feasibility", (ValidatedKleenean(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ValidatedVector, ExactVector)const) &FeasibilityCheckerInterface::check_feasibility);

    feasibility_checker_interface_class.def("validate_feasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ValidatedVector)const) &FeasibilityCheckerInterface::validate_feasibility);
    feasibility_checker_interface_class.def("validate_feasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedVectorMultivariateFunction, ValidatedVector)const) &FeasibilityCheckerInterface::validate_feasibility);
    feasibility_checker_interface_class.def("validate_infeasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ApproximateVector,ExactVector)const) &FeasibilityCheckerInterface::validate_infeasibility);
    feasibility_checker_interface_class.def("validate_infeasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, UpperBoxType, ExactVector)const) &FeasibilityCheckerInterface::validate_infeasibility);
    feasibility_checker_interface_class.def("validate_infeasibility", (Bool(FeasibilityCheckerInterface::*)(ValidatedFeasibilityProblem, ExactVector)const) &FeasibilityCheckerInterface::validate_infeasibility);

    pybind11::class_<FeasibilityChecker,FeasibilityCheckerInterface> feasibility_checker_class(module,"FeasibilityChecker");
    feasibility_checker_class.def(pybind11::init<>());
    feasibility_checker_class.def("__str__", &__cstr__<FeasibilityChecker>);


}

Void export_optimiser_interface(pybind11::module& module)
{
    typedef ValidatedOptimiserInterface::ValidatedVector ValidatedVector;
    typedef ApproximateOptimiserInterface::ApproximateVector ApproximateVector;

    pybind11::class_<ValidatedOptimiserInterface> validated_optimiser_interface_class(module,"ValidatedOptimiserInterface");
    validated_optimiser_interface_class.def("minimise", (ValidatedVector(ValidatedOptimiserInterface::*)(ValidatedScalarMultivariateFunction, ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType)const) &ValidatedOptimiserInterface::minimise);
    validated_optimiser_interface_class.def("minimise", (ValidatedVector(ValidatedOptimiserInterface::*)(ValidatedOptimisationProblem)const) &ValidatedOptimiserInterface::minimise);
    validated_optimiser_interface_class.def("feasible", (ValidatedKleenean(ValidatedOptimiserInterface::*)( ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType)const) &ValidatedOptimiserInterface::feasible);
    validated_optimiser_interface_class.def("feasible", (ValidatedKleenean(ValidatedOptimiserInterface::*)(ValidatedFeasibilityProblem)const) &ValidatedOptimiserInterface::feasible);

    pybind11::class_<ApproximateOptimiserInterface> approximate_optimiser_interface_class(module,"ApproximateOptimiserInterface");
    approximate_optimiser_interface_class.def("minimise", (ApproximateVector(ApproximateOptimiserInterface::*)(ApproximateScalarMultivariateFunction, ApproximateBoxType, ApproximateVectorMultivariateFunction, ApproximateBoxType)const) &ApproximateOptimiserInterface::minimise);
    approximate_optimiser_interface_class.def("minimise", (ApproximateVector(ApproximateOptimiserInterface::*)(ApproximateOptimisationProblem)const) &ApproximateOptimiserInterface::minimise);
//    approximate_optimiser_interface_class.def("feasible", (ApproximateKleenean(ApproximateOptimiserInterface::*)(ApproximateBoxType, ApproximateVectorMultivariateFunction, ApproximateBoxType)const) &ApproximateOptimiserInterface::feasible);
    approximate_optimiser_interface_class.def("feasible", (ApproximateKleenean(ApproximateOptimiserInterface::*)(ApproximateFeasibilityProblem)const) &ApproximateOptimiserInterface::feasible);
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
    using YA = ApproximateNumber;
    using YB = ValidatedNumber;
    using XA = InteriorPointOptimiser::ApproximateNumberType;
    //using XB = KarushKuhnTuckerOptimiser::ValidatedNumberType;
    using VYA = Vector<YA>;
    using VYB = Vector<YB>;
    using VXA = Vector<XA>;
    //using VXB = Vector<XB>;
    using AB = ApproximateBoxType;
    //using EB = ExactBoxType;

    pybind11::class_<ValuePrimalDualData<YA>> value_primal_dual_data_class(module,"ApproximateValuePrimalDualData");
    value_primal_dual_data_class.def("value",&ValuePrimalDualData<YA>::value);
    value_primal_dual_data_class.def("primal",&ValuePrimalDualData<YA>::primal);
    value_primal_dual_data_class.def("dual",&ValuePrimalDualData<YA>::dual);

    pybind11::class_<FeasiblePrimalDualData<YA>> feasible_primal_dual_data_class(module,"ApproximateFeasiblePrimalDualData");
    feasible_primal_dual_data_class.def("is_feasible",&FeasiblePrimalDualData<YA>::is_feasible);
    feasible_primal_dual_data_class.def("primal",&FeasiblePrimalDualData<YA>::primal);
    feasible_primal_dual_data_class.def("dual",&FeasiblePrimalDualData<YA>::dual);

    {
        using Self=InteriorPointOptimiserBase;

        pybind11::class_<InteriorPointOptimiserBase,ApproximateOptimiserInterface> interior_point_optimiser_base_class(module,"InteriorPointOptimiserBase");

        interior_point_optimiser_base_class.def("minimise_hotstarted", (ValuePrimalDualData<YA>(Self::*)(const AOP&, const VXA&, const VXA&)const) &Self::minimise_hotstarted);
        interior_point_optimiser_base_class.def("feasible_hotstarted", (FeasiblePrimalDualData<YA>(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::feasible_hotstarted);

        interior_point_optimiser_base_class.def("compute_dual", (VXA(Self::*)(const AB&, const VXA&, const XA&)const) &Self::compute_dual);

        interior_point_optimiser_base_class.def("compute_w", (VXA(Self::*)(const AFP&, const VXA&, const XA&)const) &Self::compute_w);
        interior_point_optimiser_base_class.def("compute_z", (VXA(Self::*)(const AFP&, const VXA&, const XA&)const) &Self::compute_z);
        interior_point_optimiser_base_class.def("compute_mu", (XA(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::compute_mu);

        interior_point_optimiser_base_class.def("compute_tz", (Void(Self::*)(const AFP&, VXA&, XA&, VXA&)const) &Self::compute_tz);
    }

    {
        using Self=PrimalDualInteriorPointOptimiser;

        pybind11::class_<Self,InteriorPointOptimiserBase> interior_point_optimiser_class(module,"PrimalDualInteriorPointOptimiser");
        interior_point_optimiser_class.def(pybind11::init<>());
        interior_point_optimiser_class.def("__str__", &__cstr__<Self>);

#warning
        //        interior_point_optimiser_class.def("feasibility_step", (Void(Self::*)(const AFP&, Self::StepData&)const) &Self::feasibility_step);

        interior_point_optimiser_class.def("initial_step_data", (Self::StepData*(Self::*)(const AFP&)const) &Self::initial_step_data);
        interior_point_optimiser_class.def("initial_step_data_hotstarted", (Self::StepData*(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::initial_step_data_hotstarted);
        interior_point_optimiser_class.def("minimisation_step", (Void(Self::*)(const AOP&, Self::StepData&)const) &Self::minimisation_step);
    }

    {
        using Self=PrimalDualComplementaryInteriorPointOptimiser;

        pybind11::class_<Self,InteriorPointOptimiserBase> interior_point_optimiser_class(module,"PrimalDualComplementaryInteriorPointOptimiser");
        interior_point_optimiser_class.def(pybind11::init<>());
        interior_point_optimiser_class.def("__str__", &__cstr__<Self>);

        module.attr("InteriorPointOptimiser") = module.attr("PrimalDualComplementaryInteriorPointOptimiser");

        interior_point_optimiser_class.def("minimisation_step", (Void(Self::*)(const AOP&, VXA&, VXA&, VXA&, XA&)const) &Self::minimisation_step);
        interior_point_optimiser_class.def("feasibility_step", (Void(Self::*)(const AFP&, VXA&, VXA&, VXA&)const) &Self::feasibility_step);

        interior_point_optimiser_class.def("initial_step_data", (Self::StepData*(Self::*)(const AFP&)const) &Self::initial_step_data);
        interior_point_optimiser_class.def("initial_step_data_hotstarted", (Self::StepData*(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::initial_step_data_hotstarted);
        interior_point_optimiser_class.def("minimisation_step", (Void(Self::*)(const AOP&, Self::StepData&)const) &Self::minimisation_step);

        pybind11::class_<Self::StepData> interior_point_optimiser_step_data_class(module,"PrimalDualComplementaryInteriorPointOptimiserStepData");
        interior_point_optimiser_step_data_class.def_readwrite("x", &Self::StepData::x);
        interior_point_optimiser_step_data_class.def_readwrite("y", &Self::StepData::y);
        interior_point_optimiser_step_data_class.def_readwrite("z", &Self::StepData::z);
        interior_point_optimiser_step_data_class.def("__str__", &__cstr__<Self::StepData>);
    }

    {
        using Self = SlackPrimalDualComplementaryInteriorPointOptimiser;

        pybind11::class_<Self,InteriorPointOptimiserBase> interior_point_optimiser_class(module,"SlackPrimalDualComplementaryInteriorPointOptimiser");
        interior_point_optimiser_class.def(pybind11::init<>());
        interior_point_optimiser_class.def("__str__", &__cstr__<Self>);

        module.attr("InfeasibleInteriorPointOptimiser") = module.attr("SlackPrimalDualComplementaryInteriorPointOptimiser");

        interior_point_optimiser_class.def("minimisation_step", (XA(Self::*)(const AOP&, VXA&, VXA&, VXA&, VXA&, XA&)const) &Self::minimisation_step);
        interior_point_optimiser_class.def("feasibility_step", (Void(Self::*)(const AFP&, VXA&, VXA&, VXA&, VXA&)const) &Self::feasibility_step);

        interior_point_optimiser_class.def("initial_step_data", (Self::StepData*(Self::*)(const AFP&)const) &Self::initial_step_data);
        interior_point_optimiser_class.def("initial_step_data_hotstarted", (Self::StepData*(Self::*)(const AFP&, const VXA&, const VXA&)const) &Self::initial_step_data_hotstarted);
        interior_point_optimiser_class.def("minimisation_step", (XA(Self::*)(const AOP&, Self::StepData&)const) &Self::minimisation_step);

        pybind11::class_<Self::StepData> interior_point_optimiser_step_data_class(module,"SlackPrimalDualComplementaryInteriorPointOptimiserStepData");
        interior_point_optimiser_step_data_class.def(pybind11::init<VXA,VXA,VXA,VXA,XA>());
        interior_point_optimiser_step_data_class.def(pybind11::init<Self::StepData>());
        interior_point_optimiser_step_data_class.def_readwrite("w", &Self::StepData::w);
        interior_point_optimiser_step_data_class.def_readwrite("x", &Self::StepData::x);
        interior_point_optimiser_step_data_class.def_readwrite("y", &Self::StepData::y);
        interior_point_optimiser_step_data_class.def_readwrite("z", &Self::StepData::z);
        interior_point_optimiser_step_data_class.def_readwrite("mu", &Self::StepData::mu);
        interior_point_optimiser_step_data_class.def("assemble", &Self::StepData::assemble);
        interior_point_optimiser_step_data_class.def("__str__", &__cstr__<Self::StepData>);

        using MXA=Matrix<XA>; using SMXA=SymmetricMatrix<XA>; using DMXA=DiagonalMatrix<XA>;
        pybind11::class_<SlackPrimalDualComplementaryMatrix<XA>> interior_point_optimiser_matrix_class(module,"SlackPrimalDualComplementaryMatrix");
#warning
        //interior_point_optimiser_matrix_class.def(pybind11::init<Self::StepData>());
        interior_point_optimiser_matrix_class.def(pybind11::init<SMXA,MXA,DMXA,DMXA,DMXA,DMXA>());
        interior_point_optimiser_matrix_class.def("assemble", &SlackPrimalDualComplementaryMatrix<XA>::assemble);
    }

    {
        using Self=SlackPrimalSplitDualComplementaryInteriorPointOptimiser;

        pybind11::class_<Self,InteriorPointOptimiserBase> interior_point_optimiser_class(module,"SlackPrimalSplitDualComplementaryInteriorPointOptimiser");
        interior_point_optimiser_class.def(pybind11::init<>());
        interior_point_optimiser_class.def("__str__", &__cstr__<Self>);

        module.attr("SplitInfeasibleInteriorPointOptimiser") = module.attr("SlackPrimalSplitDualComplementaryInteriorPointOptimiser");

        interior_point_optimiser_class.def("initial_step_data_hotstarted", &Self::initial_step_data_hotstarted);
        interior_point_optimiser_class.def("minimisation_step", &Self::minimisation_step);

        pybind11::class_<Self::StepData> interior_point_optimiser_step_data_class(module,"SlackPrimalSplitDualComplementaryInteriorPointOptimiserStepData");
        interior_point_optimiser_step_data_class.def_readwrite("w", &Self::StepData::w);
        interior_point_optimiser_step_data_class.def_readwrite("x", &Self::StepData::x);
        interior_point_optimiser_step_data_class.def_readwrite("yl", &Self::StepData::yl);
        interior_point_optimiser_step_data_class.def_readwrite("yu", &Self::StepData::yu);
        interior_point_optimiser_step_data_class.def_readwrite("zl", &Self::StepData::zl);
        interior_point_optimiser_step_data_class.def_readwrite("zu", &Self::StepData::zu);
        interior_point_optimiser_step_data_class.def_readwrite("mu", &Self::StepData::mu);
        interior_point_optimiser_step_data_class.def("__str__", &__cstr__<Self::StepData>);
    }

    {
        using Self=KarushKuhnTuckerOptimiser;

        pybind11::class_<KarushKuhnTuckerOptimiser,ValidatedOptimiserInterface> karush_kuhn_tucker_optimiser_class(module,"KarushKuhnTuckerOptimiser");
        karush_kuhn_tucker_optimiser_class.def(pybind11::init<>());
        //karush_kuhn_tucker_optimiser_class.def("__str__", &__cstr__<KarushKuhnTuckerOptimiser>);

        karush_kuhn_tucker_optimiser_class.def("minimise_hotstarted", (ValuePrimalDualData<YB>(Self::*)(const VOP&, const VXA&, const VXA&)const) &Self::minimise_hotstarted);
        karush_kuhn_tucker_optimiser_class.def("feasible_hotstarted", (FeasiblePrimalDualData<YB>(Self::*)(const VFP&, const VXA&, const VXA&)const) &Self::feasible_hotstarted);
    }

    {
        using Self = InfeasibleKarushKuhnTuckerOptimiser;

        pybind11::class_<InfeasibleKarushKuhnTuckerOptimiser,ValidatedOptimiserInterface> infeasible_karush_kuhn_tucker_optimiser_class(module,"InfeasibleKarushKuhnTuckerOptimiser");
        infeasible_karush_kuhn_tucker_optimiser_class.def(pybind11::init<>());
        //infeasible_karush_kuhn_tucker_optimiser_class.def("_str_", &_cstr_<InfeasibleKarushKuhnTuckerOptimiser>);

        infeasible_karush_kuhn_tucker_optimiser_class.def("minimise_hotstarted", (ValuePrimalDualData<YB>(Self::*)(const VOP&, const VXA&, const VXA&)const) &Self::minimise_hotstarted);
        infeasible_karush_kuhn_tucker_optimiser_class.def("feasible_hotstarted", (FeasiblePrimalDualData<YB>(Self::*)(const VFP&, const SlackPrimalDualData<XA>&)const) &Self::feasible_hotstarted);
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
    pybind11::class_<IntervalOptimiser,KarushKuhnTuckerOptimiser> interval_optimiser_class(module,"IntervalOptimiser");
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
    export_constraints(module);

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
