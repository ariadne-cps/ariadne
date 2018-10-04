/***************************************************************************
 *            solver_submodule.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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

#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "solvers/solver_interface.hpp"
#include "solvers/solver.hpp"
#include "function/taylor_function.hpp"

#include "solvers/integrator_interface.hpp"
#include "solvers/integrator.hpp"
#include "solvers/runge_kutta_integrator.hpp"

using namespace Ariadne;

namespace Ariadne {

typedef Vector<ValidatedNumericType> ValidatedVectorType;
typedef Vector<ApproximateNumericType> ApproximateVectorType;

class SolverWrapper
  : public pybind11::wrapper< SolverInterface >
{
  public:
    SolverInterface* clone() const { return this->get_override("clone")(); }
    Void set_maximum_error(RawFloatDP me) { this->get_override("set_maximum_error")(me); }
    FloatDPValue maximum_error() const { return this->get_override("maximum_error")(); }
    Void set_maximum_number_of_steps(Nat ns) { this->get_override("set_maximum_number_of_steps")(ns); }
    Nat maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    ValidatedVectorType zero(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("zero")(f,bx); }
    ValidatedVectorType fixed_point(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("fixed_point")(f,bx); }
    ValidatedVectorType solve(const ValidatedVectorFunction& f, const ValidatedVectorType& pt) const {
        return this->get_override("solve")(f,pt); }
    ValidatedVectorType solve(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve")(f,bx); }
    ValidatedVectorFunctionModelDP implicit(const ValidatedVectorFunction& f, const ExactBoxType& pd, const ExactBoxType& bx) const {
        return this->get_override("implicit")(f,pd,bx); }
    ValidatedScalarFunctionModelDP implicit(const ValidatedScalarFunction& f, const ExactBoxType& pd, const ExactIntervalType& ivl) const {
        return this->get_override("implicit")(f,pd,ivl); }
    ValidatedVectorFunctionModelDP continuation(const ValidatedVectorFunction& f, const ApproximateVectorType& a, const ExactBoxType& X,  const ExactBoxType& A) const {
        return this->get_override("continuation")(f,a,X,A); }
    Set< ValidatedVectorType > solve_all(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve_all")(f,bx); }
    Void write(OutputStream& os) const { this->get_override("write")(os); }
};


class IntegratorWrapper
  : public pybind11::wrapper<IntegratorInterface>
{
  public:
    IntegratorInterface* clone() const {
        return this->get_override("clone")(); }
    Void set_temporal_order(uint to) {
        this->get_override("set_temporal_order")(to); }
    Void set_maximum_error(double me) {
        this->get_override("set_maximum_error")(me); }
    double maximum_error() const {
        return this->get_override("maximum_error")(); }
    Pair<StepSizeType,UpperBoxType> flow_bounds(const ValidatedVectorFunction& vf, const ExactBoxType& D, const StepSizeType& h) const {
        return this->get_override("flow_bounds")(vf,D,h); }
    ValidatedVectorFunctionModelDP flow_step(const ValidatedVectorFunction& vf, const ExactBoxType& D, StepSizeType& h) const {
        return this->get_override("flow_step")(vf,D,h); }
    ValidatedVectorFunctionModelDP flow_step(const ValidatedVectorFunction& vf, const ExactBoxType& D, const StepSizeType& h, const UpperBoxType& B) const {
        return this->get_override("flow_step")(vf,D,h,B); }
    ValidatedVectorFunctionModelDP flow_to(const ValidatedVectorFunction& vf ,const ExactBoxType& D, const Real& tf) const {
        return this->get_override("flow_to")(vf,D,tf); }
    List<ValidatedVectorFunctionModelDP> flow(const ValidatedVectorFunction& vf, const ExactBoxType& D, const Real& t0, const Real& tf) const {
        return this->get_override("flow")(vf,D,t0,tf); }
    List<ValidatedVectorFunctionModelDP> flow(const ValidatedVectorFunction& vf, const ExactBoxType& D, const Real& tf) const {
        return this->get_override("flow")(vf,D,tf); }
    Void write(OutputStream& os) const {
        this->get_override("write")(os); }
};

} // namespace Ariadne


Void export_solvers(pybind11::module& module)
{
    pybind11::class_<SolverInterface,SolverWrapper> solver_interface_class(module,"SolverInterface");
    solver_interface_class.def("solve", (Vector<ValidatedNumericType>(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&)const) &SolverInterface::solve);
    solver_interface_class.def("implicit",(ValidatedVectorFunctionModelDP(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const ExactBoxType&)const) &SolverInterface::implicit);
    solver_interface_class.def("implicit",(ValidatedScalarFunctionModelDP(SolverInterface::*)(const ValidatedScalarFunction&,const ExactBoxType&,const ExactIntervalType&)const) &SolverInterface::implicit);
    solver_interface_class.def("solve_all",(Set< Vector<ValidatedNumericType> >(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&)const) &SolverInterface::solve_all);
    solver_interface_class.def("__str__",&__cstr__<SolverInterface>);

    pybind11::class_<IntervalNewtonSolver, SolverInterface> interval_newton_solver_class(module,"IntervalNewtonSolver");
    interval_newton_solver_class.def(pybind11::init<double,unsigned int>());

    pybind11::class_<KrawczykSolver, SolverInterface> krawczyk_solver_class(module,"KrawczykSolver");
    krawczyk_solver_class.def(pybind11::init<double,unsigned int>());
}



Void export_integrators(pybind11::module& module)
{
    pybind11::class_<IntegratorInterface,IntegratorWrapper> integrator_interface_class(module,"IntegratorInterface");
    integrator_interface_class.def("flow_bounds",(Pair<StepSizeType,UpperBoxType>(IntegratorInterface::*)(const ValidatedVectorFunction&, const ExactBoxType&, const RawFloatDP&)const)&IntegratorInterface::flow_bounds);
    integrator_interface_class.def("flow_step",(ValidatedVectorFunctionModelDP(IntegratorInterface::*)(const ValidatedVectorFunction&, const ExactBoxType&, StepSizeType&)const)&IntegratorInterface::flow_step);
    integrator_interface_class.def("flow_step",(ValidatedVectorFunctionModelDP(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const StepSizeType&,const UpperBoxType&)const)&IntegratorInterface::flow_step);
    integrator_interface_class.def("flow_to",(ValidatedVectorFunctionModelDP(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&IntegratorInterface::flow_to);
    integrator_interface_class.def("flow",(List<ValidatedVectorFunctionModelDP>(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&IntegratorInterface::flow);
    integrator_interface_class.def("__str__", &__cstr__<IntegratorInterface>);

    pybind11::class_<TaylorPicardIntegrator,IntegratorInterface> taylor_picard_integrator_class(module,"TaylorPicardIntegrator");
    taylor_picard_integrator_class.def(pybind11::init<double>());

    pybind11::class_<TaylorSeriesIntegrator,IntegratorInterface> taylor_series_integrator_class(module,"TaylorSeriesIntegrator");
    taylor_series_integrator_class.def(pybind11::init<double>());
    taylor_series_integrator_class.def("maximum_spacial_order",&TaylorSeriesIntegrator::maximum_spacial_order);
    taylor_series_integrator_class.def("maximum_temporal_order",&TaylorSeriesIntegrator::maximum_temporal_order);
    taylor_series_integrator_class.def("maximum_error",&TaylorSeriesIntegrator::maximum_error);
    taylor_series_integrator_class.def("maximum_step_size",&TaylorSeriesIntegrator::maximum_step_size);
    taylor_series_integrator_class.def("set_maximum_spacial_order",&TaylorSeriesIntegrator::set_maximum_spacial_order);
    taylor_series_integrator_class.def("set_maximum_temporal_order",&TaylorSeriesIntegrator::set_maximum_temporal_order);
    taylor_series_integrator_class.def("set_maximum_error",&TaylorSeriesIntegrator::set_maximum_error);
    taylor_series_integrator_class.def("set_maximum_step_size",&TaylorSeriesIntegrator::set_maximum_step_size);

    pybind11::class_<RungeKutta4Integrator> runge_kutta_4_integrator_class(module,"RungeKutta4Integrator");
    runge_kutta_4_integrator_class.def(pybind11::init<double>());
    runge_kutta_4_integrator_class.def("step", &RungeKutta4Integrator::step);
    runge_kutta_4_integrator_class.def("evolve", &RungeKutta4Integrator::evolve);
}


Void solver_submodule(pybind11::module& module)
{
    export_solvers(module);
    export_integrators(module);
}


