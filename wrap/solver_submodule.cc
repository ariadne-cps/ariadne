/***************************************************************************
 *            solver_submodule.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include "boost_python.h"
#include "utilities.h"

#include <boost/python.hpp>

#include "function.h"
#include "solver_interface.h"
#include "solver.h"
#include "taylor_function.h"

#include "integrator_interface.h"
#include "integrator.h"
#include "runge_kutta_integrator.h"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {

typedef Vector<ValidatedNumber> ValidatedPointType;
typedef Vector<ApproximateNumber> ApproximatePointType;

class SolverWrapper
  : public SolverInterface, public wrapper< SolverInterface >
{
  public:
    SolverInterface* clone() const { return this->get_override("clone")(); }
    void set_maximum_error(double) { this->get_override("set_maximum_error")(); }
    double maximum_error() const { return this->get_override("maximum_error")(); }
    void set_maximum_number_of_steps(uint) { this->get_override("set_maximum_number_of_steps")(); }
    uint maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    ValidatedPointType zero(const ValidatedVectorFunction& f, const ExactBox& bx) const {
        return this->get_override("zero")(); }
    ValidatedPointType fixed_point(const ValidatedVectorFunction& f, const ExactBox& bx) const {
        return this->get_override("fixed_point")(); }
    ValidatedPointType solve(const ValidatedVectorFunction& f, const ValidatedPointType& pt) const {
        return this->get_override("solve")(); }
    ValidatedPointType solve(const ValidatedVectorFunction& f, const ExactBox& bx) const {
        return this->get_override("solve")(); }
    ValidatedVectorFunctionModel implicit(const ValidatedVectorFunction& f, const ExactBox& pd, const ExactBox& bx) const {
        return this->get_override("implicit")(); }
    ValidatedScalarFunctionModel implicit(const ValidatedScalarFunction& f, const ExactBox& pd, const ExactInterval& ivl) const {
        return this->get_override("implicit")(); }
    ValidatedVectorFunctionModel continuation(const ValidatedVectorFunction& f, const ApproximatePointType& a, const ExactBox& X,  const ExactBox& A) const {
        return this->get_override("continuation")(); }
    Set< ValidatedPointType > solve_all(const ValidatedVectorFunction& f, const ExactBox& bx) const {
        return this->get_override("solve_all")(); }
    void write(std::ostream&) const { this->get_override("write")(); }
};


class IntegratorWrapper
  : public IntegratorInterface, public wrapper< IntegratorInterface >
{
  public:
    IntegratorInterface* clone() const {
        return this->get_override("clone")(); }
    void set_temporal_order(uint) {
        this->get_override("set_temporal_order")(); }
    void set_maximum_error(double) {
        this->get_override("set_maximum_error")(); }
    double maximum_error() const {
        return this->get_override("maximum_error")(); }
    Pair<ExactFloat,UpperBox> flow_bounds(const ValidatedVectorFunction&,const ExactBox&,const RawFloat&) const {
        return this->get_override("flow_bounds")(); }
    ValidatedVectorFunctionModel flow_step(const ValidatedVectorFunction&,const ExactBox&,RawFloat&) const {
        return this->get_override("flow_step")(); }
    ValidatedVectorFunctionModel flow_step(const ValidatedVectorFunction&,const ExactBox&,const ExactFloat&,const UpperBox&) const {
        return this->get_override("flow_step")(); }
    ValidatedVectorFunctionModel flow_to(const ValidatedVectorFunction& vector_field,const ExactBox&,const Real&) const {
        return this->get_override("flow_to")(); }
    List<ValidatedVectorFunctionModel> flow(const ValidatedVectorFunction&,const ExactBox&,const Real&,const Real&) const {
        return this->get_override("flow")(); }
    List<ValidatedVectorFunctionModel> flow(const ValidatedVectorFunction&,const ExactBox&,const Real&) const {
        return this->get_override("flow")(); }
    void write(std::ostream&) const {
        this->get_override("write")(); }
};


} // namespace Ariadne


void export_solver()
{
    class_<SolverWrapper, boost::noncopyable> solver_wrapper_class("SolverInterface");
    solver_wrapper_class.def("solve",pure_virtual((Vector<ValidatedNumber>(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBox&)const) &SolverInterface::solve));
    solver_wrapper_class.def("implicit",pure_virtual((ValidatedVectorFunctionModel(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBox&,const ExactBox&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("implicit",pure_virtual((ValidatedScalarFunctionModel(SolverInterface::*)(const ValidatedScalarFunction&,const ExactBox&,const ExactInterval&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("solve_all",pure_virtual((Set< Vector<ValidatedNumber> >(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBox&)const) &SolverInterface::solve_all));
    //solver_wrapper_class.def(self_ns::str(self));

    class_<IntervalNewtonSolver, bases<SolverInterface> > interval_newton_solver_class("IntervalNewtonSolver",init<double,unsigned int>());
    class_<KrawczykSolver, bases<SolverInterface> > krawczyk_solver_class("KrawczykSolver",init<double,unsigned int>());
}



void export_integrator()
{
    class_<IntegratorWrapper, boost::noncopyable> integrator_wrapper_class("IntegratorInterface");
    integrator_wrapper_class.def("flow_bounds",(Pair<ExactFloat,UpperBox>(IntegratorInterface::*)(const ValidatedVectorFunction&, const ExactBox&, const RawFloat&)const)&IntegratorInterface::flow_bounds);
    integrator_wrapper_class.def("flow_step",(ValidatedVectorFunctionModel(IntegratorInterface::*)(const ValidatedVectorFunction&, const ExactBox&, RawFloat&)const)&IntegratorInterface::flow_step);
    integrator_wrapper_class.def("flow_step",(ValidatedVectorFunctionModel(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBox&,const ExactFloat&,const UpperBox&)const)&IntegratorInterface::flow_step);
    integrator_wrapper_class.def("flow_to",(ValidatedVectorFunctionModel(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBox&,const Real&)const)&IntegratorInterface::flow_to);
    integrator_wrapper_class.def("flow",(List<ValidatedVectorFunctionModel>(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBox&,const Real&)const)&IntegratorInterface::flow);


    class_<TaylorPicardIntegrator > taylor_picard_integrator_class("TaylorPicardIntegrator",init<double>());
    taylor_picard_integrator_class.def("flow_bounds",(Pair<Float,UpperBox>(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBox&,const Float&)const)&TaylorPicardIntegrator::flow_bounds);
    taylor_picard_integrator_class.def("flow_step", (ValidatedVectorFunctionModel(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,RawFloat&)const)&TaylorPicardIntegrator::flow_step);
    taylor_picard_integrator_class.def("flow_step", (ValidatedVectorFunctionModel(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,const ExactFloat&,const UpperBox&)const)&TaylorPicardIntegrator::flow_step);
    taylor_picard_integrator_class.def("flow_to",(ValidatedVectorFunctionModel(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,const Real&)const)&TaylorPicardIntegrator::flow_to);
    taylor_picard_integrator_class.def("flow",(List<ValidatedVectorFunctionModel>(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,const Real&)const)&TaylorPicardIntegrator::flow);

    class_<TaylorSeriesIntegrator > taylor_series_integrator_class("TaylorSeriesIntegrator",init<double>());
    taylor_series_integrator_class.def("maximum_spacial_order",&TaylorSeriesIntegrator::maximum_spacial_order);
    taylor_series_integrator_class.def("maximum_temporal_order",&TaylorSeriesIntegrator::maximum_temporal_order);
    taylor_series_integrator_class.def("maximum_error",&TaylorSeriesIntegrator::maximum_error);
    taylor_series_integrator_class.def("maximum_step_size",&TaylorSeriesIntegrator::maximum_step_size);
    taylor_series_integrator_class.def("set_maximum_spacial_order",&TaylorSeriesIntegrator::set_maximum_spacial_order);
    taylor_series_integrator_class.def("set_maximum_temporal_order",&TaylorSeriesIntegrator::set_maximum_temporal_order);
    taylor_series_integrator_class.def("set_maximum_error",&TaylorSeriesIntegrator::set_maximum_error);
    taylor_series_integrator_class.def("set_maximum_step_size",&TaylorSeriesIntegrator::set_maximum_step_size);
    taylor_series_integrator_class.def("flow_bounds",(Pair<ExactFloat,UpperBox>(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBox&,const RawFloat&)const)&TaylorSeriesIntegrator::flow_bounds);
    taylor_series_integrator_class.def("flow_step", (ValidatedVectorFunctionModel(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,RawFloat&)const)&TaylorSeriesIntegrator::flow_step);
    taylor_series_integrator_class.def("flow_step", (ValidatedVectorFunctionModel(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,const ExactFloat&,const UpperBox&)const)&TaylorSeriesIntegrator::flow_step);
    taylor_series_integrator_class.def("flow_to",(ValidatedVectorFunctionModel(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,const Real&)const)&TaylorSeriesIntegrator::flow_to);
    taylor_series_integrator_class.def("flow",(List<ValidatedVectorFunctionModel>(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBox&,const Real&)const)&TaylorSeriesIntegrator::flow);

    class_<RungeKutta4Integrator > runge_kutta_4_integrator_class("RungeKutta4Integrator",init<double>());
    runge_kutta_4_integrator_class.def("step", &RungeKutta4Integrator::step);
    runge_kutta_4_integrator_class.def("evolve", &RungeKutta4Integrator::evolve);
}


void solver_submodule()
{
    to_python_list< Set< ExactBox > >();
    to_python< List< ExactBox > >();
    to_python< Pair< Float, ExactBox > >();

    export_solver();
    export_integrator();

}


