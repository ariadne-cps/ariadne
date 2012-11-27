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

#include <boost/python.hpp>

#include "function.h"
#include "solver_interface.h"
#include "solver.h"
#include "taylor_function.h"

#include "integrator_interface.h"
#include "integrator.h"
#include "runge_kutta_integrator.h"

#include "utilities.h"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {


class SolverWrapper
  : public SolverInterface, public wrapper< SolverInterface >
{
  public:
    SolverInterface* clone() const { return this->get_override("clone")(); }
    void set_maximum_error(double) { this->get_override("set_maximum_error")(); }
    double maximum_error() const { return this->get_override("maximum_error")(); }
    void set_maximum_number_of_steps(uint) { this->get_override("set_maximum_number_of_steps")(); }
    uint maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    IntervalVector zero(const IntervalVectorFunction& f, const IntervalVector& bx) const {
        return this->get_override("zero")(); }
    IntervalVector fixed_point(const IntervalVectorFunction& f, const IntervalVector& bx) const {
        return this->get_override("fixed_point")(); }
    IntervalVector solve(const IntervalVectorFunction& f, const IntervalVector& bx) const {
        return this->get_override("solve")(); }
    IntervalVectorFunctionModel implicit(const IntervalVectorFunction& f, const IntervalVector& pd, const IntervalVector& bx) const {
        return this->get_override("implicit")(); }
    IntervalScalarFunctionModel implicit(const IntervalScalarFunction& f, const IntervalVector& pd, const Interval& ivl) const {
        return this->get_override("implicit")(); }
    IntervalVectorFunctionModel continuation(const IntervalVectorFunction& f, const Vector<Float>& a, const Vector<Interval>& X,  const Vector<Interval>& A) const {
        return this->get_override("continuation")(); }
    Set< IntervalVector > solve_all(const IntervalVectorFunction& f, const IntervalVector& bx) const {
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
    Pair<Float,IntervalVector> flow_bounds(const IntervalVectorFunction&,const IntervalVector&,const Float&) const {
        return this->get_override("flow_bounds")(); }
    IntervalVectorFunctionModel flow_step(const IntervalVectorFunction&,const IntervalVector&,Float&) const {
        return this->get_override("flow_step")(); }
    IntervalVectorFunctionModel flow_step(const IntervalVectorFunction&,const IntervalVector&,const Float&,const IntervalVector&) const {
        return this->get_override("flow_step")(); }
    IntervalVectorFunctionModel flow_to(const IntervalVectorFunction& vector_field,const IntervalVector&,const Real&) const {
        return this->get_override("flow_to")(); }
    List<IntervalVectorFunctionModel> flow(const IntervalVectorFunction&,const IntervalVector&,const Real&,const Real&) const {
        return this->get_override("flow")(); }
    List<IntervalVectorFunctionModel> flow(const IntervalVectorFunction&,const IntervalVector&,const Real&) const {
        return this->get_override("flow")(); }
    void write(std::ostream&) const {
        this->get_override("write")(); }
};


} // namespace Ariadne


void export_solver()
{
    class_<SolverWrapper, boost::noncopyable> solver_wrapper_class("SolverInterface");
    solver_wrapper_class.def("solve",pure_virtual((IntervalVector(SolverInterface::*)(const IntervalVectorFunction&,const IntervalVector&)const) &SolverInterface::solve));
    solver_wrapper_class.def("implicit",pure_virtual((IntervalVectorFunctionModel(SolverInterface::*)(const IntervalVectorFunction&,const IntervalVector&,const IntervalVector&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("implicit",pure_virtual((IntervalScalarFunctionModel(SolverInterface::*)(const IntervalScalarFunction&,const IntervalVector&,const Interval&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("solve_all",pure_virtual((Set< IntervalVector >(SolverInterface::*)(const IntervalVectorFunction&,const IntervalVector&)const) &SolverInterface::solve_all));
    //solver_wrapper_class.def(self_ns::str(self));

    class_<IntervalNewtonSolver, bases<SolverInterface> > interval_newton_solver_class("IntervalNewtonSolver",init<double,unsigned int>());
    class_<KrawczykSolver, bases<SolverInterface> > krawczyk_solver_class("KrawczykSolver",init<double,unsigned int>());
}



void export_integrator()
{
    class_<IntegratorWrapper, boost::noncopyable> integrator_wrapper_class("IntegratorInterface");
    integrator_wrapper_class.def("flow_bounds",(Pair<Float,IntervalVector>(IntegratorInterface::*)(const IntervalVectorFunction&,const IntervalVector&,const Float&)const)&IntegratorInterface::flow_bounds);
    integrator_wrapper_class.def("flow_step",(IntervalVectorFunctionModel(IntegratorInterface::*)(const IntervalVectorFunction&,const IntervalVector&,Float&)const)&IntegratorInterface::flow_step);
    integrator_wrapper_class.def("flow_step",(IntervalVectorFunctionModel(IntegratorInterface::*)(const IntervalVectorFunction&,const IntervalVector&,const Float&,const IntervalVector&)const)&IntegratorInterface::flow_step);
    integrator_wrapper_class.def("flow_to",(IntervalVectorFunctionModel(IntegratorInterface::*)(const IntervalVectorFunction&,const IntervalVector&,const Real&)const)&IntegratorInterface::flow_to);
    integrator_wrapper_class.def("flow",(List<IntervalVectorFunctionModel>(IntegratorInterface::*)(const IntervalVectorFunction&,const IntervalVector&,const Real&)const)&IntegratorInterface::flow);


    class_<TaylorPicardIntegrator > taylor_picard_integrator_class("TaylorPicardIntegrator",init<double>());
    taylor_picard_integrator_class.def("flow_bounds",(Pair<Float,IntervalVector>(IntegratorInterface::*)(const IntervalVectorFunction&,const IntervalVector&,const Float&)const)&TaylorPicardIntegrator::flow_bounds);
    taylor_picard_integrator_class.def("flow_step", (IntervalVectorFunctionModel(TaylorPicardIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,Float&)const)&TaylorPicardIntegrator::flow_step);
    taylor_picard_integrator_class.def("flow_step", (IntervalVectorFunctionModel(TaylorPicardIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,const Float&,const IntervalVector&)const)&TaylorPicardIntegrator::flow_step);
    taylor_picard_integrator_class.def("flow_to",(IntervalVectorFunctionModel(TaylorPicardIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,const Real&)const)&TaylorPicardIntegrator::flow_to);
    taylor_picard_integrator_class.def("flow",(List<IntervalVectorFunctionModel>(TaylorPicardIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,const Real&)const)&TaylorPicardIntegrator::flow);

    class_<TaylorSeriesIntegrator > taylor_series_integrator_class("TaylorSeriesIntegrator",init<double>());
    taylor_series_integrator_class.def("flow_bounds",(Pair<Float,IntervalVector>(IntegratorInterface::*)(const IntervalVectorFunction&,const IntervalVector&,const Float&)const)&TaylorSeriesIntegrator::flow_bounds);
    taylor_series_integrator_class.def("flow_step", (IntervalVectorFunctionModel(TaylorSeriesIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,Float&)const)&TaylorSeriesIntegrator::flow_step);
    taylor_series_integrator_class.def("flow_step", (IntervalVectorFunctionModel(TaylorSeriesIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,const Float&,const IntervalVector&)const)&TaylorSeriesIntegrator::flow_step);
    taylor_series_integrator_class.def("flow_to",(IntervalVectorFunctionModel(TaylorSeriesIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,const Real&)const)&TaylorSeriesIntegrator::flow_to);
    taylor_series_integrator_class.def("flow",(List<IntervalVectorFunctionModel>(TaylorSeriesIntegrator::*)(const IntervalVectorFunction&,const IntervalVector&,const Real&)const)&TaylorSeriesIntegrator::flow);

    class_<RungeKutta4Integrator > runge_kutta_4_integrator_class("RungeKutta4Integrator",init<double>());
    runge_kutta_4_integrator_class.def("step", &RungeKutta4Integrator::step);
    runge_kutta_4_integrator_class.def("evolve", &RungeKutta4Integrator::evolve);
}


void solver_submodule()
{
    to_python_list< Set< IntervalVector > >();
    to_python< List< IntervalVector > >();
    to_python< Pair< Float, IntervalVector > >();

    export_solver();
    export_integrator();

}


