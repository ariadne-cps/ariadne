/***************************************************************************
 *            optimization_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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


template<class X>
pybind11::tuple
python_compute_basis(const Matrix<X>& A) {
    Array<SizeType> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);
    pybind11::list l;
    for(SizeType i=0; i!=p.size(); ++i) {
        l.append(p[i]);
    }
    return pybind11::make_tuple(l,B);
}

template<class T> T get(const Array<T>& ary, SizeType i) { return ary[i]; }
template<class T> Void set(Array<T>& ary, SizeType i, const T& t) { ary[i]=t; }

template<class T>
Void export_internal_array(pybind11::module& module, const char* name)
{
    pybind11::class_< Array<T> > array_class(module,name);
    array_class.def("__len__", &Array<T>::size);
    array_class.def("__getitem__",&get<T>);
    array_class.def("__str__", &__cstr__<Array<T>>);
}


Void export_variable_type(pybind11::module& module)
{
    pybind11::enum_<Slackness> variable_enum(module,"Slackness");
    variable_enum.value("BASIS", Slackness::BASIS);
    variable_enum.value("LOWER", Slackness::LOWER);
    variable_enum.value("UPPER", Slackness::UPPER);
}

Void export_constraint(pybind11::module& module)
{
    pybind11::class_<EffectiveConstraint> effective_nonlinear_constraint_class(module,"EffectiveConstraint");
    effective_nonlinear_constraint_class.def(pybind11::init<Real,EffectiveScalarFunction,EffectiveNumericType>());
    effective_nonlinear_constraint_class.def("__str__",&__cstr__<EffectiveConstraint>);
    
    pybind11::class_<ValidatedConstraint> validated_nonlinear_constraint_class(module,"ValidatedConstraint");
    validated_nonlinear_constraint_class.def(pybind11::init<ValidatedNumericType,ValidatedScalarFunction,ValidatedNumericType>());
    validated_nonlinear_constraint_class.def(pybind11::init<ValidatedConstraint>());
    validated_nonlinear_constraint_class.def(pybind11::init<EffectiveConstraint>());
    validated_nonlinear_constraint_class.def("lower_bound", &ValidatedConstraint::lower_bound);
    validated_nonlinear_constraint_class.def("upper_bound", &ValidatedConstraint::upper_bound);
    validated_nonlinear_constraint_class.def("function", (const ValidatedScalarFunction&(ValidatedConstraint::*)()const) &ValidatedConstraint::function);
    validated_nonlinear_constraint_class.def("__str__",&__cstr__<ValidatedConstraint>);
}

Void export_interior_point_solver(pybind11::module& module)
{
//    to_python< Ariadne::Tuple< Vector<FloatDP>, Vector<FloatDP>, Vector<FloatDP> > >();

    pybind11::class_<InteriorPointSolver> interior_point_solver_class(module,"InteriorPointSolver");
    interior_point_solver_class.def(pybind11::init<>());
    interior_point_solver_class.def("minimise", &InteriorPointSolver::minimise);
    interior_point_solver_class.def("feasible", (ValidatedKleenean(InteriorPointSolver::*)(const Vector<FloatDP>&,const Vector<FloatDP>&, const Matrix<FloatDP>&,const Vector<FloatDP>&)const) &InteriorPointSolver::feasible);
    interior_point_solver_class.def("validate_feasibility", &InteriorPointSolver::validate_feasibility);
}


Void export_constraint_solver(pybind11::module& module)
{
    pybind11::class_<ConstraintSolver> constraint_solver_class(module,"ConstraintSolver");
    constraint_solver_class.def(pybind11::init<>());
    constraint_solver_class.def("hull_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarFunction&,const ExactIntervalType&)const) &ConstraintSolver::hull_reduce);
    constraint_solver_class.def("box_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarFunction&,const ExactIntervalType&,Nat)const) &ConstraintSolver::box_reduce);
    constraint_solver_class.def("monotone_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarFunction&,const ExactIntervalType&,Nat)const) &ConstraintSolver::monotone_reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const List<ValidatedConstraint>&)const) &ConstraintSolver::reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedVectorFunction&,const ExactBoxType&)const) &ConstraintSolver::reduce);
}


template<class X>
Void export_simplex_solver(pybind11::module& module)
{
    typedef Array<SizeType> SizeArray;

//    to_python< std::pair< Array<SizeType>, Matrix<X> > >();

    pybind11::class_< SimplexSolver<X> > simplex_solver_class(module,(class_name<X>()+"SimplexSolver").c_str());
    simplex_solver_class.def(pybind11::init<>());
    simplex_solver_class.def("lpstep",(Bool(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,Array<Slackness>& ,SizeArray&,Matrix<X>&,Vector<X>&)const) &SimplexSolver<X>::lpstep);


    simplex_solver_class.def("feasible",(ValidatedKleenean(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&)const) &SimplexSolver<X>::feasible);

    simplex_solver_class.def("verify_feasibility",(ValidatedKleenean(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,const Array<Slackness>&)const) &SimplexSolver<X>::verify_feasibility);

    simplex_solver_class.def("compute_basis",(std::pair< SizeArray, Matrix<X> >(SimplexSolver<X>::*)(const Matrix<X>&)const) &SimplexSolver<X>::compute_basis);

}


Void optimization_submodule(pybind11::module& module) {
    export_variable_type(module);
    export_constraint(module);
    export_array<SizeType>(module,"SizeArray");
    export_internal_array<Slackness>(module,"SlacknessArray");
    export_simplex_solver<FloatDP>(module);
    export_simplex_solver<Rational>(module);
    export_interior_point_solver(module);
    export_constraint_solver(module);
}
