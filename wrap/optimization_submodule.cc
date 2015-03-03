/***************************************************************************
 *            optimization_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "config.h"

#include <boost/python.hpp>

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "geometry/box.h"

#include "function/function.h"
#include "function/constraint.h"
#include "solvers/linear_programming.h"
#include "solvers/nonlinear_programming.h"
#include "solvers/constraint_solver.h"

using namespace boost::python;
using namespace Ariadne;


template<class X>
boost::python::tuple
python_compute_basis(const Matrix<X>& A) {
    Array<SizeType> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);
    boost::python::list l;
    for(SizeType i=0; i!=p.size(); ++i) {
        l.append(p[i]);
    }
    return boost::python::make_tuple(l,B);
}

template<class T> T get(const Array<T>& ary, SizeType i) { return ary[i]; }
template<class T> Void set(Array<T>& ary, SizeType i, const T& t) { ary[i]=t; }

template<class T>
Void export_internal_array(const char* name)
{
    class_< Array<T> > array_class(name,no_init);
    array_class.def("__len__", &Array<T>::size);
    array_class.def("__getitem__",&get<T>);
    array_class.def(boost::python::self_ns::str(self));
}


Void export_variable_type()
{
    typedef Array<Slackness> SlacknessArray;

    enum_<Slackness> variable_enum("Slackness");
    variable_enum.value("BASIS", BASIS);
    variable_enum.value("LOWER", LOWER);
    variable_enum.value("UPPER", UPPER);
}

Void export_constraint()
{
    class_<EffectiveConstraint> effective_nonlinear_constraint_class("EffectiveConstraint",init<Real,EffectiveScalarFunction,EffectiveNumericType>());
    effective_nonlinear_constraint_class.def(self_ns::str(self));

    class_<ValidatedConstraint> validated_nonlinear_constraint_class("ValidatedConstraint",init<ValidatedNumericType,ValidatedScalarFunction,ValidatedNumericType>());
    validated_nonlinear_constraint_class.def(init<ValidatedConstraint>());
    validated_nonlinear_constraint_class.def(init<EffectiveConstraint>());
    validated_nonlinear_constraint_class.def("lower_bound", &ValidatedConstraint::lower_bound, return_value_policy<copy_const_reference>());
    validated_nonlinear_constraint_class.def("upper_bound", &ValidatedConstraint::upper_bound, return_value_policy<copy_const_reference>());
    validated_nonlinear_constraint_class.def("function", (const ValidatedScalarFunction&(ValidatedConstraint::*)()const) &ValidatedConstraint::function, return_value_policy<copy_const_reference>());
    validated_nonlinear_constraint_class.def(self_ns::str(self));
}

Void export_interior_point_solver()
{
    to_python< Ariadne::Tuple< Vector<Float64>, Vector<Float64>, Vector<Float64> > >();

    class_<InteriorPointSolver> interior_point_solver_class("InteriorPointSolver",init<>());
    interior_point_solver_class.def("minimise", &InteriorPointSolver::minimise);
    interior_point_solver_class.def("feasible", (Kleenean(InteriorPointSolver::*)(const Vector<Float64>&,const Vector<Float64>&, const Matrix<Float64>&,const Vector<Float64>&)const) &InteriorPointSolver::feasible);
    interior_point_solver_class.def("validate_feasibility", &InteriorPointSolver::validate_feasibility);
}


Void export_constraint_solver()
{
    class_<ConstraintSolver> constraint_solver_class("ConstraintSolver", init<>());
    constraint_solver_class.def("hull_reduce", (Bool(ConstraintSolver::*)(UpperBox&,const ValidatedScalarFunctionInterface&,const ExactInterval&)const) &ConstraintSolver::hull_reduce);
    constraint_solver_class.def("box_reduce", (Bool(ConstraintSolver::*)(UpperBox&,const ValidatedScalarFunctionInterface&,const ExactInterval&,Nat)const) &ConstraintSolver::box_reduce);
    constraint_solver_class.def("monotone_reduce", (Bool(ConstraintSolver::*)(UpperBox&,const ValidatedScalarFunctionInterface&,const ExactInterval&,Nat)const) &ConstraintSolver::monotone_reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBox&,const List<ValidatedConstraint>&)const) &ConstraintSolver::reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBox&,const ValidatedVectorFunction&,const ExactBox&)const) &ConstraintSolver::reduce);
}



template<class X>
Void export_simplex_solver()
{
    typedef Array<SizeType> SizeArray;

    to_python< std::pair< Array<SizeType>, Matrix<X> > >();

    class_< SimplexSolver<X> > simplex_solver_class("SimplexSolver", init<>());
    simplex_solver_class.def("lpstep",(Bool(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,Array<Slackness>& ,SizeArray&,Matrix<X>&,Vector<X>&)const) &SimplexSolver<X>::lpstep);


    simplex_solver_class.def("feasible",(Kleenean(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&)const) &SimplexSolver<X>::feasible);

    simplex_solver_class.def("verify_feasibility",(Kleenean(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,const Array<Slackness>&)const) &SimplexSolver<X>::verify_feasibility);

    simplex_solver_class.def("compute_basis",(std::pair< SizeArray, Matrix<X> >(SimplexSolver<X>::*)(const Matrix<X>&)const) &SimplexSolver<X>::compute_basis);

}


Void optimization_submodule() {
    export_variable_type();
    export_constraint();
    export_array<SizeType>("SizeArray");
    export_internal_array<Slackness>("SlacknessArray");
    export_simplex_solver<Float64>();
    export_simplex_solver<Rational>();
    export_interior_point_solver();
    export_constraint_solver();
}
