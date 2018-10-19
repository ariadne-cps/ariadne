/***************************************************************************
 *            differential_inclusion.cpp
 *
 *  Copyright  2008-18  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

#include "inclusion_vector_field.hpp"
#include "../symbolic/expression_set.hpp"
#include "../function/symbolic_function.hpp"

namespace Ariadne {

BoxDomainType input_bounds_to_domain(RealVariableIntervals const& inputs) {
    List<IntervalDomainType> result;
    for (auto input : inputs) {
        result.push_back(cast_exact(IntervalDomainType(input.lower().get_d(),input.upper().get_d())));
    }
    return Vector<IntervalDomainType>(result);
}

Void InclusionVectorField::_acquire_properties() {

    List<Nat> input_indices;
    for (Nat i : range(this->dimension(),this->dimension()+this->number_of_inputs()))
        input_indices.append(i);

    const EffectiveVectorFormulaFunction* ff = dynamic_cast<const EffectiveVectorFormulaFunction*>(_function.raw_pointer());

    if (ff == nullptr)
        ARIADNE_THROW(NotFormulaFunctionException,"InclusionVectorField::_acquire_properties","The field function currently must be based on Formula");

    _is_input_affine = is_affine_in(*ff,input_indices);
    _is_input_additive = is_additive_in(*ff,input_indices);
}

InclusionVectorField::InclusionVectorField(DottedRealAssignments const& dynamics, RealVariableIntervals const& inputs)
{
    List<RealVariable> dyn_var_list = left_hand_sides(dynamics);
    List<RealVariable> inp_var_list;
    for (auto input : inputs) { inp_var_list.append(input.variable()); }

    Set<UntypedVariable> dyn_arg_set;
    for (auto expr : right_hand_sides(dynamics)) {
        dyn_arg_set.adjoin(expr.arguments());
    }

    Set<RealVariable> inp_var_set(inp_var_list);
    for (auto iv : inp_var_set) {
        if(not dyn_arg_set.contains(iv))
            ARIADNE_THROW(UnusedInputException,"InclusionVectorField(dynamics,inputs) with dynamics="<<dynamics<<" and inputs="<<inputs,"The input '" << iv << "' is not in the arguments of the dynamics");
    }

    Set<UntypedVariable> dyn_var_set(dyn_var_list);
    Set<UntypedVariable> inp_var_set_untyped(inp_var_set);
    for (auto av : dyn_arg_set) {
        if (not dyn_var_set.contains(av) and not inp_var_set_untyped.contains(av))
            ARIADNE_THROW(MissingInputException,"InclusionVectorField(dynamics,inputs) with dynamics="<<dynamics<<" and inputs="<<inputs,"The argument '" << av << "' of the dynamics is not among the inputs");
    }

    RealSpace var_spc(dyn_var_list);
    RealSpace inp_spc(inp_var_list);
    RealSpace spc = var_spc.adjoin(inp_spc);

    _function = make_function(spc,Vector<RealExpression>(right_hand_sides(dynamics)));
    _inputs = input_bounds_to_domain(inputs);
    _variable_names = variable_names(left_hand_sides(dynamics));

    _acquire_properties();
}

InclusionVectorField::InclusionVectorField(EffectiveVectorMultivariateFunction const& function, BoxDomainType const& inputs)
    : _function(function), _inputs(inputs)
{
    if (function.argument_size() != function.result_size()+inputs.size())
        ARIADNE_THROW(FunctionArgumentsMismatchException,"InclusionVectorField(function,inputs) with function="<<function<<" and inputs="<<inputs,
                "Incompatible inputs size (" << inputs.size() << ") with the provided function (R^" << function.argument_size() << " -> R^" << function.result_size() << ")");

    List<Identifier> variable_names;
    for (auto i : range(0,function.result_size()))
        variable_names.append(Identifier("x"+i));
    _variable_names = variable_names;

    _acquire_properties();
}


RealSpace InclusionVectorField::state_space() const
{
    return real_space(this->_variable_names);
}

}
