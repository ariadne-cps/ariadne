#include "system/system_function.hpp"
#include "system/hybrid_automaton_interface.hpp"
#include "function/functional.hpp"
#include "utility/container.hpp"
#include "symbolic/expression.hpp"

namespace Ariadne {

List<RealAssignment>
algebraic_sort(const List<RealAssignment>& auxiliary) {
    List<RealVariable> lhs_array=left_hand_sides(auxiliary);
    LinkedList<RealVariable> lhs_list(lhs_array.begin(),lhs_array.end());
    Map<Identifier, Set<UntypedVariable> > dependencies;
    Set<UntypedVariable> variables(lhs_array.begin(),lhs_array.end());

    for(List<RealAssignment>::ConstIterator asgn_iter=auxiliary.begin();
        asgn_iter!=auxiliary.end(); ++asgn_iter)
    {
        dependencies.insert(asgn_iter->lhs.base().name(), intersection(asgn_iter->rhs.arguments(),variables));
    }

    List<RealAssignment> sorted_auxiliary;
    sorted_auxiliary.reserve(auxiliary.size());
    if(!lhs_list.empty()) {
        Bool found=false;
        for(LinkedList<RealVariable>::Iterator iter=lhs_list.begin(); iter!=lhs_list.end(); ) {
            if(dependencies[iter->name()].empty()) {
                for(Map<Identifier, Set<UntypedVariable> >::Iterator dep_iter=dependencies.begin(); dep_iter!=dependencies.end(); ++dep_iter) {
                    dep_iter->second.remove(static_cast<UntypedVariable>(*iter));
                }
                for(List<RealAssignment>::ConstIterator asgn_iter=auxiliary.begin();
                    asgn_iter!=auxiliary.end(); ++asgn_iter)
                {
                    if(asgn_iter->lhs.name()==iter->name()) {
                        sorted_auxiliary.append(*asgn_iter);
                        break;
                    }
                }
                dependencies.erase(iter->name());
                found=true;
                lhs_list.erase(iter);
                iter=lhs_list.begin();
            } else {
                ++iter;
            }
        }
        if(!found) {
            ARIADNE_THROW(AlgebraicLoopError,"HybridAutomaton::sort(List<RealAssignment>)",
                          "Algebraic dependencies among variables "<<lhs_list<<" in auxiliary equations "<<auxiliary);
        }
    }
    ARIADNE_ASSERT(auxiliary.size()==sorted_auxiliary.size());
    return sorted_auxiliary;
}

EffectiveVectorMultivariateFunction make_auxiliary_function(
    Space<Real> const& state_space,
    List<RealAssignment> const& algebraic)
{
    List<RealAssignment> sorted_algebraic = algebraic_sort(algebraic);
    Space<Real> auxiliary_space(left_hand_sides(algebraic));
    Space<Real> sorted_auxiliary_space(left_hand_sides(sorted_algebraic));

    RealExpression default_expression;
    Vector<RealExpression> results(sorted_algebraic.size(),default_expression);
    for(SizeType j=0; j!=sorted_algebraic.size(); ++j) {
        SizeType i=auxiliary_space[sorted_auxiliary_space[j]];
        results[i]=substitute(sorted_algebraic[j].rhs,sorted_algebraic);
    }
    return make_function(state_space,results);
}

EffectiveVectorMultivariateFunction make_dynamic_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<DottedRealAssignment> const& differential)
{
    RealExpression default_expression;
    Vector<RealExpression> results(differential.size(),default_expression);
    for(SizeType i=0; i!=differential.size(); ++i) { results[space.index(differential[i].lhs.base())]=substitute(differential[i].rhs,algebraic); }

    return make_function(space,results);
}

EffectiveVectorMultivariateFunction make_reset_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<PrimedRealAssignment> const& primed)
{
    RealExpression default_expression;
    Vector<RealExpression> results(primed.size(),default_expression);
    for(SizeType i=0; i!=primed.size(); ++i) { results[i]=substitute(primed[i].rhs,algebraic); }

    return make_function(space,results);
}

EffectiveScalarMultivariateFunction make_constraint_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    ContinuousPredicate const& constraint,
    Sign sign)
{
    RealExpression constraint_expression=indicator(constraint,sign);
    return make_function(space,substitute(constraint_expression,algebraic));
}

} // namespace Ariadne
