/***************************************************************************
 *            geometry/set.cpp
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

#include "../config.hpp"

#include "../utility/tuple.hpp"
#include "../numeric/numeric.hpp"
#include "../symbolic/templates.hpp"

#include "set.decl.hpp"
#include "set_interface.hpp"
#include "set.hpp"

#include "geometry.hpp"

namespace Ariadne {
#warning
//namespace Detail {

template<class OP, class... ARGS> struct LogicalExpression;

template<class OP, class ARG1, class ARG2> struct LogicalExpression<OP,ARG1,ARG2>
    : virtual LogicalInterface, Symbolic<OP,ARG1,ARG2>
{
    using Symbolic<OP,ARG1,ARG2>::Symbolic;
    virtual LogicalValue _check(Effort eff) const { return this->_op(this->_arg1,this->_arg2,eff).repr(); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,ARG1,ARG2>const&>(*this); }
};

//} // namespace Detail

template<class OP, class ARG1, class ARG2> decltype(auto) make_shared_logical_expression(OP const& op, ARG1 const& arg1, ARG2 const& arg2) {
    return std::make_shared<LogicalExpression<OP,ARG1,ARG2>>(op,arg1,arg2); }



template<class T> LowerKleenean EffectiveOvertSet<T>::overlaps(const EffectiveOpenSet<T>& other) const {
    return LowerKleenean(make_shared_logical_expression(Overlap(),*this,other));
}

template<class T> LowerKleenean EffectiveCompactSet<T>::inside(const EffectiveOpenSet<T>& other) const {
    return LowerKleenean(make_shared_logical_expression(Inside(),*this,other));
}

template<class T> LowerKleenean EffectiveCompactSet<T>::separated(const EffectiveClosedSet<T>& other) const {
    return LowerKleenean(make_shared_logical_expression(Separated(),*this,other));
}

template<class T> Kleenean EffectiveLocatedSet<T>::overlaps(const EffectiveRegularSet<T>& other) const {
    return Kleenean(make_shared_logical_expression(Overlap(),*this,other));
}

template<class T> Kleenean EffectiveLocatedSet<T>::inside(const EffectiveRegularSet<T>& other) const {
    return Kleenean(make_shared_logical_expression(Inside(),*this,other));
}

template<class T> Kleenean EffectiveLocatedSet<T>::separated(const EffectiveRegularSet<T>& other) const {
    return Kleenean(make_shared_logical_expression(Separated(),*this,other));
}


template<class T> ValidatedLowerKleenean intersection_overlap(const OvertSetInterface<EffectiveTag,T>& ovs, const OpenSetInterface<EffectiveTag,T>& ops, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedLowerKleenean intersection_inside(const CompactSetInterface<EffectiveTag,T>& cps, const OpenSetInterface<EffectiveTag,T>& ops, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedLowerKleenean intersection_separated(const CompactSetInterface<EffectiveTag,T>& cps, const ClosedSetInterface<EffectiveTag,T>& cls, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedKleenean intersection_overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedKleenean intersection_inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedKleenean intersection_separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedKleenean intersection_overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Accuracy eps);
template<class T> ValidatedKleenean intersection_inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Accuracy eps);
template<class T> ValidatedKleenean intersection_separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Accuracy eps);





template<class T> ValidatedLowerKleenean OvertSetInterface<EffectiveTag,T>::_overlaps(BasicSetType const& bx, Effort eff) const {
    return this->_overlaps(bx);
}


template<class T> LowerKleenean OvertSetInterface<EffectiveTag,T>::overlaps(OpenSetInterface<P,T> const& other) const {
    return OvertSet<P,T>(this->_copy()).overlaps(OpenSet<P,T>(other._copy()));
}

template<class T> LowerKleenean CompactSetInterface<EffectiveTag,T>::inside(OpenSetInterface<P,T> const& other) const {
    return CompactSet<P,T>(this->_copy()).inside(OpenSet<P,T>(other._copy()));
}

template<class T> LowerKleenean CompactSetInterface<EffectiveTag,T>::separated(ClosedSetInterface<P,T> const& other) const {
    return CompactSet<P,T>(this->_copy()).separated(ClosedSet<P,T>(other._copy()));
}

template<class T> Kleenean LocatedSetInterface<EffectiveTag,T>::overlaps(RegularSetInterface<P,T> const& other) const {
    return LocatedSet<P,T>(this->_copy()).overlaps(RegularSet<P,T>(other._copy()));
}

template<class T> Kleenean LocatedSetInterface<EffectiveTag,T>::inside(RegularSetInterface<P,T> const& other) const {
    return LocatedSet<P,T>(this->_copy()).inside(RegularSet<P,T>(other._copy()));
}

template<class T> Kleenean LocatedSetInterface<EffectiveTag,T>::separated(RegularSetInterface<P,T> const& other) const {
    return LocatedSet<P,T>(this->_copy()).separated(RegularSet<P,T>(other._copy()));
}






template<class T> ValidatedLowerKleenean OvertSetInterface<ValidatedTag,T>::overlaps(OpenSetInterface<P,T> const& other) const {
    return SetOperations<T>::overlap(*this,other);
}

template<class T> ValidatedLowerKleenean CompactSetInterface<ValidatedTag,T>::inside(OpenSetInterface<P,T> const& other) const {
    return SetOperations<T>::inside(*this,other);
}

template<class T> ValidatedLowerKleenean CompactSetInterface<ValidatedTag,T>::separated(ClosedSetInterface<P,T> const& other) const {
    return SetOperations<T>::overlap(*this,other);
}


template<class T> ValidatedKleenean LocatedSetInterface<ValidatedTag,T>::overlaps(RegularSetInterface<P,T> const& other) const {
    return Ariadne::overlap(*this,other);
}

template<class T> ValidatedKleenean LocatedSetInterface<ValidatedTag,T>::inside(RegularSetInterface<P,T> const& other) const {
    return SetOperations<T>::inside(*this,other);
}

template<class T> ValidatedKleenean LocatedSetInterface<ValidatedTag,T>::separated(RegularSetInterface<P,T> const& other) const {
    return SetOperations<T>::separated(*this,other);
}




template class BoundedSetInterface<EffectiveTag,RealVector>;
template class OvertSetInterface<EffectiveTag,RealVector>;
template class OpenSetInterface<EffectiveTag,RealVector>;
template class ClosedSetInterface<EffectiveTag,RealVector>;
template class CompactSetInterface<EffectiveTag,RealVector>;
template class RegularSetInterface<EffectiveTag,RealVector>;
template class LocatedSetInterface<EffectiveTag,RealVector>;
template class RegularLocatedSetInterface<EffectiveTag,RealVector>;

template class BoundedSet<EffectiveTag,RealVector>;
template class OvertSet<EffectiveTag,RealVector>;
template class OpenSet<EffectiveTag,RealVector>;
template class ClosedSet<EffectiveTag,RealVector>;
template class CompactSet<EffectiveTag,RealVector>;
template class RegularSet<EffectiveTag,RealVector>;
template class LocatedSet<EffectiveTag,RealVector>;
template class RegularLocatedSet<EffectiveTag,RealVector>;

template class BoundedSet<ValidatedTag,RealVector>;
template class OvertSet<ValidatedTag,RealVector>;
template class OpenSet<ValidatedTag,RealVector>;
template class ClosedSet<ValidatedTag,RealVector>;
template class CompactSet<ValidatedTag,RealVector>;
template class RegularSet<ValidatedTag,RealVector>;
template class LocatedSet<ValidatedTag,RealVector>;
template class RegularLocatedSet<ValidatedTag,RealVector>;




} // namespace Ariadne

