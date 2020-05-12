/***************************************************************************
 *            geometry/geometry.cpp
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
#include "geometry.hpp"
#include "set_interface.hpp"
#include "set.hpp"

namespace Ariadne {

template<class T> ValidatedLowerKleenean intersection_overlap(const OvertSetInterface<EffectiveTag,T>& ovs, const OpenSetInterface<EffectiveTag,T>& ops, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedLowerKleenean intersection_inside(const CompactSetInterface<EffectiveTag,T>& cps, const OpenSetInterface<EffectiveTag,T>& ops, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedLowerKleenean intersection_separated(const CompactSetInterface<EffectiveTag,T>& cps, const ClosedSetInterface<EffectiveTag,T>& cls, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedKleenean intersection_overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedKleenean intersection_inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff);
template<class T> ValidatedKleenean intersection_separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff);




#warning
/*
template<class T> LowerKleenean overlap(const EffectiveOvertSetInterface<T>& ovs, const EffectiveOpenSetInterface<T>& ops) {
    return LowerKleenean(make_shared_logical_expression(Overlap(),EffectiveOvertSet<T>(ovs),EffectiveOpenSet<T>(ops)));
}

template<class T> LowerKleenean inside(const EffectiveCompactSetInterface<T>& cps, const EffectiveOpenSetInterface<T>& ops) {
    return LowerKleenean(make_shared_logical_expression(Inside(),EffectiveCompactSet<T>(cps),EffectiveOpenSet<T>(ops)));
}

template<class T> LowerKleenean separated(const EffectiveCompactSetInterface<T>& cps, const EffectiveClosedSetInterface<T>& cls) {
    return LowerKleenean(make_shared_logical_expression(Separated(),EffectiveCompactSet<T>(cps),EffectiveClosedSet<T>(cls)));
}

template<class T> Kleenean overlap(const EffectiveLocatedSetInterface<T>& ls, const EffectiveRegularSetInterface<T>& rs) {
    return Kleenean(make_shared_logical_expression(Overlap(),EffectiveLocatedSet<T>(ls),EffectiveRegularSet<T>(rs)));
}

template<class T> Kleenean inside(const EffectiveLocatedSetInterface<T>& ls, const EffectiveRegularSetInterface<T>& rs) {
    return Kleenean(make_shared_logical_expression(Inside(),EffectiveLocatedSet<T>(ls),EffectiveRegularSet<T>(rs)));
}

template<class T> Kleenean separated(const EffectiveLocatedSetInterface<T>& ls, const EffectiveRegularSetInterface<T>& rs) {
    return Kleenean(make_shared_logical_expression(Separated(),EffectiveLocatedSet<T>(ls),EffectiveRegularSet<T>(rs)));
}
*/


inline FloatDPExactBox cast_basic_set(Vector<EffectiveNumber> const& pt, Effort eff) {
    return cast_exact_box(
        FloatDPUpperBox(pt.size(),[&](SizeType i){return FloatDPUpperInterval(pt[i].get(BoundedTag(),double_precision));}));
}
inline Box<Interval<FloatDPValue>> cast_basic_set(Vector<EffectiveNumber> const& pt, Accuracy acc) {
    return cast_exact_box(
        FloatDPUpperBox(pt.size(),[&](SizeType i){return FloatDPUpperInterval(pt[i].get(BoundedTag(),double_precision));}));
}
/*
inline Box<Interval<FloatMPValue>> const& cast_exact_box(Box<Interval<FloatMPUpperBound>> const& bx) {
    return reinterpret_cast<Box<Interval<FloatMPValue>>const&>(bx);
}
inline Box<Interval<FloatMPValue>> cast_basic_set(Vector<EffectiveNumber> const& pt, Accuracy acc) {
    MultiplePrecision prec(acc.bits());
    return cast_exact_box(
        Box<FloatMPUpperInterval>(pt.size(),[&](SizeType i){return FloatMPUpperInterval(FloatMPBounds(pt[i].get(prec)));}));
}
*/

template<class T> ValidatedLowerKleenean SetOperations<T>::contains(EffectiveOpenSetInterface<T> const& ops, ElementType const& pt, Effort eff) {
    return ops.covers(cast_basic_set(pt,eff));
}

template<class T> ValidatedUpperKleenean SetOperations<T>::contains(EffectiveClosedSetInterface<T> const& cls, ElementType const& pt, Effort eff) {
    return not cls.separated(cast_basic_set(pt,eff));
}

template<class T> ValidatedKleenean SetOperations<T>::contains(EffectiveRegularSetInterface<T> const& rs, ElementType const& pt, Effort eff) {
    auto bs=cast_basic_set(pt,eff);
    if (definitely(rs.covers(bs))) { return true; }
    else if (definitely(rs.separated(bs))) { return false; }
    else { return indeterminate; }
}

template<class T> ValidatedKleenean SetOperations<T>::contains(EffectiveRegularSetInterface<T> const& rs, ElementType const& pt, Accuracy acc) {
    auto bs=cast_basic_set(pt,acc);
    if (definitely(rs.covers(bs))) { return true; }
    else if (definitely(rs.separated(bs))) { return false; }
    else { return indeterminate; }
}


template<class T> ValidatedLowerKleenean SetOperations<T>::overlap(OvertSetInterface<EffectiveTag,T> const& ovs, OpenSetInterface<EffectiveTag,T> const& ops, Effort eff) {
    FloatDPValue val(pow(two,eff.work()),double_precision);
    ExactBoxType box(ovs.dimension(),ExactIntervalType(-val,+val));
    return intersection_overlap(ovs,ops,box,eff);
}

template<class T> ValidatedLowerKleenean SetOperations<T>::inside(CompactSetInterface<EffectiveTag,T> const& cps, OpenSetInterface<EffectiveTag,T> const& ops, Effort eff) {
    return intersection_inside(cps,ops,cast_exact_box(cps.bounding_box()),eff);
}

template<class T> ValidatedLowerKleenean SetOperations<T>::separated(CompactSetInterface<EffectiveTag,T> const& cps, ClosedSetInterface<EffectiveTag,T> const& cls, Effort eff) {
    return intersection_separated(cps,cls,cast_exact_box(cps.bounding_box()),eff);
}

template<class T> ValidatedKleenean SetOperations<T>::overlap(LocatedSetInterface<EffectiveTag,T> const& self, RegularSetInterface<EffectiveTag,T> const& other, Effort eff) {
    return intersection_overlap(self,other,cast_exact_box(self.bounding_box()),eff);
}

template<class T> ValidatedKleenean SetOperations<T>::inside(LocatedSetInterface<EffectiveTag,T> const& self, RegularSetInterface<EffectiveTag,T> const& other, Effort eff) {
    return intersection_inside(self,other,cast_exact_box(self.bounding_box()),eff);
}

template<class T> ValidatedKleenean SetOperations<T>::separated(LocatedSetInterface<EffectiveTag,T> const& self, RegularSetInterface<EffectiveTag,T> const& other, Effort eff) {
    return intersection_separated(self,other,cast_exact_box(self.bounding_box()),eff);
}

template<class T> ValidatedKleenean
SetOperations<T>::overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, Accuracy acc)
{
    ExactBoxType bb=cast_exact_box(ls.bounding_box());
    if(definitely(bb.is_empty())) { return false; }
    return intersection_overlap(ls,rs,bb,acc);
}


template<class T> ValidatedKleenean
SetOperations<T>::inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, Accuracy acc)
{
    ExactBoxType bb=cast_exact_box(ls.bounding_box());
    if(definitely(bb.is_empty())) { return true; }
    return intersection_inside(ls,rs,bb,acc);
}

template<class T> ValidatedKleenean
SetOperations<T>::separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, Accuracy acc)
{
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    BasicSetType bb=cast_exact_box(ls.bounding_box());
    if(definitely(bb.is_empty())) { return true; }
    return intersection_separated(ls,rs,bb,acc);
}




Int log10floor(Dyadic);
#warning Check BoxSet

template<class T> ValidatedKleenean
intersection_overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff)
{
    return intersection_overlap(ls,rs,bx,Accuracy(two^(-eff.work())));
}


template<class T> ValidatedKleenean
intersection_inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff)
{
    return intersection_inside(ls,rs,bx,Accuracy(two^(-eff.work())));
}

template<class T> ValidatedKleenean
intersection_separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Effort eff)
{
    return intersection_separated(ls,rs,bx,Accuracy(two^(-eff.work())));
}


template<class T> ValidatedKleenean
intersection_overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Accuracy eps)
{
    if(definitely(ls.separated(bx))) {
        return false;
    }
    if(definitely(rs.separated(bx))) {
        return false;
    }
    else if(definitely(rs.covers(bx))) {
        return true;
    }
    else if(definitely(bx.radius()<eps.error())) {
        return indeterminate;
    } else {
        ExactBoxType bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(definitely(ls.separated(bx1))) {
            return intersection_overlap(ls,rs,bx2,eps);
        } else if(definitely(ls.separated(bx2))) {
            return intersection_overlap(ls,rs,bx1,eps);
        } else {
            return intersection_overlap(ls,rs,bx1,eps) || intersection_overlap(ls,rs,bx2,eps);
        }
    }
}

template<class T> ValidatedKleenean
intersection_inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Accuracy eps)
{
    if(definitely(ls.separated(bx) || rs.separated(bx))) {
        return true;
    } else if(decide(bx.radius()<eps.error())) {
        return indeterminate;
    } else {
        ExactBoxType bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(definitely(ls.separated(bx1))) {
            return intersection_inside(ls,rs,bx2,eps);
        } else if(definitely(ls.separated(bx2))) {
            return intersection_inside(ls,rs,bx1,eps);
        } else {
            return intersection_inside(ls,rs,bx1,eps) && intersection_inside(ls,rs,bx2,eps);
        }
    }
}


template<class T> ValidatedKleenean
intersection_separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const ExactBoxType& bx, Accuracy eps)
{
    if(definitely(ls.separated(bx) || rs.separated(bx))) {
        return true;
    } else if(definitely(bx.radius()<eps.error())) {
        return indeterminate;
    } else {
        ExactBoxType bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(definitely(ls.separated(bx1))) {
            return intersection_separated(ls,rs,bx2,eps);
        } else if(definitely(ls.separated(bx2))) {
            return intersection_separated(ls,rs,bx1,eps);
        } else {
            return intersection_separated(ls,rs,bx1,eps) && intersection_separated(ls,rs,bx2,eps);
        }
    }
}


#warning Should not need cast of interminate to ValidatedLowerKleenean
template<class T> ValidatedLowerKleenean
intersection_overlap(const OvertSetInterface<EffectiveTag,T>& ovs, const OpenSetInterface<EffectiveTag,T>& ops, const ExactBoxType& bx, Effort eff)
{
    Accuracy eps(pow(two,-eff.work()));
    if(definitely(ovs.overlaps(bx,eff))) {
        if(definitely(ops.covers(bx))) {
            return true;
        } else if(decide(bx.radius()<eps.error())) {
            return ValidatedKleenean(indeterminate);
        } else {
            ExactBoxType bx1,bx2;
            make_lpair(bx1,bx2)=split(bx);
            if(definitely(intersection_overlap(ovs,ops,bx1,eff))) {
                return true;
            } else {
                return intersection_overlap(ovs,ops,bx2,eff);
            }
        }
    } else {
        return ValidatedKleenean(indeterminate);
    }
}

template<class T> ValidatedLowerKleenean
intersection_inside(const CompactSetInterface<EffectiveTag,T>& cps, const OpenSetInterface<EffectiveTag,T>& ops, const ExactBoxType& bx, Effort eff)
{
    Accuracy eps(pow(two,-eff.work()));
#warning
//    if(definitely(cps.separated(bx,eff) || ops.covers(bx,eff))) {
    if(definitely(cps.separated(bx) || ops.covers(bx))) {
        return true;
    } else if(decide(bx.radius()<eps.error())) {
        return ValidatedKleenean(indeterminate);
    } else {
        ExactBoxType bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(definitely(intersection_inside(cps,ops,bx1,eff))) {
            return intersection_inside(cps,ops,bx2,eff);
        } else {
            return ValidatedKleenean(indeterminate);
        }
    }
}


template<class T> ValidatedLowerKleenean
intersection_separated(const CompactSetInterface<EffectiveTag,T>& cps1, const ClosedSetInterface<EffectiveTag,T>& cls2, const ExactBoxType& bx, Effort eff)
{
    Accuracy eps(pow(two,-eff.work()));
    if(definitely(cps1.separated(bx) || cls2.separated(bx))) {
        return true;
    } else if(decide(bx.radius()<eps.error())) {
        return ValidatedKleenean(indeterminate);
    } else {
        ExactBoxType bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(definitely(intersection_separated(cps1,cls2,bx1,eff))) {
            return intersection_separated(cps1,cls2,bx2,eff);
        } else {
            return ValidatedKleenean(indeterminate);
        }
    }
}


template<class T> ValidatedLowerKleenean
SetOperations<T>::overlap(const OvertSetInterface<ValidatedTag,T>& ovs, const OpenSetInterface<ValidatedTag,T>& ops) {
    ARIADNE_NOT_IMPLEMENTED;
}


template class SetOperations<RealVector>;


} // namespace Ariadne

