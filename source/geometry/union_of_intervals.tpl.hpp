/***************************************************************************
 *            union_of_intervals.tpl.hpp
 *
 *  Copyright  2020-21  Pieter Collins
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

#include "union_of_intervals.hpp"

namespace Ariadne {

namespace {
inline Bool decide_cast_bool(Boolean const& bl) {
    return bl; }
inline Bool decide_cast_bool(ValidatedLowerKleenean const& lk) {
    return definitely(lk); }
inline Bool decide_cast_bool(ValidatedUpperKleenean const& uk) {
    return possibly(uk); }
inline Bool decide_cast_bool(ApproximateKleenean const& ak) {
    return probably(ak); }
}


template<class U> UnionOfIntervals<U> UnionOfIntervals<U>::_intersection(UnionOfIntervals<U> const& ivls1, UnionOfIntervals<U> const& ivls2) {
    // TODO: Make more efficient
    List<Interval<U>> rivls;
    auto iter1=ivls1.begin();
    auto iter2=ivls2.begin();
    while (iter1!=ivls1.end() && iter2!=ivls2.end()) {
        if (decide_cast_bool(iter1->lower_bound() >= iter2->upper_bound())) {
            ++iter1;
        } else if (decide_cast_bool(iter2->lower_bound() >= iter1->upper_bound())) {
            ++iter2;
        } else {
            rivls.append(Interval<U>(max(iter1->lower_bound(),iter2->lower_bound()),min(iter1->upper_bound(),iter2->upper_bound())));
            if (decide_cast_bool(iter1->upper_bound()<iter2->upper_bound())) { ++iter1; } else { ++iter2; }
        }
    }
    return UnionOfIntervals<U>(rivls);
}

template<class U> UnionOfIntervals<U> UnionOfIntervals<U>::_join(UnionOfIntervals<U> const& ivls1, UnionOfIntervals<U> const& ivls2) {
    // TODO: Make more efficient
    return UnionOfIntervals<U>(catenate(ivls1._data,ivls2._data));
}

template<class U> Bool UnionOfIntervals<U>::_intersect(UnionOfIntervals<U> const& ivls1, UnionOfIntervals<U> const& ivls2) {
    auto iter1=ivls1.begin();
    auto iter2=ivls2.begin();
    while (iter1!=ivls1.end() && iter2!=ivls2.end()) {
        if (decide_cast_bool(iter1->lower_bound() >= iter2->upper_bound())) {
            ++iter2;
        } else if (decide_cast_bool(iter1->upper_bound() <= iter2->lower_bound())) {
            ++iter1;
        } else {
            return true;
        }
    }
    return false;
}

template<class U> Bool UnionOfIntervals<U>::_subset(UnionOfIntervals<U> const& ivls1, UnionOfIntervals<U> const& ivls2) {
    auto iter1=ivls1.begin();
    auto iter2=ivls2.begin();
    while (iter1!=ivls1.end() && iter2!=ivls2.end()) {
        if (decide_cast_bool(iter1->lower_bound() >= iter2->upper_bound())) {
            ++iter2;
        } else if (decide_cast_bool(iter1->lower_bound() >= iter2->lower_bound() && iter1->upper_bound() <= iter2->upper_bound())) {
            ++iter1;
        } else {
            return false;
        }
    }
    if (iter1==ivls1.end()) {
        return true;
    } else {
        return false;
    }
}


template<class U> auto UnionOfIntervals<U>::measure() const -> MeasureType {
    MeasureType r(nul(_data[0].width()));
    for (auto ivl : this->_data) {
        r+=ivl.width();
    }
    return r;
}

template<class U> Interval<U> UnionOfIntervals<U>::bounding_box() const {
    ARIADNE_ASSERT(!_data.empty());
    return Interval<U>(this->_data.front().lower_bound(),this->_data.back().upper_bound());
}

template<class U> Void UnionOfIntervals<U>::_sort() {
    auto cmp = [](Interval<U>const& ivl1, Interval<U> const& ivl2) {
        return decide_cast_bool(ivl1.lower_bound()<ivl2.lower_bound() or (ivl1.lower_bound()==ivl2.lower_bound() and ivl1.upper_bound()<ivl2.upper_bound())); };
    std::sort(this->_data.begin(),this->_data.end(),cmp);
}




template<class U> Void UnionOfIntervals<U>::_combine() {
    ConstIterator start=begin();
    Iterator curr=_data.begin();
    ConstIterator next=curr;
    while (next!=end()) {
        *curr=*next;
        ++next;
        while (next!=end() && decide_cast_bool(next->lower_bound() <= curr->upper_bound())) {
            curr->set_upper_bound(max(curr->upper_bound(),next->upper_bound()));
            ++next;
        }
        ++curr;
    }

    _data.resize(static_cast<SizeType>(curr-start));
}

} // namespace Ariadne
