/***************************************************************************
 *            regular_expression.tpl.hpp
 *
 *  Copyright  2018  Pieter Collins
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

#include "regular_expression.hpp"

namespace Ariadne {

template<class T> RegularExpression<T> simplify(RegularExpression<T> re) {
    if (re.template holds_alternative<Catenation<T>>()) {
        Catenation<T> const& cat = re.template extract_alternative<Catenation<T>>();
        Catenation<T> res({});
        for(auto elmt : cat.elements()) {
            elmt = simplify(elmt);
            if (elmt.template holds_alternative<Empty<T>>()) { return Empty<T>(); }
            res.append(elmt);
        }
        if (res.elements().size()==1) { return res.elements().front(); }
        return res;
    } else if (re.template holds_alternative<Alternation<T>>()) {
        Alternation<T> const& alt = re.template extract_alternative<Alternation<T>>();
        Alternation<T> res({});
        for(auto poss : alt.possibilities()) {
            poss = simplify(poss);
            if (not poss.template holds_alternative<Empty<T>>()) { res.adjoin(poss); }
        }
        if (res.possibilities().size()==1) { return res.possibilities().front(); }
        return re;
        //return RegularExpression<T>( res );
    } else {
        return re;
    }
}

} // namespace Ariadne

