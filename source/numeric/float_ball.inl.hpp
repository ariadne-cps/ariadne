/***************************************************************************
 *            float_ball.tpl.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

namespace Ariadne {


template<class F, class FE> Ball<F,FE> Operations<Ball<F,FE>>::_rec(Ball<F,FE> const& x) {
    static const bool use_midpoint = true;
    if constexpr (use_midpoint) {// Use this code to find value same as reciprocal value
        auto rv=rec(approx,x._v);
        auto ru=rec(up,sub(down,x._v,x._e));
        auto rl=rec(down,add(up,x._v,x._e));
        auto re=max(sub(up,ru,rv),sub(up,rv,rl));
        return Ball<F,FE>(rv,_make_error<FE>(re));
    } else {
        // Use this code to get same result as interval computation
        auto ru=rec(up,sub(down,x._v,x._e));
        auto rl=rec(down,add(up,x._v,x._e));
        auto re=hlf(sub(up,ru,rl));
        auto rv=hlf(add(near,rl,ru));
        return Ball<F,FE>(rv,re);
    }
}

template<class F, class FE> Ball<F,FE> Operations<Ball<F,FE>>::_add(Ball<F,FE> const& x, Ball<F,FE> const& y) {
    auto rv=add(near,x._v,y._v);
    auto ru=add(up,x._v,y._v);
    auto rl=add(down,x._v,y._v);
    auto ae=_make_error<FE>(hlf(sub(up,ru,rl)));
    auto re=add(up,ae,add(up,x._e,y._e));
    return Ball<F,FE>(rv,re);
}

template<class F, class FE> Ball<F,FE> Operations<Ball<F,FE>>::_sub(Ball<F,FE> const& x, Ball<F,FE> const& y) {
    auto rv=sub(near,x._v,y._v);
    auto ru=sub(up,x._v,y._v);
    auto rl=sub(down,x._v,y._v);
    auto ae=_make_error<FE>(hlf(sub(up,ru,rl)));
    auto re=add(up,ae,add(up,x._e,y._e));
    return Ball<F,FE>(rv,re);
}

template<class F, class FE> Ball<F,FE> Operations<Ball<F,FE>>::_mul(Ball<F,FE> const& x, Ball<F,FE> const& y) {
    auto rv=mul(near,x._v,y._v);
    auto ru=mul(up,x._v,y._v);
    auto rl=mul(down,x._v,y._v);
    auto re0=_make_error<FE>(hlf(sub(up,ru,rl)));
    auto re1=add(up,re0,mul(up,x._e,y._e));
    auto re2=add(up,mul(up,_make_error<FE>(abs(x._v)),y._e),mul(up,x._e,_make_error<FE>(abs(y._v))));
    auto re=add(up,re1,re2);
    return Ball<F,FE>(rv,re);
}

template<class F, class FE> Ball<F,FE> Operations<Ball<F,FE>>::_div(Ball<F,FE> const& x, Ball<F,FE> const& y) {
    return mul(x,rec(y));
}


} // namespace Ariadne
