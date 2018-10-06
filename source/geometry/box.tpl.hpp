/***************************************************************************
 *            geometry/box.tpl.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file geometry/box.tpl.hpp
 *  \brief
 */



#include "box.hpp"

namespace Ariadne {


template<class I> decltype(declval<I>().is_empty()) Box<I>::is_empty() const
{
    const Box<I>& bx=*this;
    decltype(declval<I>().is_empty()) res=false;
    for(SizeType i=0; i!=bx.size(); ++i) {
        res=res || bx[i].is_empty();
        if(definitely(res)) { return true; }
    }
    return res;
}


template<class I> decltype(declval<I>().is_bounded()) Box<I>::is_bounded() const
{
    const Box<I>& bx=*this;
    decltype(declval<I>().is_bounded()) res=true;
    for(SizeType i=0; i!=bx.size(); ++i) {
        res=res && bx[i].is_bounded();
    }
    return res;
}



template<class I> SizeType irmax(const Box<I>& bx) {
    SizeType imw(0);
    auto mw=bx[0].width();
    for(SizeType i=1; i!=bx.size(); ++i) {
        if(decide(bx[i].width()>mw)) { imw=i; mw=bx[i].width(); }
    }
    return imw;
}


template<class I> Box<I> Box<I>::split(SizeType k, SplitPart lmu) const
{
    const Box<I>& bx=*this;
    ARIADNE_ASSERT(k<bx.size());
    Box<I> r(bx);
    r[k]=Ariadne::split(bx[k],lmu);
    return r;
}

template<class I> Pair< Box<I>, Box<I> > Box<I>::split(SizeType k) const
{
    const Box<I>& bx=*this;
    ARIADNE_ASSERT(k<bx.size());
    Pair< Box<I>, Box<I> > r(bx,bx);
    auto c=bx[k].midpoint();
    r.first[k].set_upper(c);
    r.second[k].set_lower(c);
    return r;
}

template<class I> Box<I> Box<I>::split(SplitPart lmu) const
{
    return this->split(irmax(*this),lmu);
}

template<class I> Pair< Box<I>, Box<I> > Box<I>::split() const
{
    return this->split(irmax(*this));
}

template<class I> Box<I> Box<I>::bounding_box() const {
    return *this;
}

template<class I> typename Box<I>::MidpointType Box<I>::midpoint() const
{
    const Box<I>& bx=*this;
    MidpointType r(bx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        r[i]=bx[i].midpoint();
    }
    return r;
}

template<class I> typename Box<I>::CentreType Box<I>::centre() const {
    const Box<I>& bx=*this;
    CentreType r(bx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        r[i]=bx[i].centre();
    }
    return r;
}

template<class I> typename Box<I>::VertexType Box<I>::lower_bounds() const
{
    const Box<I>& bx=*this;
    VertexType r(bx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        r[i]=bx[i].lower();
    }
    return r;
}

template<class I> typename Box<I>::VertexType Box<I>::upper_bounds() const
{
    const Box<I>& bx=*this;
    VertexType r(bx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        r[i]=bx[i].upper();
    }
    return r;
}

template<class I> typename Box<I>::RadiusType Box<I>::radius() const
{
    const Box<I>& bx=*this;
    ARIADNE_ASSERT(bx.dimension()>0);
    decltype(declval<I>().radius()) r=bx[0].radius();
    for(SizeType i=1; i!=bx.size(); ++i) {
        r=max(r,bx[i].radius());
    }
    return r;
}

template<class I> typename Box<I>::RadiusType Box<I>::lengths() const
{
    const Box<I>& bx=*this;
    ARIADNE_ASSERT(bx.dimension()>0);
    decltype(declval<I>().width()) r=bx[0].width();
    for(SizeType i=1; i!=bx.dimension(); ++i) {
        r += bx[i].width();
    }
    return r;
}

template<class I> typename Box<I>::RadiusType Box<I>::measure() const
{
    const Box<I>& bx=*this;
    ARIADNE_ASSERT(bx.dimension()>0);
    decltype(declval<I>().width()) r=bx[0].width();
    for(SizeType i=1; i!=bx.dimension(); ++i) {
        r*=bx[i].width();
    }
    return r;
}


template<class I> Box<I> Box<I>::_project(const Box<I>& bx, const Array<SizeType>& rng)
{
    Box<I> res(rng.size());
    for(SizeType i=0; i!=rng.size(); ++i) {
        res[i]=bx[rng[i]];
    }
    return std::move(res);
}

template<class I> Box<I> Box<I>::_product(const Box<I>& bx1, const Box<I>& bx2)
{
    Box<I> r(bx1.size()+bx2.size());
    for(SizeType i=0; i!=bx1.size(); ++i) {
        r[i]=bx1[i];
    }
    for(SizeType i=0; i!=bx2.size(); ++i) {
        r[i+bx1.size()]=bx2[i];
    }
    return r;
}

template<class I> Box<I> Box<I>::_product(const Box<I>& bx1, const I& ivl2)
{
    Box<I> r(bx1.size()+1u);
    for(SizeType i=0; i!=bx1.size(); ++i) {
        r[i]=bx1[i];
    }
    r[bx1.size()]=ivl2;
    return r;
}

template<class I> Box<I> Box<I>::_product(const Box<I>& bx1, const Box<I>& bx2, const Box<I>& bx3)
{
    Box<I> r(bx1.size()+bx2.size()+bx3.size());
    for(SizeType i=0; i!=bx1.size(); ++i) {
        r[i]=bx1[i];
    }
    for(SizeType i=0; i!=bx2.size(); ++i) {
        r[i+bx1.size()]=bx2[i];
    }
    for(SizeType i=0; i!=bx3.size(); ++i) {
        r[i+bx1.size()+bx2.size()]=bx3[i];
    }
    return r;
}


}
