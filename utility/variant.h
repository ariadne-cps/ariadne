/***************************************************************************
 *            utility/variant.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file utility/variant.h
 *  \brief 
 */



#ifndef ARIADNE_VARIANT_H
#define ARIADNE_VARIANT_H

#include "metaprogramming.h"

namespace Ariadne {

//! \brief A class which can hold a value of one of two types, with fully-checked access.
template<class T1, class T2> class Variant {
    static const SizeType N=std::max(sizeof(T1),sizeof(T2));
    char _data[N];
    char _which;
    template<class T> struct Return { };
  public:
    ~Variant() { this->_destroy(); }
    Variant(const Variant<T1,T2>& other) { this->_create(other); }
    Variant<T1,T2>& operator=(const Variant<T1,T2>& other) { if(this!=&other) { this->_destroy(); this->_create(other); } return *this; }
    Variant(const T1& t1) { this->_create(t1); }
    Variant(const T2& t2) { this->_create(t2); }
    template<class T> T get() const { return this->_get(Return<T>()); }
    bool operator==(const Variant<T1,T2>& other) const { 
        if(this->_which==1) { return this->get<T1>()==other.get<T1>(); }
        if(this->_which==2) { return this->get<T2>()==other.get<T2>(); }
        if(this->_which!=other->_which) { return false; }
    }
    OutputStream& write(OutputStream& os) const { 
        if(this->_which==1) { return os << this->get<T1>(); }
        if(this->_which==2) { return os << this->get<T1>(); }
        return os;
    }
  private:
    T1 _get(Return<T1>) const { assert(_which==1); return *static_cast<const T1*>(_data); }
    T2 _get(Return<T2>) const { assert(_which==2); return *static_cast<const T2*>(_data); }
    Void _destroy() { 
        if(_which==1) { delete static_cast<T1*>(_data); } else if(_which==2) { delete static_cast<T1*>(_data); } _which=0; }
    Void _create(const Variant<T1,T2>& other) { 
        if(other._which==1) { this->_create(other.get<T1>()); } else if(other._which==2) { this->_create(other.get<T2>()); }  }
    Void _create(const T1& t1) { _which=1; new (_data) T1(t1); }
    Void _create(const T2& t2) { _which=2; new (_data) T2(t2); }
};
template<class T1, class T2> inline OutputStream& operator<<(OutputStream& os, Variant<T1,T2> const& t1t2) {
    return t1t2.write(os);
}
    
} // namespace Ariadne

#endif
