/***************************************************************************
 *            Iterator.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file Iterator.hpp
 *  \brief Iterator support, similar to boost::Iterator package.
 */

#ifndef ARIADNE_ITERATOR_HPP
#define ARIADNE_ITERATOR_HPP

#include <iterator>

namespace Ariadne {

typedef OutputStream OutputStream;

struct RandomAccessTraversalTag { };
struct ForwardTraversalTag { };

class IteratorCoreAccess { };

template<class I, class Val, class Cat, class Ref=Val&> class IteratorFacade;

template<class I, class Val, class Ref> class IteratorFacade<I,Val,RandomAccessTraversalTag,Ref> {
    typedef Val* Ptr;
  public:
    // Standard typedefs
    typedef typename std::remove_const<Val>::type value_type;
    typedef Ref reference;
    typedef std::ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef Val* pointer;

    // Ariadne typedefs
    typedef typename std::remove_const<Val>::type ValueType;
    typedef Ref Reference;

    template<class II> bool operator==(const II& other) const { return static_cast<const I&>(*this).equal(other); }
    template<class II> bool operator!=(const II& other) const { return !static_cast<const I&>(*this).equal(other); }
    bool operator< (const I& other) const { return static_cast<const I&>(*this).distance_to(other)> 0; }
    bool operator<=(const I& other) const { return static_cast<const I&>(*this).distance_to(other)>=0; }
    bool operator> (const I& other) const { return static_cast<const I&>(*this).distance_to(other)< 0; }
    bool operator>=(const I& other) const { return static_cast<const I&>(*this).distance_to(other)<=0; }
    std::ptrdiff_t operator-(const I& other) const { return -static_cast<const I&>(*this).distance_to(other); }
    I operator+(std::ptrdiff_t n) const { I result(static_cast<const I&>(*this)); result.advance(n); return result; }
    I operator-(std::ptrdiff_t n) const { I result(static_cast<const I&>(*this)); result.advance(-n); return result; }
    I& operator+=(std::ptrdiff_t n) { static_cast<I&>(*this).advance(n); return static_cast<I&>(*this); }
    I& operator-=(std::ptrdiff_t n) { static_cast<I&>(*this).advance(-n); return static_cast<I&>(*this); }
    I& operator++() { static_cast<I&>(*this).advance(1); return static_cast<I&>(*this); }
    I& operator--() { static_cast<I&>(*this).advance(-1); return static_cast<I&>(*this); }
    Ref operator*() { return static_cast<I&>(*this).dereference(); }
    Ptr operator->() { return &static_cast<I&>(*this).dereference(); }
};

template<class I, class Val, class Ref> class IteratorFacade<I,Val,ForwardTraversalTag, Ref> {
    typedef Val* Ptr;
  public:
    // Standard typedefs
    typedef typename std::remove_const<Val>::type value_type;
    typedef Ref reference;
    typedef std::ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef Val* pointer;

    // Ariadne typedefs
    typedef typename std::remove_const<Val>::type ValueType;
    typedef Ref Reference;

    bool operator==(const I& other) const { return static_cast<const I&>(*this).equal(other); }
    bool operator!=(const I& other) const { return !static_cast<const I&>(*this).equal(other); }
    I& operator++() { static_cast<I&>(*this).increment(); return static_cast<I&>(*this); }
    I operator++(int) { I r=static_cast<const I&>(*this); static_cast<I&>(*this).increment(); return r; }
    Ref operator*() { return static_cast<I&>(*this).dereference(); }
    Ptr operator->() { return &static_cast<I&>(*this).dereference(); }
};

template<class I1, class I2> class PairIterator
//    : public IteratorFacade<PairIterator<I1,I2>, Monomial<X>, RandomAccessTag, MonomialReference<X>>
{
    typedef std::pair<I1*,I2*> Reference;
  public:
    I1 _iter1;
    I2 _iter2;
  public:
    typedef std::ptrdiff_t DifferenceType;

    PairIterator(I1 iter1, I2 iter2);
    template<class II1, class II2> PairIterator(const PairIterator<II1,II2>&);
  public:
    template<class II1, class II2> bool equal(const PairIterator<II1,II2>&) const;
    template<class II1, class II2> DifferenceType distance_to(const PairIterator<II1,II2>&) const;
    void advance(DifferenceType k);
    Reference dereference();
};
template<class I1, class I2> inline PairIterator<I1,I2>::PairIterator(I1 iter1, I2 iter2)
    : _iter1(iter1), _iter2(iter2) { }
template<class I1, class I2> template<class II1, class II2> inline PairIterator<I1,I2>::PairIterator(const PairIterator<II1,II2>& iter)
    : _iter1(iter._iter1), _iter2(iter._iter2) { }
template<class I1, class I2> template<class II1, class II2> inline bool PairIterator<I1,I2>::equal(const PairIterator<II1,II2>& other) const {
    return this->_iter1 == other._iter1 && this->_iter2 == other._iter2; }
template<class I1, class I2> inline void PairIterator<I1,I2>::advance(DifferenceType k) {
    _iter1+=(k), _iter2+=(k); }
template<class I1, class I2> inline auto PairIterator<I1,I2>::dereference() -> Reference {
    return std::make_pair(*this->_iter1,*this->_iter2); }
template<class I1, class I2> OutputStream& operator<<(OutputStream& os, const PairIterator<I1,I2>& e) {
    return os << "{" << e._iter1 << "," << e._iter2 << "}"; }

} // namespace Ariadne

#endif /* ARIADNE_ITERATOR_HPP */
