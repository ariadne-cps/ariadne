/***************************************************************************
 *            list_set.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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
 
/*! \file list_set.h
 *  \brief Sets which are lists of simple shapes.
 */

#ifndef ARIADNE_LIST_SET_H
#define ARIADNE_LIST_SET_H


#include <vector>
#include "stlio.h"
#include "macros.h"

#include "graphics_interface.h"
#include "discrete_location.h"

#include "box.h"

typedef unsigned int uint;

namespace Ariadne {

class Box;
template<class ES> class ListSet;
class DiscreteLocation;

// Declare template specialisation for hybrid list set
template<class ES> class HybridBasicSet;
template<class ES> class ListSet< HybridBasicSet<ES> >;

struct ListSetSummary { uint size, dimension; };

template<class BS>
class ListSet
    : public DrawableInterface
{
  private:
    std::vector<BS> _data;

  public:
    typedef typename std::vector<BS>::const_iterator const_iterator;
    typedef typename std::vector<BS>::iterator iterator;
    typedef BS value_type;

    virtual ~ListSet() { }

    ListSet() { };
    explicit ListSet(uint d) { };
    explicit ListSet(const BS& bs) { this->adjoin(bs); }
    template<class BST> ListSet(const ListSet<BST>& ls) {
        this->_data.insert(this->end(),ls.begin(),ls.end()); }
    ListSet(const std::vector<BS>& lst) {
        this->_data.insert(this->end(),lst.begin(),lst.end()); }
    template<class Iter> ListSet(Iter first, Iter last) {
        this->_data.insert(this->end(),first,last); };

    ListSet<BS>* clone() const { return new ListSet<BS>(*this); }

    /*! \brief Tests if the sets are equal. Tests equality of sequences, 
     *  including ordering. */   
    bool operator==(const ListSet<BS>& other) const { return this->_data==other._data; }


    /*! \brief Returns true if the list is empty. */
    bool empty() const { return this->_data.empty(); }

    /*! \brief Returns the number of basic sets forming this object. */
    size_t size() const { return this->_data.size(); }

    /*! \brief Accesses the i-th BasicSet. */
    const BS& operator[](size_t i) const { return this->_data[i]; };
  
    /*! \brief Make the set empty. */
    void clear() { this->_data.clear(); }
      
    /*! \brief A constant iterator to the beginning of the list of basic sets. */
    const_iterator begin() const { return this->_data.begin(); }
  
    /*! \brief A constant iterator to the end of the list of basic sets. */
    const_iterator end() const { return this->_data.end(); };

    /*! \brief A iterator to the beginning of the list of basic sets. */
    iterator begin() { return this->_data.begin(); }
  
    /*! \brief A iterator to the end of the list of basic sets. */
    iterator end() { return this->_data.end(); };

    /*! \brief Returns the denotable set's space dimension. */
    uint dimension() const { if(this->empty()) { return 0; } else { return this->_data.back().dimension(); } }

    /*! \brief Removes a set from the list and return it. */
    BS pop() { BS result=this->_data.back(); this->_data.pop_back(); return result; }

    /*! \brief Removes a set given identified by an iterator from the list. */
    iterator erase(iterator iter) { return this->_data.erase(iter); }

    /*! \brief Pushes a basic set to the end of the list. */
    void push_back(const BS& bs) { this->_data.push_back(bs); }

    /*! \brief Adjoins (makes union with) a basic set. */
    void adjoin(const BS& bs) { 
        ARIADNE_ASSERT(this->empty() || this->_data.back().dimension()==bs.dimension());
        this->_data.push_back(bs); }

    /*! \brief Adjoins (makes union with) another list set. */
    void adjoin(const ListSet<BS>& ls) { 
        ARIADNE_ASSERT(this->empty() || ls.empty() || this->_data.back().dimension()==ls._data.back().dimension());
        this->_data.insert(this->_data.end(),ls.begin(),ls.end()); }

    /*! \brief compute a list of the bounding boxes of the set elements. */
    ListSet<Box> bounding_boxes() const {
        ListSet<Box> result(this->dimension());
        for(uint i=0; i!=this->size(); ++i) { 
            result.adjoin((*this)[i].bounding_box());
        }
        return result;
    }

    /*! \brief A bounding box for the whole set. */
    Box bounding_box() const {
        if(this->size()==0) { return Box(this->dimension()); }
        Box result((*this)[0].bounding_box());
        for(uint i=1; i!=this->size(); ++i) {
            result=hull(result,(*this)[i].bounding_box()); }
        return result;
    }

    /*! \brief Draw on a canvas.
     *
     *  The ListSet template does not implement the DrawableInterface to avoid a dependency on the box.h header file.
     */
    void draw(CanvasInterface& c) const {
        for(uint i=0; i!=this->size(); ++i) {
            (*this)[i].draw(c);
        }
    }

    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const {
        return os << *this; 
    }

    /*! \brief Returns a summary of the size and dimension. */
    ListSetSummary summary() const { ListSetSummary summary = { this->size(),this->dimension() }; return summary; }

};
  
inline std::ostream& operator<<(std::ostream& os, const ListSetSummary& lss) {
    return os << "ListSet( s="<<lss.size<<", d="<<lss.dimension<<")"; }

template<class BS>
std::ostream& 
operator<<(std::ostream& os, const ListSet<BS>& ls)
{
    os << "ListSet";
    if(!ls.empty()) { for(uint i=0; i!=ls.size(); ++i) { os << (i==0?'(':',') << ls[i]; } }
    return os << ")";
}



} // namespace Ariadne

#endif /* ARIADNE_LIST_SET_H */
