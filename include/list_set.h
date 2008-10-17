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

namespace Ariadne {


template<class BS>
class ListSet 
{
 private:
  uint _dimension;
  std::vector<BS> _data;

 public:
  typedef typename std::vector<BS>::const_iterator const_iterator;

  virtual ~ListSet() { }

  ListSet() { };
  explicit ListSet(uint) { };
  explicit ListSet(const BS& bs) { this->adjoin(bs); }
  template<class BST> ListSet(const ListSet<BST>& ls) {
    this->_data.insert(this->end(),ls.begin(),ls.end()); }
  template<class Iter> ListSet(Iter first, Iter last) {
    this->_data.insert(this->end(),first,last); };



  /*! \brief Returns the number of basic sets forming this object. */
  uint size() const { return this->_data.size(); }

  /*! \brief Accesses the i-th BasicSet. */
  const BS& operator[](uint i) const { return this->_data[i]; };
  
  /*! \brief Make the set empty. */
  void clear() { this->_data.clear(); }
      
  /*! \brief A constant iterator to the beginning of the list of basic sets. */
  const_iterator begin() const { return this->_data.begin(); }
  
  /*! \brief A constant iterator to the end of the list of basic sets. */
  const_iterator end() const { return this->_data.end(); };

  /*! \brief Returns the denotable set's space dimension. */
  uint dimension() const { return this->_dimension; }

  /*! \brief Removes a set from the list and return it. */
  BS pop() { BS result=this->_data.back(); this->_data.pop_back(); return result; }

  /*! \brief Pushes a basic set to the end of the list. */
  void push_back(const BS& bs) { this->_data.push_back(bs); }

  /*! \brief Adjoins (makes union with) a basic set. */
  void adjoin(const BS& bs) { this->_data.push_back(bs); }

  /*! \brief Adjoins (makes union with) another list set. */
  void adjoin(const ListSet<BS>& ls) { this->_data.insert(this->_data.end(),ls.begin(),ls.end()); }
};      
  

} // namespace Ariadne

#endif /* ARIADNE_LIST_SET_H */
