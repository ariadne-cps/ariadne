/***************************************************************************
 *            base/binary_tree.h
 *
 *  Copyright  2007  Pieter Collins  pieter.collins@cwi.nl
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


#include <boost/iterator/indirect_iterator.hpp>



#ifndef ARIADNE_BASE_BINARY_TREE_H
#define ARIADNE_BASE_BINARY_TREE_H

namespace Ariadne {
  namespace Base {

    template<class T> class binary_tree;

    template<class T>
    class _binary_tree_node
    {
      typedef _binary_tree_node<T> _self;
      friend class binary_tree<T>;
     private:
      T _data;
      _self* _parent;
      _self* _children[2];
    };
      
    
    template<class T>
    class _binary_tree_iterator 
    {
     private:
      _binary_tree_node _curr;
    };


    /*! \brief A binary tree based on nodes. */
    template<class T> 
    class binary_tree
    {
      typedef _binary_tree_node<T> node_type;
      typedef node_type* node_pointer;
      typedef node_type& node_reference;
      typedef node_type* node_const_reference;
     public:
      binary_tree
      iterator set_children(iterator& p, const T& t1, const T& t2) {
        assert(p->_children[0]==0);
        assert(p->_children[1]==0);
        node_pointer child=new node_type(t1);
        p->_curr._children[0]=child;
        node_pointer child=new node_type(t2);
        p->_curr._children[1]=child;
      } 
      iterator set_child(iterator& p, bool lr, const T& t) {
        assert(p->_children[lr]==0);
        node_pointer child=new node_type(t);
        p->_children[lr]=child;
      } 
      void erase(iterator) {
        if(iterator._curr._parent._children[0]==iterator._curr) {
          node._parent._children[0]=0;
        } else {
          assert(iterator._curr._parent._children[1]==iterator._curr);
          iterator._curr._parent._children[1]=0;
        }
        this->_erase(iterator._curr);
      }                      
          
     private:
      void _erase(node_type& node) {
        node_type& node=iterator;
        if(node._children[0]) {
          this->_erase(node._children[0]);
        } 
        if(node._children[1]) {
          this->_erase(node._children[1]);
        }
        delete node;
      }                      
        
     private:
      _node_type _root;
    };
 

  }
}

#endif // ARIADNE_BASE_BINARY_TREE_H
