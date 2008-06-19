/***************************************************************************
 *            grid_set_iterator.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file grid_set_iterator.h
 *  \brief Iterators for grid sets. 
 */


#ifndef ARIADNE_GRID_SET_ITERATOR_H
#define ARIADNE_GRID_SET_ITERATOR_H

namespace Ariadne {
  

    template<class Base, class Value>
    class GridSetIterator 
      : public boost::iterator_adaptor<GridSetIterator<Base,Value>,Base,Value,boost::use_default,Value>
    { 
     public:
      typedef typename Value::real_type real_type;
      // Default constructor needed for some lists of lists. 
      GridSetIterator() 
        : GridSetIterator::iterator_adaptor_(), _grid() { }
      GridSetIterator(const Grid<real_type>& g, Base i) 
        : GridSetIterator::iterator_adaptor_(i), _grid(g) { }
     private:
      friend class boost::iterator_core_access;
      Value dereference() const { return Value(_grid,*this->base_reference()); }
      Grid<real_type> _grid;
    };

  
} // namespace Ariadne

#endif // ARIADNE_GRID_SET_ITERATOR_H
