/***************************************************************************
 *            partition_tree_set.tpl
 *
 *  1 July 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include "geometry/partition_tree_set.h"
#include "geometry/partition_tree_operations.h"

#include "geometry/rectangle.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"

#include <vector>

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    PartitionScheme<R>::PartitionScheme(const Rectangle<R>& bb) 
      : _bounding_box(bb), 
        _subdivisions(default_subdivision_coordinates(bb.dimension())) 
      { }

    template<typename R>
    PartitionTreeCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result.set_lower_bound( i,
          _bounding_box.lower_bound(i)+_cell.lower_bound(i)
            *(_bounding_box.upper_bound(i)-_bounding_box.lower_bound(i)) ); 
        result.set_upper_bound( i,
          _bounding_box.lower_bound(i)+_cell.upper_bound(i)
            *(_bounding_box.upper_bound(i)-_bounding_box.lower_bound(i)) ); 
      }
      return result;
    }

    template<typename R>
    PartitionTreeSet<R>::PartitionTreeSet(const GridMaskSet<R>& gms) 
      : _bounding_box(gms.bounding_box()),
        _unit_set(gms._unit_set)
    { }

    template<typename R>
    PartitionTreeSet<R>::operator ListSet<R,Rectangle>() const 
    {
      ListSet<R,Rectangle> res(this->dimension());
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(Rectangle<R>(*iter));
      }
      return res;
    }

/*
    template<typename R>
    PartitionTreeSet<R>::operator GridRectangleListSet<R>() const 
    {
      std::cerr << "PartitionTreeSet<R>::operator GridRectangleListSet<R>() const" << std::endl;
      FiniteGrid<R>* grid_ptr=new FiniteGrid<R>(this->bounding_box(),this->subdivisions());
      GridRectangleListSet<R> res(*grid_ptr);
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(GridRectangle<R>(res.grid(),Rectangle<R>(*iter)));
      }
      return res;
    }
*/

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionScheme<R>& g)
    {
      os << "PartitionScheme<" << name<R>() << ">(\n";
      os << "  bounding_box=" << g.bounding_box() << ",\n";
      os << "  subdivision_coordinates=" << g.subdivisions() << "\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTree<R>& t)
    {
      os << "PartitionTree<" << name<R>() << ">(\n";
      os << "  bounding_box=" << t.bounding_box() << ",\n";
      os << "  subdivisions=" << t.subdivisions() << "\n";
      os << "  words=";
      write_sequence(os, t.binary_tree().begin(), t.binary_tree().end());
      os << "\n";
      os << "  rectangles=[ " << Rectangle<R>(*t.begin());
      for(typename PartitionTree<R>::const_iterator ptree_iter=++t.begin() ; ptree_iter!=t.end(); ++ptree_iter) {
        os << ", " << Rectangle<R>(*ptree_iter);
      }
      os << " ]\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeCell<R>& c)
    {
      os << "PartitionTreeCell<" << name<R>() << ">(\n";
      os << "  bouding_box=" << c.bounding_box() << ",\n";
      os << "  subdivisions=" << c.subdivisions() << ",\n";
      os << "  word=" << c.word() << ",\n";
      os << "  bounds=" << c.bounds() << ",\n";
      os << "  rectangle=" << Rectangle<R>(c) << "\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeSet<R>& set)
    {
      os << "PartitionTreeSet<" << name<R>() << ">(\n";
      os << "  bounding_box=" << set.bounding_box() << ",\n";
      os << "  subdivisions=" << set.subdivisions() << ",\n";
      BinarySubtreeIterator bstb(set.binary_tree().begin(),set.mask().begin(),set.mask().end());
      BinarySubtreeIterator bste(set.binary_tree().end(),set.mask().end(),set.mask().end());
      os << "  words="; write_sequence(os, bstb, bste); os << ",\n";
      os << "  cells=["; 
      for(typename PartitionTreeSet<R>::const_iterator i=set.begin(); i!=set.end(); ++i) {
        os << (*i).bounds() << ",";
      }
      os << "]\n";
      os << ")\n";
      return os;
    }

  }
}
